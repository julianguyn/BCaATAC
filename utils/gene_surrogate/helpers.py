import pandas as pd
import numpy as np
from sklearn.model_selection import KFold, GridSearchCV
from sklearn import linear_model
from sklearn.metrics import root_mean_squared_error, mean_absolute_error, make_scorer
from scipy.stats import spearmanr, pearsonr
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
newcmp = mcolors.LinearSegmentedColormap.from_list("custom", ["blue", "red"])

def important_features(train, dt, thres):

  # initialize dictionary to hold correlation results
  corr_dict = {}

  # correlate exp of each gene to drug response
  for feature in train.columns:
    corr_dict[feature] = train[feature].corr(dt)
  correlations = pd.DataFrame.from_dict(corr_dict, orient='index', columns=['Correlation'])

  # count number of univariable associations that meet the threshold
  num_pos = (correlations['Correlation'] > thres).sum()
  num_neg = (correlations['Correlation'] < -thres).sum()

  print('Selected threshold:', thres)
  print('Num features with positive correlation > threshold:', num_pos)
  print('Num features with negative correlation > threshold:', num_neg)
  print('Total number of features:', str(num_pos + num_neg))

  # identify features that pass selected threshold
  keep = correlations[correlations['Correlation'].abs() > thres].index

  # subset training dataframe to only genes of interest
  X_subset = train[keep]

  return X_subset

def corr_features(X_subset, thres):

    # correlate exp of remaining genes
    corr_mat = X_subset.corr(method='pearson', min_periods=1)
    vals = corr_mat.values.copy()
    np.fill_diagonal(vals, 0)
    corr_mat.iloc[:, :] = vals  

    corr_pairs = (corr_mat.abs() > thres)
    correlated = set() # initialize set to store correlated features

    # loop through correlated pairs
    for i in range(corr_pairs.shape[0]):
        for j in range(i + 1, corr_pairs.shape[1]):

            # if True (highly correlated)
            if corr_pairs.iloc[i, j]:

                #print(corr_mat.columns[i], corr_mat.columns[j])

                # add one of the genes to the set
                correlated.add(corr_mat.columns[i])

    #print('Highly correlated features:', correlated)
    print('\nNum correlated features:', len(correlated))

    print('Original number of feaures:', X_subset.shape[1])

    # remove correlated genes
    X_subset = X_subset.drop(columns=list(correlated))
    print('Number of features remaining:', X_subset.shape[1])

    return X_subset

def spearman_scorer(y_true, y_pred):
    """
    Helper function to create scorer
    """
    return spearmanr(y_true, y_pred).correlation


def run_elastic_net(X, y, path, full_cohort, repeats=10, preds=True, features=True):

    # initialize variables to store results
    fold_preds = []
    results = []
    hyperparams_list = []
    feature_counter = Counter()
    total_models = 0

    # format full cohort
    full_cohort = full_cohort.loc[:, full_cohort.columns.isin(X.columns)]
    full_cohort = full_cohort.reindex(columns=X.columns)

    # repeats of 5foldCV
    for repeat in range(repeats):

        # initialize outer folds (5 folds, 80% train, 20% test)
        outer_cv = KFold(n_splits=5, shuffle=True, random_state=101 + repeat)
        fold = 1

        # loop through each of the outer five folds
        for train_index, test_index in outer_cv.split(X):

            # initialize inner folds (5 folds, 80% train, 20% test)
            inner_cv = KFold(n_splits=5, shuffle=True, random_state=fold)

            # split train and test
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]

            # initialize ElasticNet model
            en = linear_model.ElasticNet()

            # specify parameters for optimization
            parameters = {
                'alpha': [0.1, 1, 10, 100],
                'l1_ratio': [0.2, 0.5, 0.8],
                'max_iter': [1000, 5000, 7500]
            }

            # identify optimal parameters
            reg = GridSearchCV(
                estimator = en,
                param_grid = parameters,
                cv=inner_cv,
                scoring=make_scorer(spearman_scorer, greater_is_better=True),
                n_jobs=-1
            )

            # fit model
            reg.fit(X_train, y_train)

            # get best model & parameters
            reg_best = reg.best_estimator_

            best_params = reg.best_params_
            hyperparams_list.append(best_params)

            # get selected features
            selected_features = X.columns[reg_best.coef_ != 0]
            feature_counter.update(selected_features)
            total_models += 1

            # get predicted values for test data
            y_pred = pd.Series(reg_best.predict(X_test), index=y_test.index)
            fold_preds.append(pd.DataFrame({
                "sample": y_test.index,
                "y_true": y_test.values,
                "y_pred": y_pred.values,
                "repeat": repeat+1,
                "fold": fold
            }))

            # compute metrics
            s_corr = spearmanr(y_test, y_pred).correlation
            p_corr = pearsonr(y_test, y_pred)[0]
            rmse = root_mean_squared_error(y_test, y_pred)
            mae = mean_absolute_error(y_test, y_pred)

            # save model correlation and features
            results.append({
                "repeat": repeat+1,
                "fold": fold,
                "spearman": s_corr,
                "pearson": p_corr,
                "rmse": rmse,
                "mae": mae,
                **best_params
            })

            fold += 1

    # save results to dataframe
    results_df = pd.DataFrame(results)
    results_df.to_csv(("data/results/data/"+path+"en_folds.csv"), index=False)

    # get most common hyperparameters
    params_df = pd.DataFrame(hyperparams_list)
    final_params = params_df.mode().iloc[0].to_dict()
    final_params['max_iter'] = int(final_params['max_iter'])

    # fit final model
    final_model = linear_model.ElasticNet(**final_params)
    final_model.fit(X, y)

    y_pred = pd.Series(final_model.predict(full_cohort), index=full_cohort.index)
    full_preds = pd.DataFrame({
        "sample": full_cohort.index,
        "y_pred": y_pred.values
    })
    full_preds.to_csv(("data/results/data/"+path+"full_cohort_predictions.csv"), index=False)

    if features:
        feature_freq = pd.DataFrame.from_dict(feature_counter, orient='index', columns=['count'])
        feature_freq['frequency'] = feature_freq['count'] / total_models
        feature_freq = feature_freq.sort_values(by='frequency', ascending=False)
        feature_freq.to_csv(f"data/results/data/{path}en_feature_stability.csv")

    if preds:
        fold_preds_df = pd.concat(fold_preds)
        fold_preds_df.to_csv(f"data/results/data/{path}en_predictions.csv")