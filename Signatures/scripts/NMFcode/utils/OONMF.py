'''
class: NMFobject

functions:


__init__  -  initiate NMF instance with basic attributes

matrix_input_name  - set the filename if reading is needed

read_matrix_input - read the input matrices defined in matrix_input_name

performNMF - actually do the deed. Sets values for Basis and Mixture. will replace read matrix from previous step 

build_reconstruction - just take dot product of Basis and Mixture. Not done unless requested since this can take up a lot of memory. 

normalize_matrices - create NormedBasis and NormedMixture

compute_reweighted_matrices - computes ReweightedBasis and Reweighted Mixture. This is a specific reweighting method to attempt to attribute the elements of one matrix by understanding how much they contribute to the other. I.e. figure out how many DHSs are accounted for by the C1 in sample dimensions. 

normalize_reweighted_matrices - normalize the above matrices

writeNMF - write numpy binary files of Basis and Mixture. Mixture is not transposed in this case, preserving the NC x NDHS dimensionality

writeNMF_CSV - write CSV file for Basis and Mixture

define_colors - this sets the color scheme that we use for visualization

make_stacked_bar_plot - make our signature stacked bar plot. Should this really be part of the default library? I don't know but that's how I've decided to arrange things

make_anatomy_key - make quick visual showing the labels we associate with each NMF color in ENCODE project

make_standard_heatmap_plot - make a more traditional (matrix heatmap) visualization. 

precision_recall_curve - only works when objective matrix is known, and consists of entries 0/1. Compares reconstruction to the original data, using sort of precision/recall mechanics for samples. 

quick_precision_recall_curve - same as above, but only uses three threshold values - 0.3, 0.35, 0.4. Found to be the ideal choices.

precision_recall_curveDHS - uses the same method, but now computes precision/recall per DHS rather than per sample.

'''

ClusterMode = False
import sys
import numpy as np
import pandas as pd
if (ClusterMode):
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.decomposition import NMF
import OONMFhelpers as OH

today = OH.get_today()


class NMFobject:
    def __init__(self, theNcomps):
        self.Basis = []
        self.Mixture = []
        self.Ncomps = theNcomps
        
        self.BasisD = 0
        self.MixtureD = 0
        
        self.Basis_Names = []
        self.Mixture_Names = []
        
        self.Reconstruction = []
        
        self.ReweightedBasis = []
        self.ReweightedMixture = []
        
        self.NormedBasis = []
        self.NormedMixture = []

        self.ReweightedNormedBasis = []
        self.ReweightedNormedMixture = []

    
    def matrix_input_name(self, Basis_finname='', Mixture_finname=''):
        if (len(Basis_finname) < 1):
            print('syntax: read_matrix(Basis_finname, Mixture_finname)')
            sys.exit()
        self.Basis_finname = Basis_finname
        self.Mixture_finname = Mixture_finname

    def read_matrix_input(self, compressed = True):
        import gzip
        fn = self.Basis_finname
        if compressed: 
            fn = gzip.GzipFile(self.Basis_finname+'.gz', "r")

        self.Basis = np.load(fn)
        self.BasisD = self.Basis.shape[0]

        if (len(self.Mixture_finname)>0):
            fn = self.Mixture_finname
            if compressed: 
                fn = gzip.GzipFile(self.Mixture_finname+'.gz', "r")

            self.Mixture = np.load(fn)
            self.MixtureD = self.Mixture.shape[1]
        
        
        
    def performNMF(self, data, randomseed=0, theinit='random', thesolver='cd', thebetaloss='frobenius'):
        if(len(self.Basis) > 0):
            print('you are overwriting the Basis',self.Basis)
            cont = input('are you sure?')
            if (cont == 'n'):
                return
        print("max iter set to 10000 here")
        model = NMF(n_components=self.Ncomps, init=theinit, random_state=randomseed, solver=thesolver, beta_loss=thebetaloss,max_iter=10000)
        print('starting NMF at ', OH.mytime(), flush=True)
        self.Basis = model.fit_transform(data) 
        print('done with NMF at ', OH.mytime(), flush=True)
        self.Mixture = model.components_
        self.BasisD = self.Basis.shape[0]
        self.MixtureD = self.Mixture.shape[1]
        print('returning reconstruction error')
        return model.reconstruction_err_

    def performNMF_KL(self, data, randomseed=0):
        if(len(self.Basis) > 0):
            print('you are overwriting the Basis',self.Basis)
            cont = input('are you sure?')
            if (cont == 'n'):
                return
        model = NMF(n_components=self.Ncomps, init='random', random_state=randomseed, solver='mu', beta_loss ='kullback-leibler')
        print('starting NMF at ', OH.mytime(), flush=True)
        self.Basis = model.fit_transform(data) 
        print('done with NMF at ', OH.mytime(), flush=True)
        self.Mixture = model.components_
        self.BasisD = self.Basis.shape[0]
        self.MixtureD = self.Mixture.shape[1]

    def performNMF_MU(self, data, randomseed=0):
        if(len(self.Basis) > 0):
            print('you are overwriting the Basis',self.Basis)
            cont = input('are you sure?')
            if (cont == 'n'):
                return
            
        model = NMF(n_components=self.Ncomps, init='random', random_state=randomseed, solver='mu', beta_loss ='frobenius')
        print('starting NMF at ', OH.mytime(), flush=True)
        self.Basis = model.fit_transform(data) 
        print('done with NMF at ', OH.mytime(), flush=True)
        self.Mixture = model.components_
        self.BasisD = self.Basis.shape[0]
        self.MixtureD = self.Mixture.shape[1]


    def build_reconstruction(self):
        self.Reconstruction = np.dot(self.Basis, self.Mixture)
                

    def normalize_matrices(self):
        print(self.Mixture.shape)
        m = self.Mixture[0:self.Ncomps,np.sum(self.Mixture, axis=0)>0]
        print("filtering out rows with 0 sum")
        print(m.shape)
        self.NormedMixture = m / np.sum(m, axis=0)
        #self.NormedMixture =   self.Mixture / np.sum(self.Mixture, axis=0)
        self.NormedBasis =   (self.Basis.T / np.sum(self.Basis.T, axis=0)).T

    def compute_reweighted_matrices(self):
    
        bigAllDHSSum_ar = []
        bigAllSampleSum_ar = []
        
        
        for i in range(self.Ncomps):
            bongo = np.copy(self.Basis)
            for k in range(self.Ncomps):
                if (k!=i):
                    bongo[:,k]*=0
            sansvar = np.dot(bongo, self.Mixture)
            bigAllDHSSum_ar.append(np.sum(sansvar[:,0:], axis=1))
            bigAllSampleSum_ar.append(np.sum(sansvar[:,0:], axis=0))
            del(sansvar)
        
        self.ReweightedBasis = np.array(bigAllDHSSum_ar).T
        self.ReweightedMixture = np.array(bigAllSampleSum_ar)
            

    def normalize_reweighted_matrices(self):
        self.ReweightedNormedMixture =   self.ReweightedMixture / np.sum(self.ReweightedMixture, axis=0)
        self.ReweightedNormedBasis =   (self.ReweightedBasis.T / np.sum(self.ReweightedBasis.T, axis=0)).T


    def writeNMF(self, Basis_foutname, Mixture_foutname, compressed=True):
        import gzip
        f = gzip.GzipFile(Basis_foutname+'.gz', "w")
        np.save(f, self.Basis)
        f.close()

        #very confusing but it must be Mixture here for internal self-consistency. Can be Mixture.T for CSV files
        f = gzip.GzipFile(Mixture_foutname+'.gz', "w")
        np.save(f, self.Mixture)
        f.close()
        
        
    def writeNMF_CSV(self, Basis_foutname, Mixture_foutname, compressed=True):
        pd.DataFrame(self.Basis).to_csv(Basis_foutname+'.gz', compression='infer')
        pd.DataFrame(self.Mixture.T).to_csv(Mixture_foutname+'.gz', compression='infer')


    def define_colors(self, reordercolors=False):

        maxassigned = 16
        self.Comp_colors = ['#FFE500', '#FE8102', '#FF0000', '#07AF00', '#4C7D14', '#414613', '#05C1D9', '#0467FD', '#009588', '#BB2DD4', '#7A00FF', '#4A6876', '#08245B', '#B9461D', '#692108', '#C3C3C3']
        neworder = np.array([16,10,7,11,2,12,1,8,4,15,14,5,9,6,3,13]).astype(int) - 1
        
        self.Comp_colors = list(np.array(self.Comp_colors)[neworder])
        
        if (self.Ncomps>maxassigned):
            # somewhat defunct but whatever. Adds extra "random" colors if you use more than 16 
            from matplotlib import colors as mcolors
            colornames = np.sort(list(mcolors.CSS4_COLORS.keys()))
            count = maxassigned
            np.random.seed(10)
            myrandint = np.random.randint(len(colornames))
            while (count < self.Ncomps):
                myrandint =    np.random.randint(len(colornames))
                newcolor = colornames[myrandint]
                trialcount = 0
                while ((newcolor in self.Comp_colors) and (trialcount < 100)):
                    print('what am i doing here')
                    newcolor = colornames[np.random.randint(0,len(colornames))]
                    trialcount+=1
                print('new color ',count,newcolor)
                self.Comp_colors.append(newcolor)
                count+=1



    def make_stacked_bar_plot(self, Nrelevant, BarMatrix, bargraph_out, names = [], plotClusterMode=False, barsortorder=[], clusterTopLabels=[], plot_title='', official_order = False, no_axis=False, figdim1 = 150, figdim2 = 40):
    
        # define barsortorder if one isn't provided 
        if len(barsortorder)<1:
            barsortorder = np.arange(Nrelevant)
        
        #define names if none are provided
        if len(names) < 1:
            names = [str(i) for i in range(Nrelevant)]
            names = np.array(names)
            
        #Make a set of x coordinates for ticks
        Xpositions = np.arange(Nrelevant)
        
        # start and end matrices for each matrix. This is to ensure that you can plot only Nrelevant vectors from the matrix if that is what you want
        start = 0
        end = Nrelevant
        
        self.define_colors()

        if official_order:
          #  WSO = np.array([7,5,15,9,12,14,3,8,13,2,4,6,16,11,10,1]).astype(int) - 1
            WSO = np.array([8, 19, 16, 13, 14, 10, 18,  9,  2,  1, 17,  5,  4,  3, 15, 12,  7, 0,  6, 11]).astype(int) - 1
            BarMatrix = BarMatrix[WSO]
        else:
            print("else WSO")
            WSO = np.arange(self.Ncomps)

        
        #this is really the only meaty part
        plt.clf()
        plt.figure(figsize=(figdim1,figdim2))
        
        ground_pSample = np.zeros(len(Xpositions))
        print(self.Ncomps,start, end,WSO,"tellings comps")
        print(names[1:10])
        #nms=names[barsortorder]
        nms = [x for _,x in sorted(zip(barsortorder,names))]
        with open ("rank"+str(self.Ncomps)+".orderd.matrix","w") as FH:
             FH.write("\t".join(nms))
             FH.write("\n")
             for i in range(self.Ncomps):
                 FH.write("\t".join([str(x) for x in BarMatrix[i,start:end][barsortorder]]))
                 FH.write("\n")
        print("File write done")   
        plt.bar(Xpositions[start:end],BarMatrix[i,start:end][barsortorder], width=1.0,bottom = ground_pSample, color=self.Comp_colors[WSO[i]], alpha=1)
        ground_pSample = np.sum(BarMatrix[0: i+1,start:end], axis=0)[barsortorder]
        OH.increase_axis_fontsize()
        
        #axis labels - seems highly optional
        plt.ylabel('sum of signal in matrix',fontsize=70)
        if (len(plot_title) > 0):
	        plt.title(plot_title)
	        
	    #heuristic scaling of bottom
        samplenamesize = (1/Nrelevant)**0.5 * 300
        thebottom = min([(1/Nrelevant)**0.3 * 1.2, 0.3])
        
        #i think this is largely defunct, but i guess this can make some extra labels on the top of the plot
        if(plotClusterMode):
            plt.xticks(Xpositions, Xpositions.astype(str), rotation='vertical', fontsize=samplenamesize)
            if len(clusterTopLabels) > 0:
                ax = plt.gca()
                ax2 = ax.twiny()
                ax2.set_xticks(Xpositions)
                ax2.set_xticklabels(clusterTopLabels.astype(str), rotation=90, fontsize=samplenamesize)
        
        # default behavior 
        else:
            plt.xticks(Xpositions, names[barsortorder], rotation='vertical', fontsize=samplenamesize)	
            
        #adjust it so that it fits in the fame
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=thebottom)
        if (no_axis):
            plt.axis('off')
        
        plt.savefig(bargraph_out)
        plt.show()
        plt.close()
        
    def make_anatomy_key(self, legendout=''):
        self.define_colors()
        strings_of_labels = ['Tissue invariant', 'Stromal A', 'Embryonic / primitive', 'Stromal B', 'Lymphoid', 'Renal / cancer', 'Placental','Neural','Cardiac','Organ devel. / renal','Pulmonary devel.','Musculoskeletal',\
                     'Digestive','Vascular / endothelial','Myeloid / erythroid', 'Cancer / epithelial']
        strings_of_labels = np.array(strings_of_labels)
        fig, ax = plt.subplots(1, 2, figsize = (8,4))
        ax[1].set_xlim([-5, 5])
        for i in range(8):
            ax[1].hlines(i+1, 0, 1, color=self.Comp_colors[15-i], lw=20)
            ax[1].text(1.2, i+1-0.2, strings_of_labels[15-i],fontsize=15)
        ax[1].axis('off')
        ax[0].set_xlim([-5, 5])

        for i in range(8, 16):
            ax[0].hlines(i+1, 0, 1, color=self.Comp_colors[15-i], lw=20)
            ax[0].text(1.2, i+1-0.2, strings_of_labels[15-i],fontsize=15)
        ax[0].axis('off')
        plt.savefig(legendout)
        plt.show()
        plt.close()
    

    def make_standard_heatmap_plot(self, Nrelevant, BarMatrix, bargraph_out, names = [], plotClusterMode=False, barsortorder=[], clusterTopLabels=[], plot_title='', official_order = False, no_axis=False, figdim1 = 150, figdim2 = 40):
        plt.figure(figsize=(figdim1,figdim2))
        if len(barsortorder)<1:
            print("sorting")
            barsortorder = np.arange(Nrelevant)
        
        #define names if none are provided
        if len(names) < 1:
            names = [str(i) for i in range(Nrelevant)]
            names = np.array(names)
        print("saving ordered matrix")
        df = pd.DataFrame(data=BarMatrix[barsortorder].T)
        df.columns = names.astype(str)[barsortorder]
        df.to_csv(bargraph_out+'.order.matrix', sep='\t', header=True, float_format='%.2f', index=False)
        plt.imshow(BarMatrix[barsortorder].T, cmap='Blues', aspect='auto')
        plt.xticks(np.arange(len(barsortorder)), names.astype(str)[barsortorder], rotation=90, fontsize= (1/Nrelevant)**0.5 * 300 )
        plt.savefig(bargraph_out, bbox_inches='tight', transparent=True)
        plt.show()	
        
        
    def precision_recall_curve(self, data, names=[], write_verbose=False, filename_addon='', verbose=True):
        #only works when objective matrix is known, and consists of entries 0/1. 
        if (len(self.Reconstruction)<1):
            self.build_reconstruction()
        if verbose:
            print(data.shape, 'data')
            print(self.Reconstruction.shape, 'reconstruction')
        if (data.shape[0] != self.Reconstruction.shape[0] or data.shape[1] != self.Reconstruction.shape[1]):
            print('error! data and reconstruciton dont have matching shapes', data.shape, self.Reconstruction.shape )
            return
        if (np.max(data) > 1 or np.min(data) < 0):
            print('error, precision-recall curve only works for data between 0 and 1')
            return

        customthreshes = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,0.95]
        recall_ar = []
        precision_ar = []

        if len(names) < 1:
            if verbose:
                print('filling in names')
            names = np.arange(self.BasisD).astype(str)
            
        sample_based_table = []
        total_PR_talbe = []
        for customthresh in customthreshes:
            F1_ar = []
            if(write_verbose):
                f = open(filename_addon+'SampleCSthresh'+str(customthresh)+'.txt', 'w')
            count = 0

            totalTP = 0
            totalTN = 0
            totalFP = 0
            totalFN = 0
            for sample in range(self.BasisD):
                DHSar_cut = data[sample]>customthresh
                predDHSar_cut = self.Reconstruction[sample] > customthresh


                TP = len(self.Reconstruction[sample][DHSar_cut * predDHSar_cut]) 
                FP = len(self.Reconstruction[sample][predDHSar_cut * np.invert( DHSar_cut)])
                TN = len(self.Reconstruction[sample][np.invert(predDHSar_cut) * np.invert(DHSar_cut)])
                FN = len(self.Reconstruction[sample][np.invert(predDHSar_cut) * (DHSar_cut)])

                if ((TP + FN) > 0 ):
                    recall = TP / (TP + FN)
                else: 
                    recall = 0
                if ((TP + FP) > 0 ):
                    precision = TP / (TP + FP)
                else:
                    precision=0
                accuracy = (TP + TN) /(len(self.Reconstruction[sample]))
                if (precision + recall) == 0:
                    F1 = 0
                else:
                    F1 = 2*(precision*recall)/(precision+recall)
                F1_ar.append(F1)
                totalTP += TP
                totalTN += TN
                totalFN += FN
                totalFP += FP
                
                if (write_verbose):
                    print(sample, names[count], TP, FP, TN, FN, recall, precision, accuracy, F1, file=f)
                sample_based_table.append([customthresh, sample, names[count], TP, FP, TN, FN, recall, precision, accuracy, F1])
                
                count +=1
            if verbose:
                print('Ncomps ',self.Ncomps, 'thresh ', customthresh, ' mean F1 score ',np.mean(np.array(F1_ar)))
            total_PR_talbe.append([customthresh, totalTP, totalFP, totalTN, totalFN])
        
        total_PR_tableDF = pd.DataFrame(total_PR_talbe, columns=['threshold', 'TP', 'FP', 'TN', 'FN'])
        total_PR_tableDF.to_csv(filename_addon+'TotalPR'+'.txt', sep='\t', index=False)
        
        sample_based_tableDF = pd.DataFrame(sample_based_table, columns=['threshold', 'sample_number', 'sample_name', 'TP', 'FP', 'TN', 'FN', 'recall', 'precision', 'accuracy', 'F1'])
        return [sample_based_tableDF,total_PR_tableDF]
        

    def quick_precision_recall_curve(self, data, names=[], write_verbose=False, filename_addon=''):
        if (len(self.Reconstruction)<1):
            self.build_reconstruction()
        print(data.shape, 'data')
        print(self.Reconstruction.shape, 'reconstruction')
        if (data.shape[0] != self.Reconstruction.shape[0] or data.shape[1] != self.Reconstruction.shape[1]):
            print('error! data and reconstruciton dont have matching shapes', data.shape, self.Reconstruction.shape )
            return
        if (np.max(data) > 1 or np.min(data) < 0):
            print('error, precision-recall curve only works for data between 0 and 1')
            return

        customthreshes = [0.3, 0.35, 0.4]
        recall_ar = []
        precision_ar = []

        if len(names) < 1:
            print('filling in names')
            names = np.arange(self.BasisD).astype(str)
        sample_based_table = []
        for customthresh in customthreshes:
            F1_ar = []
            if(write_verbose):
                f = open(filename_addon+'SampleCSthresh'+str(customthresh)+'.txt', 'w')

            count = 0

            for sample in range(self.BasisD):
                DHSar_cut = data[sample]>customthresh
                predDHSar_cut = self.Reconstruction[sample] > customthresh


                TP = len(self.Reconstruction[sample][DHSar_cut * predDHSar_cut]) 
                FP = len(self.Reconstruction[sample][predDHSar_cut * np.invert( DHSar_cut)])
                TN = len(self.Reconstruction[sample][np.invert(predDHSar_cut) * np.invert(DHSar_cut)])
                FN = len(self.Reconstruction[sample][np.invert(predDHSar_cut) * (DHSar_cut)])

                if ((TP + FN) > 0 ):
                    recall = TP / (TP + FN)
                else: 
                    recall = 0
                if ((TP + FP) > 0 ):
                    precision = TP / (TP + FP)
                else:
                    precision=0
                accuracy = (TP + TN) /(len(self.Reconstruction[sample]))
                if (precision + recall) == 0:
                    F1 = 0
                else:
                    F1 = 2*(precision*recall)/(precision+recall)
                F1_ar.append(F1)
                if (write_verbose):
                    print(sample, names[count].strip(' '), TP, FP, TN, FN, recall, precision, accuracy, F1, file=f)
                sample_based_table.append([customthresh, sample, names[count].strip(' '), TP, FP, TN, FN, recall, precision, accuracy, F1])
                count +=1
            print('Ncomps ',self.Ncomps, 'thresh ', customthresh, ' mean F1 score ',np.mean(np.array(F1_ar)))
        return pd.DataFrame(sample_based_table, columns=['threshold', 'sample_number', 'sample_name', 'TP', 'FP', 'TN', 'FN', 'recall', 'precision', 'accuracy', 'F1'])


    def precision_recall_curveDHS(self, data, names=[], write_verbose=False, filename_addon=''):
        #only works when objective matrix is known, and consists of entries 0/1. 
        if (len(self.Reconstruction)<1):
            self.build_reconstruction()
        print(data.shape, 'data')
        print(self.Reconstruction.shape, 'reconstruction')
        if (data.shape[0] != self.Reconstruction.shape[0] or data.shape[1] != self.Reconstruction.shape[1]):
            print('error! data and reconstruciton dont have matching shapes', data.shape, self.Reconstruction.shape )
            return
        if (np.max(data) > 1 or np.min(data) < 0):
            print('error, precision-recall curve only works for data between 0 and 1')
            return

        customthreshes = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35]
        recall_ar = []
        precision_ar = []

        if len(names) < 1:
            print('filling in names')
            names = np.arange(self.MixtureD).astype(str)
        sample_based_table = []
        for customthresh in customthreshes:
            F1_ar = []
            if(write_verbose):
                f = open(filename_addon+'DHSCSthresh'+str(customthresh)+'.txt', 'w')

            count = 0

            for DHS in range(self.MixtureD):
                Sample_ar_cut = data[:,DHS]>customthresh
                predSamplear_cut = self.Reconstruction[:,DHS] > customthresh


                TP = len(self.Reconstruction[:,DHS][Sample_ar_cut * predSamplear_cut]) 
                FP = len(self.Reconstruction[:,DHS][predSamplear_cut * np.invert( Sample_ar_cut)])
                TN = len(self.Reconstruction[:,DHS][np.invert(predSamplear_cut) * np.invert(Sample_ar_cut)])
                FN = len(self.Reconstruction[:,DHS][np.invert(predSamplear_cut) * (Sample_ar_cut)])

                if ((TP + FN) > 0 ):
                    recall = TP / (TP + FN)
                else: 
                    recall = 0
                if ((TP + FP) > 0 ):
                    precision = TP / (TP + FP)
                else:
                    precision=0
                accuracy = (TP + TN) /(len(self.Reconstruction[:,DHS]))
                if (precision + recall) == 0:
                    F1 = 0
                else:
                    F1 = 2*(precision*recall)/(precision+recall)
                F1_ar.append(F1)
                if (write_verbose):
                    print(DHS, names[count], TP, FP, TN, FN, recall, precision, accuracy, F1, file=f)
                sample_based_table.append([customthresh, DHS, names[count], TP, FP, TN, FN, recall, precision, accuracy, F1])
                count +=1
            print('Ncomps ',self.Ncomps, 'thresh ', customthresh, ' mean F1 score ',np.mean(np.array(F1_ar)))
        return pd.DataFrame(sample_based_table, columns=['threshold', 'DHS_number', 'DHS_name', 'TP', 'FP', 'TN', 'FN', 'recall', 'precision', 'accuracy', 'F1'])

        
