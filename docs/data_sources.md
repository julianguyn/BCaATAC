# Data Sources

## Raw Data

Includes all data from `data/rawdata/`

### data/rawdata/tcga

```markdown
## BCa_binary.2.matrix (Internal)

- **Name**: Binary peak matrix from 75 BCa TCGA tumours (input into NMF)
- **Creation Date**: 2023-10-12
- **Creation Method**: Lupien Lab binary matrix pipeline
- **Input Data**: TCGA BCa tumour summit and narrowpeak files
- **Processing Scripts**: todo
```

```markdown
## Human__TCGA_BRCA__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct

- **Name**: TCGA BCa Tumours Gene Counts Matrix
- **Version**: Dataset created 2016-01-28
- **URL**: https://linkedomics.org/data_download/TCGA-BRCA/
- **Access Method**: Download from url
```

```markdown
## Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi

- **Name**: TCGA BCa Tumours metadata
- **Version**: Dataset created 2016-01-28
- **URL**: https://linkedomics.org/data_download/TCGA-BRCA/
- **Access Method**: Download from url
```

```markdown
## TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx

- **Name**: TCGA BCa Tumours sample ID mapping file
- **Version**: 
- **URL**: https://portal.gdc.cancer.gov/v1/annotations
- **Access Method**: Download from url
```

```markdown
## TCGA_sourcefiles.csv

- **Name**: TCGA BCa Tumours sources of files 
- **Notes**: Manually mapped ATAC-Seq and SNV filenames to BCa TCGA samples
```

```markdown
## /maf and /idat

- **Name**: .maf and .idat files
- **Version**: 
- **URL**: https://portal.gdc.cancer.gov/projects/TCGA-BRCA
```

### data/rawdata/ccls

```markdown
## bca_sign.Zscore.txt (Internal)

- **Name**: BCa cell line ARCHE scores
- **Creation Date**: 
- **Creation Method**: Lupien Lab consensus peaks and scoring pipeline
- **Input Data**: BCa cell line summit and narrowpeak files
- **Processing Scripts**: [Github](https://github.com/julianguyn/bca_arche_scoring)
```

```markdown
## bcacell_lines (Internal) - not used anymore

- **Name**: Binary peak matrix from BCa cell lines
- **Creation Date**: 
- **Creation Method**: Lupien Lab binary matrix pipeline
- **Input Data**: BCa cell line summit and narrowpeak files
- **Processing Scripts**: todo
```

### data/rawdata/psets

```markdown
## CCLE.rds

- **Name**: CCLE PSet
- **Version**: Dataset created 2020-06-24
- **Dataset DOI**: 10.5281/zenodo.3905461
- **URL**: https://orcestra.ca/pset/5ef3659e785cc307d861e79c
- **Access Method**: Download from url
- **Citation**: See url
- **License**: [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/)
```

```markdown
## CTRP.rds

- **Name**: CTRPv2 PSet
- **Version**: Dataset created 2020-06-24
- **Dataset DOI**: 10.5281/zenodo.7826870
- **URL**: https://orcestra.ca/pset/5ef3659e785cc307d861e79b
- **Access Method**: Download from url
- **Citation**: See url
- **License**: [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/)
```

```markdown
## gCSI.rds

- **Name**: gCSI PSet
- **Version**: Dataset created 2021-06-11
- **Dataset DOI**: 10.5281/zenodo.7829857
- **URL**: https://orcestra.ca/pset/60c3dc783940cf1de1bbc298
- **Access Method**: Download from url
- **Citation**: See url
- **License**: [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/)
```

```markdown
## GDSC2-8.2.rds

- **Name**: GDSC_2020(v2-8.2) PSet
- **Version**: Dataset created 2021-12-16
- **Dataset DOI**: 10.5281/zenodo.5787145
- **URL**: https://orcestra.ca/pset/61bb751a308ac5003a648fbe
- **Access Method**: Download from url
- **Citation**: See url
- **License**: 
```

```markdown
## PharmacoSet.RDS (Internal)

- **Name**: BHKLab Unreleased UHNBreast2 PSet
- **Creation Date**: 2024-03-18
- **Creation Method**: RNA-Seq processed using Kallisto
- **Input Data**: Drug response data from archive/PharmacoGxCards, RNA-Seq data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73526)
- **Processing Scripts**: todo
```

```markdown
## PSet_GRAY2017.rds

- **Name**: GRAY_2017
- **Version**: Dataset created 2021-02-23
- **Dataset DOI**: 10.5281/zenodo.7826847
- **URL**: https://orcestra.ca/pset/5ef3659e785cc307d861e795
- **Access Method**: Download from url
- **Citation**: See url
- **License**: [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/)
```

```markdown
## PSet_GRAY2017.rds

- **Name**: GRAY_2017
- **Version**: Dataset created 2021-02-23
- **Dataset DOI**: 10.5281/zenodo.7826847
- **URL**: https://orcestra.ca/pset/5ef3659e785cc307d861e795
- **Access Method**: Download from url
- **Citation**: See url
- **License**: [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/)
```

```markdown
## TDXd Cell Line Response Data.xlsx (Internal)
- **Name**: TDXd cell line IC50 values
- **Creation Date**: 2025-02-21
- **Source**: Cescon Lab (M.E.), via google drive
```

### data/rawdata/pdx

```markdown
## gene_tpm_normalized_matrix.xlsx (Internal)

- **Name**: PDX TPM counts
- **Access Date**: 2025-11-28
- **Source**: Cescon Lab (M.E.)
```

```markdown
## chromvar-Nov132025/ (Internal)

- **Name**: PDX ARCHE scoring outputs on 20k sites
- **Creation Date**: 2025-11-13
- **Source**: BHK Lab (J.N.)
```

```markdown
## auc.xlsx (Internal - no longer used)

- **Name**: PDX drug response data (AUC)
- **Creation Date**: 2025-07-23
- **Source**: Cescon Lab (G.F.)
```

```markdown
## mRECIST.xlsx (Internal - no longer used)

- **Name**: PDX drug response data (mRECIST)
- **Creation Date**: 2025-07-23
- **Source**: Cescon Lab (G.F.)
```

```markdown
## PDX_BR_BAR_grouped_Julia_8205_TAX_selectedmodels.xlsx (Internal - no longer used)

- **Name**: PDX drug response data (BR and BAR)
- **Creation Date**: 2025-07-23
- **Source**: Cescon Lab (G.F.)
```

```markdown
## DrugResponse_PDX.xlsx (Internal - no longer used)

- **Name**: PDX drug response data 
- **Creation Date**: 
- **Source**: Cescon Lab (S.E.G.)
```

```markdown
## PDX_Response_Sept2025.csv (Internal)

- **Name**: PDX drug response data 
- **Creation Date**: 2025-09-??
- **Source**: Cescon Lab (M.E.)
```

### data/rawdata/cfdna

```markdown
## CICADA-unfiltered/ (Internal)

- **Name**: CICADA ARCHE scores (rerun with new ARCHE subsets)
- **Creation Date**: 2025-11-28
- **Creation Method**: Bratman Lab 6Base-Seq Griffin pipeline
- **Input Data**: 6Base-seq bam files
- **Processing Scripts**: /cluster/projects/bhklab/projects/BCaATAC/Griffin
```

```markdown
## CICADA-BloodvsER-25/ (Internal)

- **Name**: CICADA ARCHE scores after filtering for 25% overlap of BloodvsER sites
- **Creation Date**: 2025-11-28
- **Creation Method**: Bratman Lab 6Base-Seq Griffin pipeline
- **Input Data**: 6Base-seq bam files
- **Processing Scripts**: /cluster/projects/bhklab/projects/BCaATAC/Griffin
```

```markdown
## CICADA-ARCHE/ (Internal)

- **Name**: CICADA ARCHE scores (old, archive)
- **Creation Date**: 2025-09-09
- **Creation Method**: Bratman Lab 6Base-Seq Griffin pipeline (S.M.)
- **Input Data**: 6Base-seq bam files
- **Processing Scripts**: /cluster/projects/bhklab/projects/BCaATAC/Griffin
```

```markdown
## preclinical-ARCHE/ (Internal)

- **Name**: Preclinical model ARCHE scores
- **Creation Date**: 2025-09-09
- **Creation Method**: Bratman Lab Griffin pipeline
- **Processing Scripts**: /cluster/projects/bhklab/projects/BCaATAC/Griffin
```

```markdown
## Collaboration_w_Tina_Griffin_Reference/ (Internal)

- **Name**: 
- **Creation Date**: 2025-07-24
- **Creation Method**: Lupien & Bratman Lab (T.K. and S.M.)
- **Access Date**: Downloaded 2025-10-09
- **url**: in log.txt
```

### data/rawdata/gmt

```markdown
## h.all.v2025.1.Hs.symbols.gmt 

- **Name**: Hallmark Gene Set
- **Download Date**: 2025-10-20
- **URL**: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
- **Access Method**: Download from url
```

```markdown
## All_MYC_Target_Signatures.gmt (Internal)

- **Name**: Additional MYC Gene Signatures
- **Creation Date**: 2025-10-22
- **Source**: Penn Lab (P.L.)
```