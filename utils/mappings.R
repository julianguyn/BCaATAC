###########################################################
# Mappings
###########################################################

mapping_cells <- c( "AU665" = "AU565",      # misspelt in tdxd
                    "BT20" = "BT-20", 
                    "BT474" = "BT-474", 
                    "BT549" = "BT-549", 
                    "CAL120" = "CAL-120", 
                    "CAL148" = "CAL-148", 
                    "CAL51" = "CAL-51", 
                    "CAL851" = "CAL-85-1",
                    "CAMA1" = "CAMA-1", 
                    "EFM19" = "EFM-19", 
                    "EFM192A" = "EFM-192A",
                    "HBL100" = "HBL-100",
                    "HDQP1" = "HDQ-P1", 
                    "HS578T" = "Hs 578T",
                    "Hs-578T" = "Hs 578T",
                    "JIMT1" = "JIMT-1", 
                    "KPL1" = "KPL-1", 
                    "LY2" = "MCF-7/LY2", 
                    "MCF-7-LY2" = "MCF-7/LY2",
                    "MCF7" = "MCF-7", 
                    "MCF-10A" = "MCF10A",
                    "MDAMB157" = "MDA-MB-157",
                    "MDA-MB-175VII" = "MDA-MB-175-VII", 
                    "MDAMB175VII" = "MDA-MB-175-VII", 
                    "MDAMB231" = "MDA-MB-231", 
                    "MDAMB361" = "MDA-MB-361", 
                    "MDAMB436" = "MDA-MB-436", 
                    "MDAMB468" = "MDA-MB-468", 
                    "MPE600" = "600MPE",
                    "MX1" = "MX-1", 
                    "SKBR3" = "SK-BR-3", 
                    "SKBR5" = "SK-BR-5", 
                    "SKBR7" = "SK-BR-7", 
                    "SUM149" = "SUM149PT", 
                    "SUM159" = "SUM159PT",
                    "SUM52" =  "SUM52PE",
                    "T47D" = "T-47D",
                    "UACC812" = "UACC-812", 
                    "X600MPE" = "600MPE",
                    "ZR751" = "ZR-75-1")

mapping_drugs <- c("carboplatin" = "Carboplatinum", 
                    "dipyridamole" = "Dipyridamole", 
                    "epirubicin" = "Epirubicin", 
                    "eribulin" = "Eribulin", 
                    "fluvastatin" = "Fluvastatin", 
                    "herceptin" = "Trastuzumab", 
                    "lapatinib" = "Lapatinib", 
                    "paclitaxel" = "Paclitaxel",
                    "OLA" = "Olaparib",
                    "doxorubicin:navitoclax (2:1 mol/mol)" = "Doxorubicin:Navitoclax", 
                    "carboplatin:UNC0638 (2:1 mol/mol)" = "Carboplatin:UNC0638", 
                    "tanespimycin:gemcitabine (1:1 mol/mol)" = "Tanespimycin:Gemcitabine", 
                    "carboplatin:etoposide (40:17 mol/mol)" = "Carboplatin:Etoposide", 
                    "serdemetan:SCH-529074 (1:1 mol/mol)" = "Serdemetan:SCH-529074", 
                    "navitoclax:pluripotin (1:1 mol/mol)" = "Navitoclax:Pluripotin", 
                    "Carbamic acid, N,N'-(1,2,3,4-tetrahydro-6,7,8-trimethoxy-4-oxo-2-quinazolinylidene)bis-, dimethyl ester" = "Carbamic acid", 
                    "navitoclax:gemcitabine (1:1 mol/mol)" = "Navitoclax:Gemcitabine",
                    "docetaxel:tanespimycin (2:1 mol/mol)" = "Docetaxel:Tanespimycin",
                    "navitoclax:MST-312 (1:1 mol/mol)" = "Navitoclax:MST-312",
                    "decitabine:carboplatin (1:1 mol/mol)" = "Decitabine:Carboplatin"
)

mapping_pdxs <- c("X108099P1" = "108099",
                  "108099P1" = "108099",
                  "X16720" = "16720",
                  "X48602" = "48602",
                  "X53782" = "53782",
                  "BPTO95" = "BPTO.95",
                  "INS_B014" = "INSB014",
                  "INS_B019" = "INSB019",
                  "REFS032" = "REF032",
                  "REFS034P1" = "REF034",
                  "REFS034" = "REF034",
                  "REFS036" = "REF036",
                  "REFS038" = "REF038",
                  "REFB047" = "REF047",
                  "REF047_S26_L002" = "REF047",
                  "REFS036P1A" = "REF036",
                  "REF_S_038" = "REF038"
)

###########################################################
# Map Functions
###########################################################

#' Check overlap of sample IDs and report missing
#' 
#' Helper function for finding mismatched sample IDs lol
#' 
#' @param list1 string. First list of model IDs
#' @param list2 string. Second list of model IDs
#' @param label1 string. Name of first list
#' @param label2 string. Name of second list
#' 
check_sample_overlap <- function(list1, list2, label1, label2) {
    list1 <- unique(list1)
    list2 <- unique(list2)
    message(paste("\nNumber of samples from", label1, "in", label2))
    print(table(list1 %in% list2))
    message(paste("\nSamples from", label1, "NOT in", label2))
    missing <- list1[-which(list1 %in% list2)]
    print(missing[order(missing)])
    message(paste("\nNumber of samples from", label2, "in", label1))
    print(table(list2 %in% list1))
    message(paste("\nSamples from", label2, "NOT in", label1))
    missing <- list2[-which(list2 %in% list1)]
    print(missing[order(missing)])
}

#' Map BCa cell lines to mapping_cells
#' 
#' @param cells_to_map string. Vector of CCLs names to map
#' @return original vector with standardized cell line names
#' 
map_cells <- function(cells_to_map) {
    for (i in 1:length(cells_to_map)) {
        cell = cells_to_map[i]
        if (cell %in% names(mapping_cells)) cells_to_map[i] <- unname(mapping_cells[cell])
    }
    return(cells_to_map)
}

#' Shorten drug names for plotting
#'
#' Reference defined mapping.
#' @param toPlot dataframe. Output from computePC()
#' @return toPlot dataframe with re-mapped drug names.
#' 

# function to standardize mapping of sensitivity dataframe
map_drugs <- function(toPlot) {

    for (i in 1:nrow(toPlot)) {
        drug = toPlot$drug[i]
        if (drug %in% names(mapping_drugs)) toPlot$drug[i] <- unname(mapping_drugs[drug])
    }

    # redo pairs
    toPlot$pairs <- paste(toPlot$signature, toPlot$drug, sep = "_")
    return(toPlot)
}

#' Standardize PSet cell and drug naming
#'
#' Reference defined mapping to standardize PSet cell and drug names.
#' @param sen dataframe. Sensitivity dataframe obtained from get_drugsen()
#' @return Sen dataframe with standardized cell and drug names.
#' 

# function to standardize mapping of sensitivity dataframe
map_sen <- function(sen) {

    colnames(sen) <- map_cells(colnames(sen))

    for (i in 1:length(rownames(sen))) {
        cell = rownames(sen)[i]
        if (cell %in% names(mapping_drugs)) rownames(sen)[i] <- unname(mapping_drugs[cell])
    }
    return(sen)
}

#' Standardize PDX Naming
#'
#' Reference defined mapping to standardize PDX names.
#' @param models string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 

map_pdx <- function(models) {

    for (i in 1:length(models)) {
        model = models[i]
        if (model %in% names(mapping_pdxs)) models[i] <- unname(mapping_pdxs[model])
    }
    return(models)
}
