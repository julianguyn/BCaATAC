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
                    "OLA" = "Olaparib")

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
        if (cell %in% names(mapping_cells)) {cells_to_map[i] <- unname(mapping_cells[cell])}
    }
    return(cells_to_map)
}

