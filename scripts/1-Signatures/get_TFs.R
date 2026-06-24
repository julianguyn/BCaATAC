# helper script to grab TFs for chromvar

library(readxl)

# make dataframe to store results
tf <- data.frame(matrix(nrow=0, ncol=7))

for (arche in paste0("ARCHE", 1:6)) {
    label <- paste(gsub("RCHE", "", arche), "K", sep = " - ")

    file <- "data/results/data/1-Signatures/findMotifsGenome/findMotifsGenome.xlsx"
    motifs <- read_excel(file, sheet = label) |> as.data.frame()

    # get variables and format
    motifs$Name <- sub("\\/.*", "", motifs$Name)
    colnames(motifs) <- gsub(" Sequences with Motif", "", colnames(motifs))
    motifs <- motifs[,c("Name", "Rank", "P-value", "% of Targets", "% of Background")]
    motifs$Rank <- factor(motifs$Rank, levels = motifs$Rank)

    # compute fold change
    motifs$FoldChange <- motifs$'% of Targets' / motifs$'% of Background'

    # split name
    motifs$family <- sub("\\)", "", sub(".*\\(", "", motifs$Name))
    motifs$Name <- sub("\\(.*", "", motifs$Name)

    motifs$ARCHE <- arche

    tf <- rbind(tf, motifs)
}

# 89 total, 59 unique

# get unique genes
genes <- unique(tf$Name)

motif_to_gene <- list(
  "Sp1"           = "SP1",
  "KLF17"         = "KLF17",
  "ETS"           = NA,             # family motif, not a single gene
  "Elk4"          = "ELK4",
  "KLF3"          = "KLF3",
  "Klf15"         = "KLF15",
  "NFY"           = NA,             # complex of NFYA/NFYB/NFYC, not a single gene
  "Elk1"          = "ELK1",
  "Sp5"           = "SP5",
  "ELF1"          = "ELF1",
  "KLF1"          = "KLF1",
  "Fli1"          = "FLI1",
  "GFY-Staf"      = "ZNF143",       # composite motif; Staf component = ZNF143, GFY unresolved
  "GABPA"         = "GABPA",
  "GRHL2"         = "GRHL2",
  "NF1-halfsite"  = "NFIX",         # half-site variant of NF1 family motif; ambiguous among NFIA/NFIB/NFIC/NFIX
  "Sox3"          = "SOX3",
  "Sox9"          = "SOX9",
  "Sox2"          = "SOX2",
  "Sox10"         = "SOX10",
  "Sox17"         = "SOX17",
  "Sox15"         = "SOX15",
  "NF1"           = NA,             # family motif; ambiguous among NFIA/NFIB/NFIC/NFIX
  "Sox4"          = "SOX4",
  "Sox21"         = "SOX21",
  "Sox6"          = "SOX6",
  "SOX1"          = "SOX1",
  "Sox7"          = "SOX7",
  "Tlx?"          = NA,             # HOMER flags as low-confidence; could be TLX1/TLX2/NR2E1
  "FOXM1"         = "FOXM1",
  "FOXA1"         = "FOXA1",
  "Foxa2"         = "FOXA2",
  "Fox:Ebox"      = NA,             # composite motif (FOX site + E-box), not a single gene
  "Foxa3"         = "FOXA3",
  "AP-2alpha"     = "TFAP2A",
  "AP-2gamma"     = "TFAP2C",
  "FoxL2"         = "FOXL2",
  "PHA-4"         = "FOXA2",        # C. elegans gene; closest mammalian ortholog (FOXA family)
  "FoxD3"         = "FOXD3",
  "Foxo3"         = "FOXO3",
  "FOXK1"         = "FOXK1",
  "Foxf1"         = "FOXF1",
  "FOXP1"         = "FOXP1",
  "FOXK2"         = "FOXK2",
  "CTCF"          = "CTCF",
  "PU.1"          = "SPI1",
  "ETS1"          = "ETS1",
  "SpiB"          = "SPIB",
  "Elf4"          = "ELF4",
  "Etv2"          = "ETV2",
  "ETV1"          = "ETV1",
  "ERG"           = "ERG",
  "Ets1-distal"   = "ETS1",         # distal genomic site variant of the ETS1 motif
  "ELF5"          = "ELF5",
  "IRF8"          = "IRF8",
  "ETV4"          = "ETV4",
  "ELF3"          = "ELF3",
  "BORIS"         = "CTCFL"
)

unique_genes <- unique(na.omit(unlist(motif_to_gene)))
writeLines(unique_genes, "data/procdata/TCGA/homer_TF_genes.txt")
