####################################################################
# Peak Set Creation Methods
####################################################################
source("utils/AnnotationGenome.R")
source("utils/ValidationUtils.R")
source("utils/LoggerUtils.R")
source("utils/GRangesUtils.R")
source("utils/HiddenUtils.R")
source("utils/HelperUtils.R")
source("utils/Utility.R")
library(magrittr)
library(data.table)

addReproduciblePeakSet <- function(
  excludeChr = c("chrM"),
  genomeSize = 2.7e9, 
  reproducibility = "2",
  extendSummits = 250,
  promoterRegion = c(1000, 100),
  plot = TRUE,
  threads = 8,
  parallelParam = NULL,
  force = FALSE,
  verbose = TRUE,
  genomeAnnotation = createGenomeAnnotation(genome_info),
  geneAnnotation = createGeneAnnotation(genome_info),
  
  logFile = createLogFile("addReproduciblePeakSet"),
  ...
){
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  .validInput(input = genomeSize, name = "genomeSize", valid = c("character", "numeric", "null"))
  .validInput(input = extendSummits, name = "extendSummits", valid = c("integer"))
  .validInput(input = reproducibility, name = "reproducibility", valid = c("character"))
  .validInput(input = promoterRegion, name = "promoterRegion", valid = c("integer"))
  .validInput(input = plot, name = "plot", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
  geneAnnotation <- .validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation)
  
   tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "ReproduciblePeakSet Args", logFile=logFile)
  

    
   
    #####################################################
    # BSgenome for Add Nucleotide Frequencies!
    #####################################################
    .requirePackage(genomeAnnotation$genome)
    .requirePackage("Biostrings",source="bioc")
    BSgenome <- eval(parse(text = genomeAnnotation$genome))
    BSgenome <- validBSgenome(BSgenome)
    
    #####################################################
    # Identify Reproducible Peaks!
    #####################################################
    .logDiffTime("Identifying Reproducible Peaks!", tstart, verbose = verbose, logFile = logFile)
    groupPeaks <- .safelapply(seq_along(outSummitList), function(i){
      prefix <- sprintf("Rep. Peaks Group (%s of %s) :", i, length(outSummitList))
      .logDiffTime(sprintf("%s Creating Reproducible Peaks", prefix), tstart, verbose = FALSE, logFile = logFile)
      peaks <- suppressMessages(.identifyReproduciblePeaks(
        summitFiles = outSummitList[[i]],
        summitNames = summitNamesList[[i]],
        reproducibility = reproducibility,
        extendSummits = extendSummits,
        blacklist = genomeAnnotation$blacklist,
        prefix = prefix,
        logFile = logFile
      ))
      .logDiffTime(sprintf("%s Annotating and Filtering Peaks", prefix), tstart, verbose = FALSE, logFile = logFile)
      peaks <- sort(sortSeqlevels(peaks))
      peaks <- subsetByOverlaps(peaks, genomeAnnotation$chromSizes, type = "within")
      peaks <- .fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion, logFile = logFile)
      peaks <- peaks[which(mcols(peaks)$N < 0.001)] #Remove N Containing Peaks
      peaks <- peaks[order(peaks$groupScoreQuantile, decreasing = TRUE)]
      #peaks <- head(peaks, groupSummary[names(outSummitList)[i],"maxPeaks"])
      mcols(peaks)$N <- NULL #Remove N Column
      print(file.path(outDir, paste0(make.names(summitNamesList[[i]]), "-reproduciblePeaks.gr.rds")))
      saveRDS(peaks, file.path(outDir, paste0(make.names(summitNamesList[[i]]), "-reproduciblePeaks.gr.rds")))
      return(peaks)
    }, threads = threads) %>% GRangesList()
    names(groupPeaks) <- names(outSummitList)
    
    #Construct Union Peak Set
    .logDiffTime("Creating Union Peak Set!", tstart, verbose = verbose, logFile = logFile)
    unionPeaks <- unlist(groupPeaks)
    unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)
    
    #Summarize Output
    peakDF <- lapply(seq_along(groupPeaks), function(x){
      data.frame(Group = names(groupPeaks)[x], table(groupPeaks[[x]]$peakType))
    }) %>% Reduce("rbind", .)
    print("print peaks")
    #peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
    peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
    peakDF$Freq <- peakDF$Freq / 1000
    metadata(unionPeaks)$PeakCallSummary <- peakDF
    fwrite(data.frame(unionPeaks),paste0(outDir,"_unionPeakSet.bed"),sep="\t")

  .logDiffTime(sprintf("Finished Creating Union Peak Set (%s)!", length(unionPeaks)), tstart, verbose = verbose, logFile = logFile)

}


#####################
# Utility Functions
#####################

.fastAnnoPeaks <- function(
  peaks = NULL, 
  BSgenome = NULL, 
  geneAnnotation = NULL, 
  promoterRegion = c(1000, 100),
  logFile = NULL
){
  
  #Validate
  peaks <- .validGRanges(peaks)
  peakSummits <- resize(peaks,1,"center")
  geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
  geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
  geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
  BSgenome <- validBSgenome(BSgenome)
  
  #First Lets Get Distance to Nearest Gene Start
  .logMessage("Annotating Peaks : Nearest Gene", logFile = logFile)
  distPeaks <- distanceToNearest(peakSummits, resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
  mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
  .logMessage("Annotating Peaks : Gene", logFile = logFile)
  promoters <- extendGR(resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
  op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
  og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
  oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
  type <- rep("Distal", length(peaks))
  type[which(og & oe)] <- "Exonic"
  type[which(og & !oe)] <- "Intronic"
  type[which(op)] <- "Promoter"
  mcols(peaks)$peakType <- type
  
  #First Lets Get Distance to Nearest TSS's
  .logMessage("Annotating Peaks : TSS", logFile = logFile)
  distTSS <- distanceToNearest(peakSummits, resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToTSS <- mcols(distTSS)$distance
  if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distPeaks)]
  }else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distPeaks)]
  }
  
  #Get NucleoTide Content
  .logMessage("Annotating Peaks : GC", logFile = logFile)
  nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  peaks
  
}

.identifyReproduciblePeaks <- function(
  summitFiles = NULL,
  summitNames = NULL,
  reproducibility = 0.51,
  extendSummits = 250,
  blacklist = NULL,
  prefix = NULL,
  logFile = NULL
){
  
  errorList <- mget(names(formals()),sys.frame(sys.nframe()))
  
  nonOverlapPassES <- tryCatch({
    
    .logMessage(paste0(prefix, " Getting Summits"), logFile = logFile)
    summits <- lapply(seq_along(summitFiles), function(x){
      
      #Read Summits!
      out <- data.table::fread(summitFiles, select = c(1,2,3,5))
      out <- GRanges(out$V1, IRanges(out$V2 + 1, out$V3), score = out$V5)
      grx <- out
      grx <- subsetByOverlaps(grx, blacklist, invert = TRUE) #Not Overlapping Blacklist!
      grx$GroupReplicate <- paste0(summitNames[x])
      grx
    })
    summits <- Reduce("c", as(summits, "GRangesList"))
    
    .logMessage(paste0(prefix, " Extending Summits"), logFile = logFile)
    extendedSummits <- resize(summits, extendSummits * 2 + 1, "center")
    extendedSummits <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
      nonES <- nonOverlappingGR(x, by = "score", decreasing = TRUE)
      nonES$replicateScoreQuantile <- round(.getQuantiles(nonES$score),3)
      nonES
    })
    extendedSummits <- Reduce("c", as(extendedSummits, "GRangesList"))
    
    .logMessage(paste0(prefix, " Creating Non-Overlapping Peaks"), logFile = logFile)
    nonOverlapES <- nonOverlappingGR(extendedSummits, by = "replicateScoreQuantile", decreasing = TRUE)
    
    .logMessage(paste0(prefix, " Identifying Reproducible Peaks"), logFile = logFile)
    overlapMat <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
      overlapsAny(nonOverlapES, x)
    }) %>% Reduce("cbind", .)
    
    if(length(summitFiles) > 1){
      nonOverlapES$Reproducibility <- rowSums(overlapMat)
      nonOverlapES$ReproducibilityPercent <- round(rowSums(overlapMat) / ncol(overlapMat) , 3)
      n <- length(summitFiles)
      minRep <- eval(parse(text=reproducibility))
      if(!is.numeric(minRep)){
        stop("Error reproducibility not numeric when evaluated!")
      }
      idxPass <- which(nonOverlapES$Reproducibility >= minRep)
      nonOverlapPassES <- nonOverlapES[idxPass]
    }else{
      nonOverlapES$Reproducibility <- rep(NA, length(nonOverlapES))
      nonOverlapPassES <- nonOverlapES
    }
    
    .logMessage(paste0(prefix, " Finalizing Peaks"), logFile = logFile)
    nonOverlapPassES$groupScoreQuantile <- round(.getQuantiles(nonOverlapPassES$replicateScoreQuantile),3)
    mcols(nonOverlapPassES) <- mcols(nonOverlapPassES)[,c("score","replicateScoreQuantile", "groupScoreQuantile", "Reproducibility", "GroupReplicate")]
    
    nonOverlapPassES
    
  }, error = function(e){
    
    .logError(e, fn = ".identifyReproduciblePeaks", info = prefix, errorList = errorList, logFile = logFile) 
    
  })
  
  return(nonOverlapPassES)
  
}


#####################################################
    # Create Output Directory
#####################################################
outDir <- file.path("results", "PeakCalls")
outSubDir <- file.path("results", "ReplicateCalls")
outBedDir <- file.path("results", "InsertionBeds")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
dir.create(outSubDir, showWarnings = FALSE, recursive = TRUE)
dir.create(outBedDir, showWarnings = FALSE, recursive = TRUE)

#################################
#Summit lists
path1="/mnt/work1/users/lupiengroup/People/ankita/3_tcga/summits"
outSummitList1 <- list.files(path1,pattern ="_summits.bed$",full.names = TRUE)
summitNamesList1 <- gsub("_summits.bed","",list.files(path1,pattern ="_summits.bed$"))


#path2="/mnt/work1/users/lupiengroup/People/ankita/3_tcga/prad"
#outSummitList2 <- list.files(path2,pattern ="_summits.bed$",full.names = TRUE)
#summitNamesList2 <- gsub("_summits.bed","",list.files(path,pattern ="_summits.bed$"))
#outSummitList<-c(outSummitList1,outSummitList2)
#summitNamesList<-c(summitNamesList1,summitNamesList2)

outSummitList<-outSummitList1
summitNamesList<-summitNamesList1

names(outSummitList)<-summitNamesList

#head(outSummitList)
head(summitNamesList)

genome_info<-"hg38"
addReproduciblePeakSet()
