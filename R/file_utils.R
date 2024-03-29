#' @importFrom basilisk BasiliskEnvironment
#' @importFrom Rsamtools ScanBamParam BamFile filterBam scanBamWhat sortBam
.subset_bams <- function(trt_bam, ctrl_bam = NULL, chrom = "chr12", start = 6522378, end = 6769097) {
  
  # Get range to subset.
  gr <- GRanges(seqnames = chrom,
                ranges = IRanges(start = c(start), end = c(end)))
  
  params <- ScanBamParam(which = gr, what = scanBamWhat())
  
  bammy <- BamFile(trt_bam)
  tmp <- tempdir()

  trt_sub <- paste0(tmp, "/trt_sub.bam")
  ctrl_sub <- NULL
  
  if (!is.null(ctrl_bam)) {
    inny <- BamFile(ctrl_bam)
    ctrl_sub <- paste0(tmp, "/ctrl_sub.bam")
    ctrl <- filterBam(inny, destination = ctrl_sub, param = params)
  }
  
  trt <- filterBam(bammy, destination = trt_sub, param = params)
  
  return(list("trt" = trt_sub, "ctrl" = ctrl_sub))
}

# Sort by query name, required by Genrich.
.qname_sort_bams <- function(trt_bam, ctrl_bam = NULL) {
  params <- ScanBamParam(what = scanBamWhat())
  
  bammy <- BamFile(trt_bam)
  tmp <- tempdir()
  
  trt_sub <- paste0(tmp, "/trt_sub.sorted.bam")
  ctrl_sub <- NULL
  
  if (!is.null(ctrl_bam)) {
    inny <- BamFile(ctrl_bam)
    ctrl_sub <- paste0(tmp, "/ctrl_sub.sorted.bam")
    ctrl <- sortBam(inny, destination = paste0(tmp, "/ctrl_sub.sorted"), param = params, byQname = TRUE)
  }
  
  trt <- sortBam(bammy, destination = paste0(tmp, "/trt_sub.sorted"), param = params, byQname = TRUE)
  
  return(list("trt" = trt_sub, "ctrl" = ctrl_sub))
}

.bamtobed <- function(trt_bam, ctrl_bam = NULL) {
  
  trt_bed <- paste0(trt_bam, ".bed")
  ctrl_bed <- NULL
  
  cl <- basiliskStart(env_sicer2)
  
  basiliskRun(cl, function(calling.args) {
    system2("bedtools", args = calling.args)
  }, calling.args=c("bamtobed", "-i", trt_bam, ">", trt_bed))
  
  if (!is.null(ctrl_bam)) {
    ctrl_bed <- paste0(ctrl_bam, ".bed")
    basiliskRun(cl, function(calling.args) {
      system2("bedtools", args = calling.args)
    }, calling.args=c("bamtobed", "-i", ctrl_bam, ">", ctrl_bed))
  }
  
  basiliskStop(cl)
  
  return(list("trt" = trt_bed, "ctrl" = ctrl_bed))
}
