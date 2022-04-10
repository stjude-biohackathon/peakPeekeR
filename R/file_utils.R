# Subset BAMs.
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
