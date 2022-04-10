# Subset and convert BAMs to BigWig.
.subset_bams <- function(trt_bam, ctrl_bam = NULL, chrom = "chr12", start = 6522378, end = 6769097) {
  
  # Get range to subset.
  gr <- GRanges(seqnames = chrom,
                ranges = IRanges(start = c(start), end = c(end)))
  
  params <- ScanBamParam(which = gr, what = scanBamWhat())
  
  bammy <- BamFile(trt_bam)
  tmp <- tempdir()
  ctrl_bw <- paste0(tmp, "/ctrl.bw")
  trt_bw <- paste0(tmp, "/trt.bw")
  
  if (!is.null(ctrl_bam)) {
    inny <- BamFile(ctrl_bam)
    ctrl <- readGAlignments(bammy, param = params)
    ctrl_cov <- coverage(ctrl)
    export.bw(ctrl_cov, con = ctrl_bw)
  }
  
  trt <- readGAlignments(bammy, param = params)
  trt_cov <- coverage(trt)
  
  export.bw(trt_cov, con = trt_bw)
  
  return(list("trt" = trt_bw, "ctrl" = ctrl_bw))
}
