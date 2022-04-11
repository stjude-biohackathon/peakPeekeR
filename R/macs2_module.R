# MACS2 peakcalling module.
macs2UI <- function(id) {
  ns <- NS(id)
  tagList(
    actionButton(ns("button"), label = "What do I do?"),
    verbatimTextOutput(ns("bams")),
    verbatimTextOutput(ns("location"))
  )
}

macs2Server <- function(id, trt_bam, ctrl_bam = NULL, chrom, start, end) {
  moduleServer(
    id,
    function(input, output, session) {

      
      output$bams <- renderText({
        paste0("Treat bam: ", trt_bam, "; Ctrl bam: ", ctrl_bam)
      })
      
      output$location <- renderText({
        paste0(chrom(), ":", start(), "-", end())
      })
    }
  )
}