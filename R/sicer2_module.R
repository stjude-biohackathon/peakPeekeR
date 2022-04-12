# SICER2 peakcalling module.
sicer2UI <- function(id) {
  ns <- NS(id)
  shiny::uiOutput(ns("sicer2"))
}
# requires bedtools to use BAM files directly 
sicer2Server <- function(id, trt_bam, ctrl_bam = NULL, chrom, start, end, trt_track, ctrl_track = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      
      peak.call <- reactive({
        # Only run when button clicked.
        input$run
        
        # Peak calling function.
        .sicer2_calling(trt_bam = trt_bam, 
                       ctrl_bam = ctrl_bam, 
                       s = isolate(input$s), 
                       rt = isolate(input$rt),
                       w = isolate(input$w), 
                       f = isolate(input$f), 
                       egf = isolate(input$egf), 
                       fdr = isolate(input$fdr), 
                       g = isolate(input$g), 
                       e = isolate(input$e), 
                       outdir = tempdir()
                       )
      })
      
      output$peaks <- renderPlot(height = 175, {
        req(peak.call())
        df <- read.delim(peak.call()$peaks, header = FALSE)
        
       
        colnames(df) <- c("chrom", "start", "end", "score")
      
        
        gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
        loc <- GRanges(seqnames = chrom(), ranges = IRanges(start = c(start()), end = c(end())))
        
        gr <- subsetByOverlaps(gr, loc)
        
        tr <- autoplot(gr, aes_string(fill = "score")) + theme_clear()

        if (!is.null(ctrl_track())) {
          tracks <- c(Control = ctrl_track(), Treat = trt_track(), SICER2 = tr)
        } else {
          tracks <- c(Treat = trt_track(), SICER2 = tr)
        }
        
        tracks(tracks, heights = c(0.4,0.8,0.4)) + theme_clear()
      })
      
      output[["sicer2"]] <- renderUI({
        ns <- session$ns
        tags$div(
          id = environment(ns)[["namespace"]],
            tagList(
             fluidRow(
               column(6, 
                wellPanel(h4("SICER2"),
                  splitLayout(
                    selectInput(ns("s"),
                                label = "-s (species)",
                                selected = "hg38",
                                choices = c("hg38", "hg19", "mm10") 
                    ),
                 
                    numericInput(ns("rt"),
                                 label = "-rt (redundancy threshold)",
                                 value = 1,
                                 min = 1
                    ),
                    numericInput(ns("w"),
                                 label = "-w (window size)",
                                 value = 200
                    ),
                    numericInput(ns("f"),
                                 label = "-f (fragment size",
                                 value = 150
                    )
                  ),
                   splitLayout(
                    numericInput(ns("egf"),
                                 label = "-egf (effective genome fraction)",
                                 value = 0.74,
                                 min = 0.00000000000001
                    ),
                    numericInput(ns("fdr"),                         
                                 label = "-fdr (false discovery rate)",
                                 value = 0.01,
                                 min = 0.000000000001
                    ),
                    numericInput(ns("g"),
                                 label = "-g (gap size, must be multiple of window size)",
                                 value = 600
                    ),
                    numericInput(ns("e"),
                                 label = "-e e value",
                                 value = 1000
                    ),
                  ),
                  splitLayout(
                    actionButton(ns("run"), label = "Run SICER2"),
                    actionButton(ns("deleteButton"),
                                 "Delete",
                                 icon = shiny::icon("times"))
                  )
                )
              ),
              column(6,
                    plotOutput(ns("peaks"), height = 175)
              )
            )
          )
        )
      })
      
    }
  )
}

# SICER2 peak calling function.
.sicer2_calling <- function(trt_bam, ctrl_bam, s, rt, w, f, egf, 
                           fdr, g, e, outdir) {
  
  args <- c("sicer", "-t", trt_bam, "-o", outdir)
  
  # Add control bam if present.
  if (!is.null(ctrl_bam)) {
    args <- c(args, "-c", ctrl_bam)
  }
  
  
  # Add args that will always have a value. 
  args <- c(args, "-s", s, "rt", rt, "w", w, "-f", f, "-egf", egf, "-fdr", fdr, "-g", g, "-e", e)
  
  
  cl <- basiliskStart(env_sicer2)
  basiliskRun(cl, function(calling.args) {
    system2("sicer2", args = calling.args)
  }, calling.args=args)
  basiliskStop(cl)
  
    outfile <- paste0(outdir, "/",trt_bam, ".", s, "-W200-G600.scoreisland")
  
  full.cmd <- paste("sicer2", args)
  
  return(list(peaks = outfile, cmd = full.cmd))
}



