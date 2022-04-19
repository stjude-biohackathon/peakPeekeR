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
        
        # Calculate egf based on bam size. Will use hg38 genome length, doesn't really matter.
        egf <- (end() - start()) / 3088286401
        
        # Peak calling function.
        .sicer2_calling(trt_bam = trt_bam, 
                        ctrl_bam = ctrl_bam, 
                        s = "hg38", 
                        rt = isolate(input$rt),
                        w = isolate(input$w), 
                        f = isolate(input$f), 
                        egf = egf, 
                        fdr = isolate(input$fdr), 
                        g = isolate(input$g), 
                        e = isolate(input$e), 
                        outdir = tempdir()
        )
      })
      
      output$peaks <- renderPlot(height = 200, {
        req(peak.call())
        df <- read.delim(peak.call()$peaks, header = FALSE)
        
        
        colnames(df) <- c("chrom", "start", "end", "score")
        
        
        gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
        loc <- GRanges(seqnames = chrom(), ranges = IRanges(start = c(start()), end = c(end())))
        
        gr <- subsetByOverlaps(gr, loc)
        
        tr <- autoplot(gr, aes_string(fill = "score")) + theme_clear()
        
        if (!is.null(ctrl_track())) {
          tracks <- c(Ctrl = ctrl_track(), Treat = trt_track(), SICER2 = tr)
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
              column(12, 
                     wellPanel(
                       fluidRow(
                         column(12,
                                div(h4("SICER2"), style = 'float:left'),
                                div(actionButton(ns("deleteButton"),
                                                 "Delete",
                                                 icon = shiny::icon("times"), class='btn-danger'), style = 'float:right;')
                         ),
                         column(12, hr()),
                         column(5,
                                plotOutput(ns("peaks"), height = 200)
                         ),
                         column(3,
                                numericInput(ns("rt"),
                                             label = "-rt (redundancy threshold)",
                                             value = 1,
                                             min = 1
                                ),
                                numericInput(ns("g"),
                                             label = "-g (must be multiple of window size)",
                                             value = 600
                                ),
                         ),
                         column(2,
                                numericInput(ns("w"),
                                             label = "-w",
                                             value = 200
                                ),
                                numericInput(ns("f"),
                                             label = "-f",
                                             value = 150
                                )
                         ),
                         column(2,
                                numericInput(ns("fdr"),                         
                                             label = "-fdr",
                                             value = 0.01,
                                             min = 0.000000000001
                                ),
                                numericInput(ns("e"),
                                             label = "-e",
                                             value = 1000
                                ),
                                div(
                                  actionButton(ns("run"), label = "Run SICER2", class='btn-success', style='margin-top:20px;'), style="float:right;"
                                )
                         )
                       )
                     )
              )
            )
          )
        )
      }
      )
    }
  )
}

# SICER2 peak calling function.
.sicer2_calling <- function(trt_bam, ctrl_bam, s, rt, w, f, egf, 
                            fdr, g, e, outdir) {
  
  args <- c("-t", paste0(trt_bam,".bed"), "-o", outdir)
  
  # Add control bam if present.
  if (!is.null(ctrl_bam)) {
    args <- c(args, "-c", paste0(ctrl_bam,".bed"))
  }
  
  
  # Add args that will always have a value. 
  args <- c(args, "-s", s, "-rt", rt, "-w", w, "-f", f, "-egf", egf, "-fdr", fdr, "-g", g, "-e", e)
  
  cl <- basiliskStart(env_sicer2)
  # bamtobed, sicer2 breaks when trying to do this itself.
  # TODO move outside this module to file utilities function.
  basiliskRun(cl, function(calling.args) {
    system2("bedtools", args = calling.args)
  }, calling.args=c("bamtobed", "-i", trt_bam, ">", paste0(trt_bam, ".bed")))
  
  if (!is.null(ctrl_bam)) {
    basiliskRun(cl, function(calling.args) {
      system2("bedtools", args = calling.args)
    }, calling.args=c("bamtobed", "-i", ctrl_bam, ">", paste0(ctrl_bam, ".bed")))
  }
  
  basiliskRun(cl, function(calling.args) {
    system2("sicer", args = calling.args)
  }, calling.args=args)
  basiliskStop(cl)
  
  outfile <- paste0(trt_bam, "-W200-G600.scoreisland")
  
  full.cmd <- paste("sicer", args)
  
  return(list(peaks = outfile, cmd = full.cmd))
}



