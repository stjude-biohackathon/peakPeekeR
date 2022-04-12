# MACS peakcalling module.
macsUI <- function(id) {
  ns <- NS(id)
  shiny::uiOutput(ns("macs"))
}

macsServer <- function(id, trt_bam, ctrl_bam = NULL, chrom, start, end, trt_track, ctrl_track = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      
      peak.call <- reactive({
        # Only run when button clicked.
        input$run
        
        # Peak calling function.
        .macs_calling(trt_bam = trt_bam, 
                       ctrl_bam = ctrl_bam, 
                       g = isolate(input$g), 
                       shiftsize = isolate(input$shiftsize), 
                       p = isolate(input$p), 
                       slocal = isolate(input$slocal), 
                       llocal = isolate(input$llocal), 
                       no.lambda = isolate(input$no.lambda), 
                       outdir = tempdir(), 
                       n = id)
      })
      
      output$peaks <- renderPlot(height = 175, {
        req(peak.call())
        df <- read.delim(peak.call()$peaks, header = FALSE)
        
        colnames(df) <- c("chrom", "start", "end", "name", "score")
        
        
        gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
        loc <- GRanges(seqnames = chrom(), ranges = IRanges(start = c(start()), end = c(end())))
        
        gr <- subsetByOverlaps(gr, loc)
        
        tr <- autoplot(gr, aes_string(fill = "score")) + theme_clear()

        if (!is.null(ctrl_track())) {
          tracks <- c(Control = ctrl_track(), Treat = trt_track(), MACS = tr)
        } else {
          tracks <- c(Treat = trt_track(), MACS = tr)
        }
        
        tracks(tracks, heights = c(0.4,0.8,0.4)) + theme_clear()
      })
      
      output[["macs"]] <- renderUI({
        ns <- session$ns
        tags$div(
          id = environment(ns)[["namespace"]],
            tagList(
             fluidRow(
               column(6, 
                wellPanel(h4("MACS"),
                  splitLayout(
                    numericInput(ns("g"),
                                 label = "-g (Genome Size)",
                                 value = 2700000000,
                                 min = 1
                    ),
                    numericInput(ns("shiftsize"),
                                 label = "--shiftsize",
                                 value = 100
                    )
                  ),
                   splitLayout(
                    numericInput(ns("p"),                         
                                 label = "-p",
                                 value = 0.05,
                                 step = 0.001
                    )
                  ),
                  splitLayout(
                    numericInput(ns("slocal"),
                                 label = "--slocal",
                                 value = 1000
                    ),
                    numericInput(ns("llocal"),
                                 label = "--llocal",
                                 value = 10000
                    ),
                    checkboxInput(ns("no.lambda"),
                                  label = "--nolambda",
                                  value = FALSE
                    )
                  ),
                  splitLayout(
                    actionButton(ns("run"), label = "Run MACS"),
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

# MACS peak calling function.
.macs_calling <- function(trt_bam, ctrl_bam, g, shiftsize, 
                           p, slocal, llocal, 
                           no.lambda, outdir, n) {
  
  args <- c("callpeak", "-t", trt_bam, "--outdir", outdir, "-n", n, "-p", p)
  
  # Add control bam if present.
  if (!is.null(ctrl_bam)) {
    args <- c(args, "-c", ctrl_bam)
  }
  
  # Add args that will always have a value. MACS revised
  args <- c(args, "-g", g, "--nomodel", "--shiftsize", shiftsize, 
            "--slocal", slocal, "--llocal", llocal)
  
  # Add optional/dependent args.
  if (no.lambda) {
    args <- c(args, "--nolambda")
  }
  
  cl <- basiliskStart(env_macs)
  basiliskRun(cl, function(calling.args) {
    system2("macs", args = calling.args)
  }, calling.args=args)
  basiliskStop(cl)
  

  outfile <- paste0(outdir, "/", n, "_peaks.bed")
  
  
  full.cmd <- paste("macs", args)
  
  return(list(peaks = outfile, cmd = full.cmd))
}



