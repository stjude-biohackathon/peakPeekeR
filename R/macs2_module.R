# MACS2 peakcalling module.
macs2UI <- function(id) {
  ns <- NS(id)
  shiny::uiOutput(ns("macs2"))
}

macs2Server <- function(id, trt_bam, ctrl_bam = NULL, chrom, start, end, trt_track, ctrl_track = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      
      peak.call <- reactive({
        # Only run when button clicked.
        input$run
        
        # Peak calling function.
        .macs2_calling(trt_bam = trt_bam, 
                       ctrl_bam = ctrl_bam, 
                       f = isolate(input$f), 
                       broad = isolate(input$broad), 
                       g = isolate(input$g), 
                       extsize = isolate(input$extsize), 
                       shift = isolate(input$shift), 
                       q = isolate(input$q), 
                       p = isolate(input$p), 
                       min.length = isolate(input$min.length), 
                       max.gap = isolate(input$max.gap), 
                       broad.cutoff = isolate(input$broad.cutoff), 
                       slocal = isolate(input$slocal), 
                       llocal = isolate(input$llocal), 
                       no.lambda = isolate(input$no.lambda), 
                       outdir = tempdir(), 
                       n = id)
      })
      
      output$peaks <- renderPlot(height = 200, {
        req(peak.call())
        df <- read.delim(peak.call()$peaks, header = FALSE)
        
        if (isolate(input$broad)) {
          colnames(df) <- c("chrom", "start", "end", "name", "score", "strand", "signal", "p", "q")
        } else {
          colnames(df) <- c("chrom", "start", "end", "name", "score", "strand", "signal", "p", "q", "peak")
        }
        
        gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
        loc <- GRanges(seqnames = chrom(), ranges = IRanges(start = c(start()), end = c(end())))
        
        gr <- subsetByOverlaps(gr, loc)
        
        tr <- autoplot(gr, aes_string(fill = "score")) + theme_clear()
        
        if (!is.null(ctrl_track())) {
          tracks <- c(Ctrl = ctrl_track(), Treat = trt_track(), MACS2 = tr)
        } else {
          tracks <- c(Treat = trt_track(), MACS2 = tr)
        }
        
        tracks(tracks, heights = c(0.4,0.8,0.5)) + theme_clear()
      })
      
      output[["macs2"]] <- renderUI({
        ns <- session$ns
        tags$div(
          id = environment(ns)[["namespace"]],
          tagList(
            fluidRow(
              column(12, 
                     wellPanel(
                       fluidRow(
                         column(12,
                                div(h4("MACS2"), style = 'float:left'),
                                div(actionButton(ns("deleteButton"),
                                                 "Delete",
                                                 icon = shiny::icon("times"), class='btn-danger'), style = 'float:right;')
                         ),
                         column(12, hr()),
                         column(5, 
                                plotOutput(ns("peaks"), height = 200)
                         ),
                         column(7,
                                column(3,
                                       numericInput(ns("g"),
                                                    label = "g",
                                                    value = 2700000000,
                                                    min = 1
                                       ),
                                       numericInput(ns("min.length"),
                                                    label = "min-length",
                                                    value = 0
                                       ),
                                       numericInput(ns("slocal"),
                                                    label = "slocal",
                                                    value = 1000
                                       ),
                                ),
                                column(3,
                                       numericInput(ns("extsize"),
                                                    label = "extsize",
                                                    value = 150
                                       ),
                                       numericInput(ns("max.gap"),
                                                    label = "max-gap",
                                                    value = 0
                                       ),
                                       numericInput(ns("llocal"),
                                                    label = "llocal",
                                                    value = 10000
                                       )
                                ),
                                
                                column(3,
                                       numericInput(ns("shift"),
                                                    label = "shift",
                                                    value = 75
                                       ),
                                       numericInput(ns("p"),                         
                                                    label = "p",
                                                    value = 1,
                                                    min = 0.000000000001
                                       ),
                                       numericInput(ns("broad.cutoff"),
                                                    label = "broad-cutoff",
                                                    value = 0.1,
                                                    min = 0.00000000000001
                                       )
                                ),
                                column(3,
                                       selectInput(ns("f"),
                                                   label = "f",
                                                   selected = "BAM",
                                                   choices = c("BAM", "BAMPE")
                                       ),
                                       
                                       numericInput(ns("q"),
                                                    label = "q",
                                                    value = 0.05,
                                                    min = 0.00000000000001
                                       ),
                                       column(6,
                                              checkboxInput(ns("no.lambda"),
                                                            label = "nolambda",
                                                            value = FALSE
                                              ),
                                              checkboxInput(ns("broad"),
                                                            label = "broad",
                                                            value = FALSE
                                              )
                                       ),
                                       column(6,
                                              div(
                                                actionButton(ns("run"), label = "Run", class='btn-success', style='margin-top:20px;')
                                                
                                              )
                                       )
                                ),
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

# MACS2 peak calling function.
.macs2_calling <- function(trt_bam, ctrl_bam, f, broad, g, extsize, shift, 
                           q, p, min.length, max.gap, broad.cutoff, slocal, llocal, 
                           no.lambda, outdir, n) {
  
  args <- c("callpeak", "-t", trt_bam, "--outdir", outdir, "-n", n)
  
  # Add control bam if present.
  if (!is.null(ctrl_bam)) {
    args <- c(args, "-c", ctrl_bam)
  }
  
  # Add args that will always have a value.
  args <- c(args, "-f", f, "-g", g, "--nomodel", "--extsize", extsize, "--shift", shift, 
            "--slocal", slocal, "--llocal", llocal, "-q", q)
  
  # Add optional/dependent args.
  if (broad) {
    args <- c(args, "--broad")
    if (broad.cutoff != 0.1) {
      args <- c(args, "--broad-cutoff", broad.cutoff)
    }
  }
  
  if (no.lambda) {
    args <- c(args, "--nolambda")
  }
  
  if (min.length != 0) {
    args <- c(args, "--min-length", min.length)
  }
  
  if (max.gap != 0) {
    args <- c(args, "--max-gap", max.gap)
  }
  
  if (p != 1) {
    args <- c(args, "-p", p)
  }
  
  cl <- basiliskStart(env_macs2)
  basiliskRun(cl, function(calling.args) {
    system2("macs2", args = calling.args)
  }, calling.args=args)
  basiliskStop(cl)
  
  if (broad) {
    outfile <- paste0(outdir, "/", n, "_peaks.broadPeak")
  } else {
    outfile <- paste0(outdir, "/", n, "_peaks.narrowPeak")
  }
  
  full.cmd <- paste("macs2", args)
  
  return(list(peaks = outfile, cmd = full.cmd))
}



