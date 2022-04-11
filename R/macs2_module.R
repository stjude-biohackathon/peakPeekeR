# MACS2 peakcalling module.
macs2UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6, 
        wellPanel(
          splitLayout(
            selectInput(ns('f'),
                      label = '-f (Paired/Single End)',
                      selected = "BAM",
                      choices = c("BAM", "BAMPE")
            ),
            checkboxInput(ns('broad'),
                        label = '--broad (call broad peaks)',
                        value = FALSE
            ),
            numericInput(ns('g'),
                         label = "-g (Genome Size)",
                         value = 2700000000,
                         min = 1
            ),
            numericInput(ns('extsize'),
                         label = "--extsize",
                         value = 150
            ),
            numericInput(ns('shift'),
                         label = "--shift",
                         value = 75
            )
          ),
          splitLayout(
            numericInput(ns('q'),
                      label = '-q',
                      value = 0.05,
                      min = 0.00000000000001
            ),
            numericInput(ns('p'),                         
                         label = "-p",
                         value = 1,
                         min = 0.000000000001
            ),
            numericInput(ns('min.length'),
                         label = "--min-length",
                         value = 0
            ),
            numericInput(ns('max.gap'),
                         label = "--max-gap",
                         value = 0
            ),
            numericInput(ns('broad.cutoff'),
                         label = '--broad-cutoff',
                         value = 0.1,
                         min = 0.00000000000001
            )
          ),
          splitLayout(
            numericInput(ns('slocal'),
                         label = '--slocal',
                         value = 1000
            ),
            numericInput(ns('llocal'),
                         label = "--llocal",
                         value = 10000
            ),
            checkboxInput(ns("no.lambda"),
                        label = "--nolambda",
                        value = FALSE
            )
          ),
          div(actionButton(ns("run"), label = "Run MACS2"), align = "right")
        )
      ),
      column(6,
        plotOutput(ns("peaks"), height = 175)
      )
    )
  )
}

macs2Server <- function(id, trt_bam, ctrl_bam = NULL, chrom, start, end, trt_track, ctrl_track = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      
      peak.call <- reactive({
        # Peak calling function.
        .macs2_calling(trt_bam = trt_bam, 
                       ctrl_bam = ctrl_bam, 
                       f = input$f, 
                       broad = input$broad, 
                       g = input$g, 
                       extsize = input$extsize, 
                       shift = input$shift, 
                       q = input$q, 
                       p = input$p, 
                       min.length = input$min.length, 
                       max.gap = input$max.gap, 
                       broad.cutoff = input$broad.cutoff, 
                       slocal = input$slocal, 
                       llocal = input$llocal, 
                       no.lambda = input$no.lambda, 
                       outdir = tempdir(), 
                       n = id)
      })
      
      output$peaks <- renderPlot(height = 175, {
        req(peak.call())
        df <- read.delim(peak.call()$peaks, header = FALSE)
        colnames(df) <- c("chrom", "start", "end", "name", "score", "strand", "signal", "p", "q", "peak")
        gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
        loc <- GRanges(seqnames = chrom(), ranges = IRanges(start = c(start()), end = c(end())))
        browser()
        tr <- autoplot(gr, aes(fill = score), which = loc) + theme_clear()

        if (!is.null(ctrl_track())) {
          tracks <- c(Control = ctrl_track(), Treat = trt_track(), MACS2 = tr)
        } else {
          tracks <- c(Treat = trt_track(), MACS2 = tr)
        }
        
        tracks(tracks, heights = c(0.4,0.8,0.4)) + theme_clear()
      })
      
    }
  )
}

# MACS2 peak calling function.
.macs2_calling <- function(trt_bam, ctrl_bam, f, broad, g, extsize, shift, 
                           q, p, min.length, max.gap, broad.cutoff, slocal, llocal, 
                           no.lambda, outdir, n) {
  
  args <- c("callpeak", "-t", trt_bam, "--outdir", outdir, '-n', n)
  
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
  
  system2("macs2", args = args)
  
  if (broad) {
    outfile <- paste0(outdir, "/", n, "_peaks.broadPeak")
  } else {
    outfile <- paste0(outdir, "/", n, "_peaks.narrowPeak")
  }
  
  full.cmd <- paste("macs", args)
  
  return(list(peaks = outfile, cmd = full.cmd))
}



