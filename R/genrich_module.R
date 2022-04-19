# Genrich peakcalling module.
genrichUI <- function(id) {
  ns <- NS(id)
  shiny::uiOutput(ns("genrich"))
}

genrichServer <- function(id, trt_bam, ctrl_bam = NULL, chrom, start, end, trt_track, ctrl_track = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      
      peak.call <- reactive({
        # Only run when button clicked.
        input$run
        
        # Peak calling function.
        .genrich_calling(t = trt_bam, 
                         c = ctrl_bam, 
                         m = isolate(input$m), 
                         s = isolate(input$s), 
                         r = isolate(input$r), 
                         y = isolate(input$y), 
                         w = isolate(input$w), 
                         x = isolate(input$x), 
                         j = isolate(input$j), 
                         d = isolate(input$d), 
                         D = isolate(input$D),
                         q = isolate(input$q), 
                         p = isolate(input$p), 
                         l = isolate(input$l), 
                         a = isolate(input$a), 
                         g = isolate(input$g), 
                         o = paste0(tempdir(), "/", id, ".narrowPeak"))
      })
      
      output$peaks <- renderPlot(height = 200, {
        req(peak.call())
        
        out <- tryCatch(
          {
            df <- read.delim(peak.call()$peaks, header = FALSE)
            
            colnames(df) <- c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
            
            gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
            loc <- GRanges(seqnames = chrom(), ranges = IRanges(start = c(start()), end = c(end())))
            
            gr <- subsetByOverlaps(gr, loc)
            
            tr <- autoplot(gr, aes_string(fill = "score")) + theme_clear()
            
            if (!is.null(ctrl_track())) {
              tracks <- c(Ctrl = ctrl_track(), Treat = trt_track(), Genrich = tr)
            } else {
              tracks <- c(Treat = trt_track(), Genrich = tr)
            }
            
            tracks(tracks, heights = c(0.4,0.8,0.4)) + theme_clear()
          },
          error=function(cond) {
            grid.newpage()
            grid.text(paste0(cond, "\n Likely unpaired sequencing, use -y option."))
          }
        )
        out
        
      })
      
      output[["genrich"]] <- renderUI({
        ns <- session$ns
        tags$div(
          id = environment(ns)[["namespace"]],
          tagList(
            fluidRow(
              column(12, 
                     wellPanel(
                       fluidRow(
                         column(12,
                                div(h4("Genrich"), style = 'float:left'),
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
                                       numericInput(ns("m"),
                                                    label = "-m",
                                                    value = 0,
                                                    step = 1
                                       ),
                                       numericInput(ns("s"),
                                                    label = "-s",
                                                    value = 0,
                                                    step = 0.1
                                       ),
                                       numericInput(ns("w"),
                                                    label = "-w",
                                                    value = 0,
                                                    step = 1
                                       ),
                                ),
                                column(3,
                                       numericInput(ns("d"),
                                                    label = "-d",
                                                    value = 100,
                                                    step = 1
                                       ),
                                       numericInput(ns("p"),
                                                    label = "-p",
                                                    value = 0.01,
                                                    step = 0.01
                                       ),
                                       numericInput(ns("q"),
                                                    label = "-q",
                                                    value = 1,
                                                    step = 0.01
                                       ),
                                ),
                                column(3,
                                       numericInput(ns("a"),
                                                    label = "-a",
                                                    value = 200.0,
                                                    step = 0.1
                                       ),
                                       numericInput(ns("l"),
                                                    label = "-l",
                                                    value = 0,
                                                    step = 1
                                       ),
                                       numericInput(ns("g"),
                                                    label = "-g",
                                                    value = 100,
                                                    step = 1
                                       )
                                ),
                                column(3,
                                       checkboxInput(ns("r"),
                                                     label = "-r",
                                                     value = FALSE
                                       ),
                                       checkboxInput(ns("y"),
                                                     label = "-y",
                                                     value = FALSE
                                       ),
                                       checkboxInput(ns("x"),
                                                     label = "-x",
                                                     value = FALSE
                                       ),
                                       checkboxInput(ns("j"),
                                                     label = "-j",
                                                     value = FALSE
                                       ),
                                       checkboxInput(ns("D"),
                                                     label = "-D",
                                                     value = FALSE
                                       ),
                                       div(
                                         actionButton(ns("run"), label = "Run Genrich", class='btn-success', style='margin-top:20px;'), style="float:right;"
                                         
                                       )
                                )
                         )
                       )
                     )
              )
            )
          )
        )
      })
      
    }
  )
}

# Genrich peak calling function.
.genrich_calling <- function(t, c, m, s, r, y, w, x, j, d, D,
                             q, p, l, a, g, o) {
  
  args <- c("-t", t, "-o", o)
  
  # Add control bam if present.
  if (!is.null(c)) {
    args <- c(args, "-c", c)
  }
  
  # Filtering options.
  if (y) {
    args <- c(args, "-y")
  }
  
  if (m != 0) {
    args <- c(args, "-m", m)
  }
  
  if (s != 0) {
    args <- c(args, "-s", s)
  }
  
  if (w != 0) {
    args <- c(args, "-w", w)
  }
  
  if (r) {
    args <- c(args, "-r")
  }
  
  # ATAC-seq options
  if (j) {
    args <- c(args, "-j")
  }
  
  if (d != 0) {
    args <- c(args, "-d", d)
  }
  
  if (D) {
    args <- c(args, "-D")
  }
  
  # Peak-calling options
  if (q != 1) {
    args <- c(args, "-q")
  }
  
  args <- c(args, "-p", p, "-l", l, "-g", g, "-a", a, "-v", "-e", "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
  
  cl <- basiliskStart(env_genrich)
  basiliskRun(cl, function(calling.args) {
    system2("Genrich", args = calling.args)
  }, calling.args=args)
  basiliskStop(cl)
  
  full.cmd <- paste(c("Genrich", args))
  
  return(list(peaks = o, cmd = full.cmd))
}



