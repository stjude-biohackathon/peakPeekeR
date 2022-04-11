#' peakPeekeR
#'
#' Does the peak peeking.
#'
#' @return A Shiny app for peak calling parameter optimization.
#' @author Jared Andrews, Lawryn Kasper, Jacob Steele
#'
#' @export
#' @importFrom ggbio autoplot tracks
#' @import shiny
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomicRanges GRanges
peakPeekeR <- function(trt_bam, ctrl_bam = NULL) {
  ui <- fluidPage(
    titlePanel("peakPeekeR"),
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    fluidRow(
      column(6, wellPanel(
        splitLayout(
          textInput('chrom',
                    label = 'Chromosome',
                    value = "chr12"
          ),
          
          numericInput('start',
                    label = "Start Position",
                    value = 6522378,
                    min = 1
          ),
          
          numericInput('end',
                       label = "End Position",
                       value = 6769097,
                       min = 2
          )
        )
      )),
      column(6,
             plotOutput("signal_plot", height = 150)
      )
    ),
    hr(),
    macs2UI("macs2")
    
  )
  
  server <- function(input, output, session) {
    
    bams <- reactive({
      .subset_bams(trt_bam = trt_bam, ctrl_bam = ctrl_bam, 
                   chrom = input$chrom, start = input$start, end = input$end)
    })
    
    #locs <- reactiveValues(chrom = input$chrom, start = input$start, end = input$end)
    
    chrom <- reactive({input$chrom})
    start <- reactive({input$start})
    end <- reactive({input$end})
    
    output$signal_plot <- renderPlot(height = 150, {
      bs <- BSgenome.Hsapiens.UCSC.hg38
      
      bamz <- bams()
      
      gr <- GRanges(seqnames = input$chrom,
                    ranges = IRanges(start = c(input$start), end = c(input$end)))
      
      # Plotting code.
      trt.track <- autoplot(bamz$trt, which = gr, bsgenome = bs)
      
      if (!is.null(bamz$ctrl)) {
        ctrl.track <- autoplot(bamz$ctrl, which = gr, bsgenome = bs)
        
        tks <- tracks(Treat = trt.track, Control = ctrl.track, 
                      heights = c(0.5, 0.5)) + theme_clear()
      } else {
        ctrl.track <- NULL
        
        tks <- tracks(Treat = trt.track,  
                      heights = c(0.5)) + theme_clear()
      }
      
      tks
    })
    
    macs2Server("macs2", trt_bam = bams()$trt, ctrl_bam = bams()$ctrl, 
                chrom = reactive(input$chrom), start = reactive(input$start), end = reactive(input$end))
    
  }
  
  shinyApp(ui, server)
}