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
peakPeekeR <- function(trt_bam, ctrl_bam = NULL) {
  ui <- fluidPage(
    titlePanel("peakPeekeR"),
    
    column(6, wellPanel(
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
    )),
    
    column(6,
           plotOutput("signal_plot")
    )
  )
  
  server <- function(input, output, session) {
    
    output$signal_plot <- renderPlot({
      bs <- BSgenome.Hsapiens.UCSC.hg38
      
      bams <- .subset_bams(trt_bam = trt_bam, ctrl_bam = ctrl_bam, 
                              chrom = input$chrom, start = input$start, end = input$end)
      
      gr <- GRanges(seqnames = input$chrom,
                    ranges = IRanges(start = c(input$start), end = c(input$end)))
      
      # Plotting code.
      trt.track <- autoplot(bams$trt, which = gr, bsgenome = bs)
      
      if (!is.null(bams$ctrl)) {
        ctrl.track <- autoplot(bams$ctrl, which = gr, bsgenome = bs)
        
        tks <- tracks(Treat = trt.track, Control = ctrl.track, 
                      heights = c(0.5, 0.5)) + theme_tracks_sunset()
      } else {
        ctrl.track <- NULL
        
        tks <- tracks(Treat = trt.track,  
                      heights = c(0.5)) + theme_tracks_sunset()
      }
      
      tks
      
    })
    
  }
  
  shinyApp(ui, server)
}