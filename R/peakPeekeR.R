#' peakPeekeR
#'
#' Does the peak peeking.
#' 
#' @param trt_bam Path for treatment sample BAM file.
#' @param ctrl_bam Optional path to control (input/Igg) BAM file.
#'
#' @return A Shiny app for peak calling parameter optimization.
#' @author Jared Andrews, Lawryn Kasper, Jacob Steele
#'
#' @export
#' @importFrom ggbio autoplot tracks theme_clear
#' @import shiny
#' @import shinyBS
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom utils read.delim
#' @importFrom ggplot2 aes_string
peakPeekeR <- function(trt_bam, ctrl_bam = NULL) {
  ui <- fluidPage(
    titlePanel("peakPeekeR"),
    tags$head(
      tags$style(HTML("
                      hr {border-top: 1px solid #000000;
                          margin: 5px 0px 15px 0px;}
                      h4 {margin: 0px;
                          padding: 0px;
                          font-weight: bold;}
                      "))
    ),
    fluidRow(
      column(6, wellPanel(
        h3("Subset Bams"),
        splitLayout(
          textInput("chrom",
                    label = "Chromosome",
                    value = "chr12"
          ),
          
          numericInput("start",
                    label = "Start Position",
                    value = 6522378,
                    min = 1
          ),
          
          numericInput("end",
                       label = "End Position",
                       value = 6769097,
                       min = 2
          )
        )
      )),
      column(6,
        wellPanel(h3("Plot Location"),
         splitLayout(
           textInput("plot.chrom",
                     label = "Chromosome",
                     value = "chr12"
           ),
           
           numericInput("plot.start",
                        label = "Start Position",
                        value = 6522378,
                        min = 1
           ),
           
           numericInput("plot.end",
                        label = "End Position",
                        value = 6769097,
                        min = 2
           )
         )
        )
      )
    ),
    hr(),
    actionButton("add", "Add Peak Caller Instance", icon = icon("plus"), class='btn-success'),
    macs2UI("macs2")
    
  )
  
  server <- function(input, output, session) {
    
    bams <- reactive({
      .subset_bams(trt_bam = trt_bam, ctrl_bam = ctrl_bam, 
                   chrom = input$chrom, start = input$start, end = input$end)
    })
    
    sig.tracks <- reactive({
      bs <- BSgenome.Hsapiens.UCSC.hg38
      
      bamz <- bams()
      
      gr <- GRanges(seqnames = input$plot.chrom,
                    ranges = IRanges(start = c(input$plot.start), end = c(input$plot.end)))
      
      # Plotting code.
      trt.track <- autoplot(bamz$trt, which = gr, bsgenome = bs)
      
      if (!is.null(bamz$ctrl)) {
        ctrl.track <- autoplot(bamz$ctrl, which = gr, bsgenome = bs)
      } else {
        ctrl.track <- NULL
      }
      
      list(trt = trt.track, ctrl = ctrl.track)
    })
    
    # observeEvent(input$add, {
    #   i <- sprintf('%04d', input$add)
    #   id <- sprintf('macs2%s', i)
    #   
    #   insertUI(
    #     selector = '#add',
    #     where = "beforeBegin",
    #     ui = macs2UI(id)
    #   )
    #   
    #   macs2Server(id, trt_bam = bams()$trt, ctrl_bam = bams()$ctrl, 
    #               chrom = reactive(input$plot.chrom), start = reactive(input$plot.start), 
    #               end = reactive(input$plot.end),
    #               trt_track = reactive(sig.tracks()$trt), ctrl_track = reactive(sig.tracks()$ctrl))
    #   
    #   observeEvent(input[[paste0(id, '-deleteButton')]], {
    #     removeUI(selector = sprintf('#%s', id))
    #     .remove_shiny_inputs(id, input)
    #   })
    # })
    
    # Return the UI for a modal dialog with caller input
    callerModal <- function() {
      modalDialog(
        selectInput("caller", "Choose peak caller",
                  choices = c("MACS2", "MACS", "SICER2", "Genrich")
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton("ok", "Add", class='btn-success')
        )
      )
    }
    
    # Show modal when button is clicked.
    observeEvent(input$add, {
      showModal(callerModal())
    })
    
    # When OK button is pressed, attempt to load the data set. If successful,
    # remove the modal. If not show another modal, but this time with a failure
    # message.
    observeEvent(input$ok, {
      removeModal()
      # Check that data object exists and is data frame.
      if (input$caller == "MACS2") {
        i <- sprintf('%04d', input$add)
        id <- sprintf('macs2%s', i)
        
        insertUI(
          selector = '#add',
          where = "beforeBegin",
          ui = macs2UI(id)
        )
        
        macs2Server(id, trt_bam = bams()$trt, ctrl_bam = bams()$ctrl, 
                    chrom = reactive(input$plot.chrom), start = reactive(input$plot.start), 
                    end = reactive(input$plot.end),
                    trt_track = reactive(sig.tracks()$trt), ctrl_track = reactive(sig.tracks()$ctrl))
        
        observeEvent(input[[paste0(id, '-deleteButton')]], {
          removeUI(selector = sprintf('#%s', id))
          .remove_shiny_inputs(id, input)
        })
      } 
    })
    
  }
  
  shinyApp(ui, server)
}