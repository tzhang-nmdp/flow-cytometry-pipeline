library(shiny)

shinyUI(fluidPage(
    titlePanel("Quality Control and Data Process for Flow Cytometry Data"),
    
    fluidRow(
        column(3,
               fileInput('fcsFiles', strong('Choose fcs file(s):'), multiple = TRUE, 
                         accept = c('text/fcs', '.fcs')),
               
               actionButton("goButton", "Submit!"),
               hr(),
               downloadButton('downloadMarkers', 'Download markers table'),
               br(),
               downloadButton('downloadFCS', 'Download new FCS files'),
               hr(),
               downloadButton('downloadreport', 'Download all plot reports'),
               hr(),
               
               ## sample limits: 50
               uiOutput("sample_select"),
               
               lapply(1:50, function(i) {
                   uiOutput(paste0('timeSlider', i))
               }),
               
               hr(),
               uiOutput("marker_select"),
               
               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ", 
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        column(9,
              # Timeline QC
               tabsetPanel(type = "pills",
                           
                           tabPanel("Timeline QC", fluidPage(
                               hr(),
                               uiOutput("tl_sample_choose"),
                               
                               hr(),
                               h4("Timeline plot:"),
                               htmlOutput("tlplot"),
                             
                               hr(),
                               h4("Expression table plot:"),
                               plotOutput("tabplot"),
                               
                               hr(),
                               fluidRow(
                                   column(4, offset = 1,
                                          numericInput("tl_binSize", "Bin size:", value = NA)
                                   ),
                                   column(4,
                                          numericInput("tl_varCut", "Variation cut:", value = 1)
                                   ) 
                               ),
                               textOutput("tl_text")
                           )),
                           
              # Compensation and transformation preprocess
              tabPanel("Compensation and transformation preprocess", 
                     fluidPage(
                     hr(),
                     fluidRow( width =12,
                     column(width = 4, 
                            box(width = 12, "Compensation matrix heatmap",plotOutput("comp_heatmap"))
                     ), 
                     column(width = 4,
                            box(width = 12, "Uncompensated marker X and Y",plotOutput("uncomp_plot"))
                     ), 
                     column(width = 4, offset = 0.5, 
                            box(width = 12, "Compensated marker X and Y",plotOutput("comp_plot"))
                     )
                     ),
                     fluidRow(
                     column(width = 4,
                            uiOutput("comp_marker_x_axis")
                     ),     
                     column(width = 4, offset = 1, 
                            uiOutput("comp_marker_y_axis")
                     )                                
                     ),
                     hr()
                     )),  

              # Singlet and viability QC            
              tabPanel("Singlet and viability QC",  fluidPage(
                            hr(),
                            uiOutput("score_sample_choose"),
                            
                            h4("QA score summary:"),
                            fluidRow( width =12,
                                   column(width = 6, 
                                          box(width = 12, "Singlet (Doublet removal)",plotOutput("scplot"))
                                   ), 
                                   column(width = 6, 
                                          box(width = 12, "viability (Dead cell removal)",plotOutput("vcplot"))
                                   )
                            ), 
                            fluidRow(
                            column(width = 2,
                                   textInput("x_gate_sc", "x-axis Singlet gate points", value = "0e0,1e1,1e7,2e7")
                            ),
                            column(width = 2,  
                                   textInput("y_gate_sc", "y-axis Singlet gate points", value = "0e0,1e1,1e7,2e7")
                            ),
                            column(width = 2, offset = 2,
                                   textInput('min_xy_axis_vc', 'Lower bound viability xy-axis:', value = "-1e5,0e0")
                            ),
                            column(width = 2, 
                                   textInput('max_xy_axis_vc', 'Upper bound viability xy-axis:', value = "1e4,4.2e6")
                            )                                
                            ) 
                            )),

              # Celluar subpopulation auto gating
              tabPanel("Celluar subpopulation auto gating", fluidPage(
                     hr(),
                     uiOutput("score_sample_choose"),
                     
                     h4("QA score summary:"),
                     fluidRow( width =12,
                     column(width = 4, 
                            box(width = 12, "Auto gating on marker X",plotOutput("sbplot1"))
                     ), 
                     column(width = 4, 
                            box(width = 12, "Auto gating on marker Y",plotOutput("sbplot2")) 
                     ), 
                     column(width = 4, 
                            box(width = 12, "Auto gating on marker X and Y",plotOutput("sbplot3"))
                     )
                     ), 
                     fluidRow(                             
                     column(width = 2, offset = 0.2, 
                            uiOutput("gate_type")
                     ),                                 
                     column(width = 2, offset = 0.2, 
                            uiOutput("marker_x_axis")
                     ),     
                     column(width = 2, offset = 0.2, 
                            uiOutput("marker_y_axis")
                     )
                     ),
                     fluidRow(
                     column(width = 2,
                            textInput("x_gate", "x-axis gate points", value = "0e0,1e1,1e7,2e7")
                     ),
                     column(width = 2, offset = 0.2, 
                            textInput("y_gate", "y-axis gate points", value = "0e0,1e1,1e7,2e7")
                     ),
                     column(width = 2, offset = 0.2, 
                            textInput("pre_gate", "pre-filter gate points", value = "1e2,1e8")
                     )
                     )
              )),
              # Advanced clustering
              tabPanel("Advanced Clustering",
                            fluidPage(
                            verticalLayout(
                                   hr(),
                                   h4("FlowSOM Clustering:"),
                                   plotOutput("clustplot1"),
                                   h4("SPADE Clustering:"),
                                   # uiOutput("pdfview"),
                                   plotOutput("clustplot2"),
                                   hr()
                            )
                            )
              )
        )
    )
)
))
