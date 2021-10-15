library(shiny)
library(plotly)
library(ggplot2)
library(shinyjs)
library(DT)
library(Seurat)
library(CARD)
library(shinydashboard)
library(shinyalert)
library(shinyFiles)
library(reticulate)
library(shinyWidgets)

shiny_one_panel = fluidPage(
    titlePanel("Spatial analysis"),
    hr(),

    fluidRow(
        ##------Sidebar---------
        column(3,
               h4('Load Data:'),
               wellPanel(
                   fileInput(inputId = 'tpmFiles',
                             label = "Load spatial dataset",
                             multiple = FALSE,
                             accept = ".h5"),
                   fileInput(inputId = 'cellAnnoFiles',
                             label = "Metadata",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain")),
                   checkboxInput(inputId = 'norm',
                                 label = "Normalise?",
                                 value = TRUE),
                   fluidRow(
                       column(6,
                              numericInput(inputId = "min.genes",
                                           label = "Min. genes",
                                           value = 200,
                                           min = 1)
                       ),
                       column(6,
                              numericInput(inputId = "min.cells",
                                           label = "Min. cells",
                                           value = 3,
                                           min = 1)
                       )
                   ),
                   textInput(inputId = "projName",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton", "Load Data", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset", "Reset Data", icon = icon("repeat"))
                       )
                   )
               ),

               ##------Plot download---------
               h4("Export to PDF:"),
               wellPanel(
                   ## Conditional panel for different plots
                   conditionalPanel(" input.QC == 'QC_panel1' && input.tabs == 'QC plots' ",
                                    actionButton("PDFa", "Download Spatial Plot in PDF", icon = icon("download"))
                   ),
                   #conditionalPanel(" input.tabs == 'Variable Gene Plot' ",
                   #                 actionButton("PDFc", "Download Dimension Reduction Plot in PDF", icon = icon("download"))
                   #),
                   conditionalPanel(" input.Pca == 'P_panel1' && input.tabs == 'Dimension Reduction' ",
                                    actionButton("PDFd", "Download PCA Plot(PDF) and co-ordinates", icon = icon("download"))
                   ),
                   conditionalPanel(" input.Pca == 'P_panel2' && input.tabs == 'Dimension Reduction' ",
                                    actionButton("PDFe", "Download PC Gene plots and PC table", icon = icon("download"))
                   ),
                   conditionalPanel(" input.diag == 'D_panel' && input.tabs == 'Deconvolution' ",
                                    actionButton("PDFg", "Download Deconvolution in PDF", icon = icon("download"))
                   ),
                   #conditionalPanel(" input.diag == 'D_panel2' && input.tabs == 'Significant PC testing' ",
                   #                 actionButton("PDFh", "Download Elbow Plot in PDF", icon = icon("download"))
                   #),
                   #conditionalPanel(" input.tabs == 'UMAP' ",
                   #                 actionButton("PDFi", "Download UMAP Plot(PDF) and co-ordinates", icon = icon("download"))
                   #),
                   #conditionalPanel(" input.tabs == 'TSNE' ",
                   #                 actionButton("PDFj", "Download TSNE Plot(PDF) and co-ordinates", icon = icon("download"))
                   #),
                   conditionalPanel(" input.tabs == 'DEGs' ",
                                    actionButton("PDFk", "Download DEG results with plots", icon = icon("download"))
                   ),
                   conditionalPanel(" input.tabs == 'DEGs Heatmap' ",
                                    actionButton("PDFl", "Download DEG heatmap", icon = icon("download"))
                   ),
                   ## ensure no spill over in button text
                   tags$head(
                       tags$style(HTML('
                                   .btn {
                                   white-space: normal;
                                   }'
                       )
                       )
                   ),
                   ## Conditional is separate from pdf options
                   hr(),
                   fluidRow(
                       column(6,
                              sliderInput(inputId="pdf_w", label = "PDF width(in):",
                                          min=3, max=20, value=8, width=100, ticks=F)
                       ),
                       column(6,
                              sliderInput(inputId="pdf_h", label = "PDF height(in):",
                                          min=3, max=20, value=8, width=100, ticks=F)
                       )),

                   actionButton("OpenDir", "Open download folder", icon = icon("folder"))
               ),

               ##------Save Data---------
               hr(),
               actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),

               hr(),
               h4(tags$a(href="mailto:Chen_Jinmiao@immunol.a-star.edu.sg?subject=[cytof-question]",
                         "Contact Us")),
               imageOutput("logo", height = "60px")
        ),
        ##------Main area---------
        column(9,
               tabsetPanel(type = "pills", id = "tabs",
                           ## add file preview tab
                           ##------QC plots---------
                           tabPanel("QC plots", fluidPage(
                               hr(),
                               tabsetPanel(id="QC",
                                           tabPanel(title="Spatial QC Plots", value = "QC_panel1",
                                                    br(),
                                                    fluidRow(
                                                        column(10,
                                                               plotOutput("nCount_SpatialPlot", width = "50%"),
                                                               br(),
                                                               plotOutput("SpatialFeaturePlot", width = "100%"),
                                                               #plotlyOutput("nCount_RNAPlot", width = "100%")
                                                        ),
                                                        #column(6,
                                                        #       verbatimTextOutput("name")
                                                        #)
                                                    )

                                           )

                               )
                           )),

                           ##------Dimension Reduction---------
                           tabPanel("Dimension Reduction", fluidPage(
                               hr(),
                               tabsetPanel(id="Pca",
                                           tabPanel(title="Dimension Reduction Plot", value="P_panel1",
                                                    br(),
                                                    fluidRow(
                                                        column(3,
                                                               actionButton("doPCA", "Run", icon = icon("hand-pointer-o"))
                                                        ),

                                                        #   column(3,
                                                        #          selectInput("x.pc",
                                                        #                      label = "X-axis PC to use",
                                                        #                      choices = 1:20,
                                                        #                      selected = 1)
                                                        #placeholder for z-axis input
                                                        #   ),
                                                        #   column(3,
                                                        #          selectInput("y.pc",
                                                        #                      label = "Y-axis PC to use",
                                                        #                      choices = 1:20,
                                                        #                      selected = 2)
                                                        #   ),
                                                        #   column(3,
                                                        #          selectInput("z.pc",
                                                        #                      label = "Z-axis PC to use",
                                                        #                      choices = 1:20,
                                                        #                      selected = 3)
                                                        #   )
                                                        # ),
                                                        # fluidRow(
                                                        #     column(3 ,
                                                        #            numericInput("pc.plot.size",
                                                        #                         label = "Point Size:",
                                                        #                         value = 1,
                                                        #                         min = 0.1,
                                                        #                         step = 0.5)
                                                        #     ),
                                                        #     column(3,
                                                        #            sliderInput("pca.plot.alpha",
                                                        #                        label = "Point Transparency",
                                                        #                        min = 0,
                                                        #                        max = 1,
                                                        #                        step = 0.1,
                                                        #                        value = 0.8)
                                                        #     ),

                                                        column(3,
                                                               numericInput("clus.res",
                                                                            label = "Resolution used",
                                                                            value = 0.6,
                                                                            min = 0.1,
                                                                            step = 0.1)
                                                        ),
                                                    ),
                                                    br(),
                                                    plotlyOutput("DimPlot", width = "100%"),
                                                    br(),
                                                    plotOutput("SpatialDimPlot", width = "100%"),
                                                    br(),
                                                    fluidRow(
                                                        column(6,
                                                               uiOutput("deg1.gene.select"),
                                                               plotOutput("Deg1.plot", width = "100%"),
                                                               #plotlyOutput("Deg1.plot", width = "100%")
                                                        ),
                                                    ),
                                           ),
                                           tabPanel(title="PC Gene Visualisation", value="P_panel2",
                                                    br(),
                                                    selectInput("select.pc",
                                                                label = "PC to plot",
                                                                choices = c(1:20)
                                                    ),
                                                    fluidRow(
                                                        column(4,
                                                               plotOutput("vizPlot", width = "100%", height = "600px")
                                                        ),
                                                        column(8,
                                                               plotOutput("PCHeatmap", width = "100%", height = "600px")
                                                        )
                                                    ),
                                                    DT::dataTableOutput("PCtable")
                                           )
                               )
                           )),

                           ##------DEGs---------
                           tabPanel("DEGs", fluidPage(
                               hr(),
                               fluidRow(
                                   #   column(4,
                                   #          uiOutput("clust1")
                                   #   ),
                                   #   column(4,
                                   #          uiOutput("clust2")
                                   #   ),

                                   column(3, selectInput("min_pct",
                                                         label = "min.pct",
                                                         choices = c("0.1", "0.25"))
                                   ),

                                   column(3, selectInput("logfc",
                                                         label = "logfc.threshold",
                                                         choices = c("0.1", "0.25"))
                                   ),

                                   column(4,
                                          actionButton("doDeg", "Run DEGs", icon = icon("hand-pointer-o"))
                                   )),
                               fluidRow(
                                   column(6,
                                          uiOutput("deg.gene.select"),
                                          plotOutput("Deg.plot", width = "100%"),
                                   ),
                                   column(6,
                                          DT::dataTableOutput("Deg.table"),
                                   )),
                               br(),
                               fluidRow(column(4,
                                               actionButton("doDegn", "Find Spatially Variable Features", icon = icon("hand-pointer-o"))
                               )),
                               br(),
                               fluidRow(
                                   column(6,
                                          uiOutput("degn.gene.select"),
                                          plotOutput("Degn.plot", width = "100%"),
                                   )),
                           )),

                           ##------Deconvolution---------

                           tabPanel("Deconvolution using Spotlight", fluidPage(
                               hr(),
                               fluidRow(
                                   column(3,
                                          fileInput(inputId = 'tpmFiles1',
                                                    label = "scRNA-seq dataset",
                                                    multiple = FALSE,
                                                    accept = ".rds")),
                                   br(),
                                   column(3,
                                          actionButton("loadButton1", "Load scRNA-seq dataset", icon = icon("hand-o-right"),
                                          ),
                                   ),
                                   #br(),
                                   fluidRow(
                                       #selectInput("dim.used",
                                       #            label = "Dimensions",
                                       #            choices = c(1:50)
                                       #),
                                       #br(),
                                       #fluidRow(
                                       column(9,
                                              br(),
                                              actionButton("process_scRNA", "Process scRNA-seq dataset", icon = icon("hand-pointer-o")),
                                              br(),
                                              #column(6,
                                              br(),
                                              plotOutput("scRNAPlot", width = "100%")),
                                       br(),
                                       column(9,
                                              br(),
                                              actionButton("vis_spRNA", "Visualise spatial dataset", icon = icon("hand-pointer-o")),
                                              br(),
                                              br(),
                                              plotOutput("spRNAPlot", width = "100%")),
                                       br(),
                                       #),
                                       column(9,
                                              br(),
                                              actionButton("get_markers", "Get Marker Genes", icon = icon("hand-pointer-o")),
                                              br(),
                                              br(),
                                              actionButton("doDeconv", "Deconvolute", icon = icon("hand-pointer-o")),
                                              br(),
                                              br(),
                                              plotOutput("DeconvPlot", width = "100%")),
                                       br(),
                                   ),
                               )
                           ))
                           ##------END---------
               )
        )
    )
)

dbHeader <- dashboardHeader(title = "cytofkit2")
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/JinmiaoChenLab',
                                           tags$img(src='https://avatars1.githubusercontent.com/u/8896007?s=400&u=b0029c2e64f405ea0a46d311239b674a430ec77c&v=4'
                                                    ,height='60',width='60', align='left')
                                           , tags$div('cytofkit2', style='color:white;font-family:arial rounded MT bold'))

dashboardPage(skin = "yellow",
              dbHeader,
              dashboardSidebar(
                  sidebarMenu(id = "sbm",
                              menuItem(tags$p(style = "display:inline;font-size: 20px;", "Seurat"), tabName = "seurat", icon = icon('cog'))


                  )# end of sidebarMenu
              ),#end of dashboardSidebar
              dashboardBody(
                  #includeCSS("www/custom.css")
                  useShinyalert()
                  , shinyjs::useShinyjs()
                  , tabItem(
                      tabName = "seurat"
                      , shiny_one_panel
                  ) # End of tabItem

              )# end of dashboard body
)# end of dashboard page
