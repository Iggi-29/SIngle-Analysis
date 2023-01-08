#Library----
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyFiles)
library(plyr)
library(stats)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggrepel)
library(DoubletFinder)
library(dyno)
library(CellChat)
library(ComplexHeatmap)
library(RECODE)
####-----

#UI ----
ui <- fluidPage(theme = "",
                div(style = "padding: 1px 0px; width: '100px'",
                    titlePanel(
                      title = ""
                    )),
                navbarPage(
                  tags$style(HTML("hr {border-top: 1px solid #000000;}")),
                  title = div("Single Analysis",
                              style = "position: relative; top: 50%; transform: translateY(-50%);"),
                  tabPanel("About Single Analysis",
                           h1("Select file to analyze"),
                           fileInput("seurat_object", label = NULL,accept = ".rds"),
                           h1("About the app"),
                           p("This app, is a Rshiny developed app that integrates different analysis that can be performed at a single sample level on scRNA-Seq data. The app is part of Ignasi Jarne Sanz's master's final thesis in bioinformatics and biostatistics."),
                           p("This app has been uploaded in a GitHub repository, hoping that it will be useful for researchers or students that are not into informatics to do their analysis. It is important to highlight that this app has not been developed with the purpose of being used at a clinical level."),
                           h1("Approach"),
                           p("Once the data has been entered into the app using one of the scripts provided in the Git-Hub repository where the app is located, the following procedures can be performed:"),
                           tags$ul(
                             tags$li(HTML("Data is firstly loaded into the system as a Seurat Object, this can be created using R codes provided below.")),
                             tags$li(HTML("After being loaded, data quality graphs will be displayed, data filtering and preprocessing practices that must be applied to the data before the downstream analysis can be performed using the app by applying different algorithms that the user can choose.")),
                             tags$li(HTML("Once the data has been filtered, the app enables its visualization by displaying UMAP and TNSE clusterization graphs. Different plots representing the expression of genes in the distinct cell groups and co-expression graphs can also be displayed, marker genes for each cluster can be also be obtained. Trajectory analysis can be performed using two different packages. Finally, cell communication inference can also be performed."))
                           ),
                           h1("Contact"),
                           p("You may contact the author of this app at: ignasijs@gmail.com or ijarnes@uoc.edu"),
                           h1("Packages used in the analysis"),
                           p("In order that the app runs the analysis correctly, the following packages should be installed in your working environment:"),
                           tags$ul(
                             tags$li(HTML("shiny")),
                             tags$li(HTML("shinythemes")),
                             tags$li(HTML("shinydashboard")),
                             tags$li(HTML("shinyFiles")),
                             tags$li(HTML("plyr")),
                             tags$li(HTML("ggplot2")),
                             tags$li(HTML("ggrepel")),
                             tags$li(HTML("RColorBrewer")),
                             tags$li(HTML("dplyr")),
                             tags$li(HTML("stats")),
                             tags$li(HTML("Matrix")),
                             tags$li(HTML("Seurat")),
                             tags$li(HTML("DoubletFinder")),
                             tags$li(HTML("dyno")),
                             tags$li(HTML("ComplexHeatmap")),
                             tags$li(HTML("CellChat"))),
                           p("ATTENTION: Docker is also required in order to enable the trajectory analysis performed by dyno !"),
                           tags$p(HTML("Once docker is isntalled follow the instructions in the following link to its propper configuration: <a href=\"https://www.digitalocean.com/community/questions/how-to-fix-docker-got-permission-denied-while-trying-to-connect-to-the-docker-daemon-socket\">https://www.digitalocean.com/community/questions/how-to-fix-docker-got-permission-denied-while-trying-to-connect-to-the-docker-daemon-socket</a>.")),
                           h1("Git-Hub"),
                           tags$p(HTML("Source code is available at: <a href=\"https://github.com/Iggi-29/Single_Analysis\">https://github.com/Iggi-29/Single_Analysis</a>.")),
                           tags$p(HTML("Code to create Seurat objects from different file types is available at: <a href=\"https://github.com/Iggi-29/Create_RDS_files\">https://github.com/Iggi-29/Create_RDS_files</a>.")),
                           tags$p(HTML("Code to install all the required packages is provided at: <a href=\"https://github.com/Iggi-29/Install_Single_Analysis\">https://github.com/Iggi-29/Install_Single_Analysis</a>."))),
                  
                  #QUALITY AND PREPROCESS----
                  tabPanel("Data quality and filtering",
                           p("First, choose whether or not you want to use RECODE to perform data denoising."),
                           p("Second, get graphs regarding the quality of the data and get some number regarding this quality."),
                           p("Third, filter the data by choosing a maximum number of genes per cell and the maximum percentage of mitocondrial genes expressed in a cell."),
                           p(),
                           tabsetPanel(type = "tabs",
                                       tabPanel("Data Quality",
                                                sidebarPanel(
                                                  h4("Data Denoising and Expression Recovery"),
                                                  radioButtons("radiorecode", choices = list(`No` = "No", `Yes` = "Yes"), label = NULL),
                                                  actionButton("clickrecode", "Display quality graphs"),
                                                  textOutput("recodetext"),
                                                  hr(),
                                                  h4("Data quality numbers"),
                                                  div(style="display:inline-block", p("Number of cells:")),
                                                  div(style="display:inline-block", textOutput("raw_cnumber")),
                                                  p(),
                                                  div(style="display:inline-block", p("Median number of genes per cell:")),
                                                  div(style="display:inline-block", textOutput("raw_ngenes")),
                                                  p(),
                                                  div(style="display:inline-block", p("Median number of UMIs per cell:")),
                                                  div(style="display:inline-block", textOutput("raw_nUMI")),
                                                  p(),
                                                  div(style="display:inline-block", p("Percentage of mitochondrial genes per cell:")),
                                                  div(style="display:inline-block", textOutput("raw_percentmito")),
                                                  hr(),
                                                  h4("Recomended quality thresholds"),
                                                  div(style="display:inline-block", p("Minimum number of genes in a cell:")),
                                                  div(style="display:inline-block", p(200)),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum number of genes in a cell:")),
                                                  div(style="display:inline-block", textOutput("recommended_maximum_genes")),
                                                  p(),
                                                  div(style="display:inline-block", p("Minimum number of UMIs in a cell:")),
                                                  div(style="display:inline-block", p(1000)),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum number of UMIs in a cell:")),
                                                  div(style="display:inline-block", textOutput("recommended_maximum_umis")),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum percentage of mitochondrial genes:")),
                                                  div(style="display:inline-block", p(12)),
                                                  hr(),
                                                  h4("Quality thresholds"),
                                                  div(style="display:inline-block", p("Minimum number of genes in a cell:")),
                                                  div(style="display:inline-block", numericInput("ngene.low", 
                                                                                                 value = 200, 
                                                                                                 label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block;", p("Maximum number of genes in a cell:")),
                                                  div(style="display:inline-block;", numericInput("ngene.high", 
                                                                                                  value = 4000, 
                                                                                                  label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block", p("Minimum number of UMIs in a cell:")),
                                                  div(style="display:inline-block", numericInput("numi.low", value = 1000, label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum number of UMIs in a cell:")),
                                                  div(style="display:inline-block", numericInput("numi.high", 
                                                                                                 value = 10000, 
                                                                                                 label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum percentage of mitochondrial genes:")),
                                                  div(style="display:inline-block", numericInput("mito.high", 
                                                                                                 value = 12, 
                                                                                                 label = NULL, step = 1, max = 100, min = 0, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block;vertical-align:top;", p("Select the species where the samples come from:")),
                                                  div(style="display:inline-block;vertical-align:top;", selectInput("mitocondrial_pattern", label = NULL, width = 100,
                                                                                                                    choices = list(`Mouse` = ("^mt-"),
                                                                                                                                   `Human` = ("^MT-"))))),
                                                mainPanel(
                                                  fluidRow(
                                                    box(title = "Number of genes - percentage of mitochondrial genes", status = "success", solidHeader = T,
                                                        plotOutput("quality_graph_1")),
                                                    box(title = "Number of UMIs - Number of genes", status = "success", solidHeader = T,
                                                        plotOutput("quality_graph_2"))),
                                                  fluidRow(
                                                    box(title = "Violin plots", status = "success",solidHeader = T,
                                                        plotOutput("quality_graph_3")),
                                                    box(title = "Violin plots", status = "success", solidHeader = T,
                                                        plotOutput("quality_graph_4"))))),
                                       tabPanel("Data Filtering",
                                                sidebarPanel(
                                                  h4("Filtering criteria"),
                                                  div(style="display:inline-block",p("Minimum number of genes in a cell:")),
                                                  div(style="display:inline-block", numericInput("minimum.genes", 
                                                                                                 value = 200, 
                                                                                                 label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum number of genes in a cell:")),
                                                  div(style="display:inline-block", numericInput("maximum.umis", 
                                                                                                 value = 100000, 
                                                                                                 label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block", p("Maximum percentage of mitochondrial genes:")),
                                                  div(style="display:inline-block", numericInput("maximum.mito",
                                                                                                 value = 12, 
                                                                                                 label = NULL, width = 100, max = 100, min = 0)),
                                                  p(),
                                                  actionButton("click", "Filter the data"),
                                                  hr(),
                                                  h4("Data metrics after filtering"),
                                                  div(style="display:inline-block", p("Number of cells")),
                                                  div(style="display:inline-block", textOutput("ncells_filtered")),
                                                  p(),
                                                  div(style="display:inline-block", p("Median number of genes per cell:")),
                                                  div(style="display:inline-block", textOutput("ngenes_filtered")),
                                                  p(),
                                                  div(style="display:inline-block",p("Median number of UMIs per cell:")),
                                                  div(style="display:inline-block",textOutput("numis_filtered")),
                                                  p(),
                                                  div(style="display:inline-block",p("Percentage of mitochondrial genes per cell:")),
                                                  div(style="display:inline-block",textOutput("mitocondrial_genes_filtered"))),
                                                
                                                mainPanel(
                                                  fluidRow(
                                                    box(title = "Number of genes - percentage of mitochondrial genes", status = "success", solidHeader = T,
                                                        plotOutput("quality_graph_1_filtered")),
                                                    box(title = "Number of UMIs - Number of genes", status = "success", solidHeader = T,
                                                        plotOutput("quality_graph_2_filtered"))),
                                                  fluidRow(
                                                    box(title = "Violin plots", status = "success",solidHeader = T,
                                                        plotOutput("quality_graph_3_filtered")),
                                                    box(title = "Violin plots", status = "success", solidHeader = T,
                                                        plotOutput("quality_graph_4_filtered"))))))),
                  
                  #----
                  #DATA PREPROCESSING AND FILTERING----
                  tabPanel("Data preprocessing and Doublet removal",
                           p("First, preprocess the data by performing: Normalization, Variable features computation, Data scaling and Dimensionality reduction."),
                           p("Once the data has been preprocessed Doublet removal can be performed, this step can be this process can be time-consuming and resource-intensive and in the case of big datasets could fail."),
                           tabsetPanel(type = "tabs",
                                       tabPanel("Data Preprocessing",
                                                (sidebarPanel(width = 2,
                                                              h4("1. Normalization"),
                                                              p("Select normalization method"),
                                                              selectInput("normalization_method", label = NULL,
                                                                          choices = list(`LogNormalize` = ("LogNormalize"),
                                                                                         `CLR` = ("CLR"),
                                                                                         `RC` = ("RC"))),
                                                              p("Set margin:"),
                                                              selectInput("margin", label = NULL,
                                                                          choices = list(`Features` = ("1"),
                                                                                         `Cell` = ("2"))),
                                                              p("Select scaling factor:"),
                                                              selectInput("scale.factor",label = NULL,
                                                                          choices = list(`Default (10.000)` = ("10000"),
                                                                                         `1.000` = ("1000"),
                                                                                         `100.000` = ("100000"),
                                                                                         `1.000.000` = ("1000000"))),
                                                              actionButton("click2",label = "Normalize data"),
                                                              p(),
                                                              textOutput("soc.normal"))),
                                                sidebarPanel(
                                                  h4("2. Highly variable features discovery"),
                                                  p("Select method to choose top variable features"),
                                                  selectInput("high_method", label = NULL,
                                                              choices = list(`vst` = ("vst"),
                                                                             `mean.var.plot` = ("mvp"),
                                                                             `dispersion` = ("disp"))),
                                                  p("Set the number of genes to be select as top:"),
                                                  numericInput("number_of_features", label = NULL, value = 2000),
                                                  p("Set the number of bins to be used:"),
                                                  numericInput("number_of_bins",label = NULL, value = 20),
                                                  actionButton("click3",label = "Find variable features"),
                                                  p(),
                                                  textOutput("soc.normal2"),
                                                  p(),
                                                  plotOutput("plot_find_variable_features")),
                                                sidebarPanel(width = 2,
                                                             h4("3. Data scaling"),
                                                             p("Scale the data:"),
                                                             actionButton("click4","Scale the data"),
                                                             p(),
                                                             textOutput("soc.normal3")),
                                                sidebarPanel(
                                                  h4("4. Dimensionality reduction"),
                                                  p("Perform PCA analysis"),
                                                  actionButton("click5","Run PCA"),
                                                  p(),
                                                  textOutput("soc.normal4"),
                                                  p(),
                                                  plotOutput("elbowplot"),
                                                )),
                                       tabPanel("Doublet Removal",
                                                br(),
                                                p("In order to perform this analysis, DoubletFinder will be used, its algorithm needs two parameters to be tunned, (pK and the amount of doublets expected to be in our dataset)."),
                                                p("The pK value is computed after clusterization, its calculation can take more than 5 minutes in large datasets."),
                                                tags$p(HTML("In order to adjust the amount of doublets in the raw dataset, please use the following proportions as a standard <a href=\"https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x\">https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x</a>. This table is provided by 10xGenomics.")),
                                                p("Once the pK value and the adjusted amount of doublets have been computed, enter the results of their respective sections in order to remove the doublets."),
                                                tags$p(HTML("For more information about DoubletFinder visit: <a href=\"https://github.com/chris-mcginnis-ucsf/DoubletFinder\">https://github.com/chris-mcginnis-ucsf/DoubletFinder</a>.")),
                                                br(),
                                                sidebarPanel(
                                                  h4("1. pK value optimization"),
                                                  p("Run preprocessing steps:"),
                                                  p("Set number of PCs to be used:"),
                                                  numericInput("pcstouse",value = 1, label = NULL),
                                                  actionButton("click6","Cluster the data"),
                                                  p(),
                                                  textOutput("soc.normal5"),
                                                  p(),
                                                  actionButton("click7","Calculate optimus pK value"),
                                                  p(),
                                                  textOutput("soc.normal6")),
                                                sidebarPanel(
                                                  h4("2. Adjust the number of expected doublets"),
                                                  p("Enter the expected proportion of doublets"),
                                                  numericInput("expected_doublets", value = 0, label = NULL),
                                                  actionButton("click8","Adjust"),
                                                  p(),
                                                  textOutput("expectdoublets_adjusted")),
                                                sidebarPanel(
                                                  h4("3. Remove Doublets"),
                                                  p("Enter the optimus pK value:"),
                                                  numericInput("pk_optim", value = 0, label = NULL),
                                                  p("Enter the number of expected doublets"),
                                                  numericInput("expect_optim", value = 0, label = NULL),
                                                  actionButton("click9","Apply doublet finder"),
                                                  p(),
                                                  tableOutput("summary_the_doublet_action"),
                                                  p(),
                                                  p("Doublet Removal"),
                                                  actionButton("click10","Remove doublets"),
                                                  textOutput("dims_final")),
                                                mainPanel(plotOutput("optimus_pk"))),
                           )),
                  tabPanel("Clusterization",
                           tabsetPanel(type = "tabs",
                                       tabPanel("UMAP and TSNE clusterization",
                                                br(),
                                                p("First, check whether if doublet removal has been performed, if not, enter the number of dimensions to reduce the dataset. In order that the following steps work correctly, display UMAP and TSNE analysis."),
                                                sidebarPanel(
                                                  h4("Have you performed doublet removal?"),
                                                  radioButtons("radio", choices = list(`Yes` = "Yes", `No` = "No"), label = NULL),
                                                  p("Set number of PCs to be used:"),
                                                  numericInput("pcstouse_noclist",value = 1, label = NULL),
                                                  hr(),
                                                  h4("Clusterization graphs"),
                                                  div(style="display:inline-block;vertical-align:top;", p("Resolution to display:")),
                                                  div(style="display:inline-block;vertical-align:top;",
                                                      selectInput("resolution_to_use", label = NULL,width = 100,
                                                                  choices = list(`0.3` = "snn_res.0.3",
                                                                                 `0.5` = "snn_res.0.5",
                                                                                 `0.7` = "snn_res.0.7",
                                                                                 `1` = "snn_res.1",
                                                                                 `1.5` = "snn_res.1.5",
                                                                                 `2` = "snn_res.2"))),
                                                  br(),
                                                  div(style="display:inline-block;vertical-align:top;",p("Display graphs:")),
                                                  div(style="display:inline-block;vertical-align:top;",actionButton("click11","UMAP")),
                                                  div(style="display:inline-block;vertical-align:top;",actionButton("click12","TSNE")),
                                                  hr(),
                                                  h4("Save Seurat Object"),
                                                  downloadButton("downloadData", "Save"),
                                                  br()),
                                                mainPanel(
                                                  splitLayout(
                                                    cellWidths = c("50%","50%"),
                                                    box(title = "UMAP", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("plot_umap", width = "90%")),
                                                    box(title = "TSNE", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("plot_tsne", width = "90%"))))),
                                       tabPanel("Expression markers",
                                                br(),
                                                p("Enter the parameters to set up the search for marker genes."),
                                                p("After looking for marker genes, if desired, plot the top markers for each cluster."),
                                                sidebarPanel(
                                                  h4("Gene markers search options"),
                                                  p("FindAllMarkers options:"),
                                                  div(style="display:inline-block;vertical-align:top;",p("LogFC threshold:")),
                                                  div(style="display:inline-block;vertical-align:top;",numericInput("logf_thres", value = 2, label = NULL, width = 100)),
                                                  p(),
                                                  div(style="display:inline-block;vertical-align:top;",p("Pct threshold:")),
                                                  div(style="display:inline-block;vertical-align:top;",numericInput("pct_thres", value = 0.25, label = NULL, width = 100)),
                                                  p(),
                                                  actionButton("click19","Look for markers"),
                                                  p(),
                                                  textOutput("soc.normal7"),
                                                  hr(),
                                                  h4("Gene markers summary"),
                                                  div(style="display:inline-block;vertical-align:top;",p("Maximum number of genes per cluster:")),
                                                  div(style="display:inline-block;vertical-align:top;",numericInput("max_genes", value = 2, label = NULL, width = 100)),
                                                  div(style="display:inline-block;vertical-align:top;",actionButton("click20","Display the data")),
                                                  p(),
                                                  hr(),
                                                  h4("Compute graphs"),
                                                  actionButton("click21","Plot the data")
                                                ),
                                                mainPanel(
                                                  box(title = "Top marker genes", status = "success", solidHeader = T,
                                                      dataTableOutput("table_marked"), style = "height:500px; width:600px; overflow-y: scroll;overflow-x: scroll;"),
                                                  box(title = "Cluster annotation", status = "success", solidHeader = T,
                                                      plotOutput("umap_marked")))))),
                  tabPanel("Gene visualizations",
                           tabsetPanel(type = "tabs",
                                       tabPanel("Gene expression",
                                                br(),
                                                p("Enter a gene symbol and display different graphs in order to visualize its expression through the different groups."),
                                                p("NOTE: In the case of the heatmap plot, more than a gene can be specified."),
                                                fluidRow(
                                                  splitLayout(cellWidths = c("50%","50%"),
                                                              box(title = "Ridge plot", status ="success", solidHeader = T, width = "45%",
                                                                  p("Enter Gene Symbol:"),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      textInput("ridgevalue",label = NULL, width = 100)),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      actionButton("click13","Plot the data")),
                                                                  br(),
                                                                  plotOutput("ridge", width = "90%")),
                                                              box(title = "Violin plot", status = "success",solidHeader = T, width = "45%",
                                                                  p("Enter Gene Symbol:"),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      textInput("violinvalue",label = NULL, width = 100)),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      actionButton("click14","Plot the data")),
                                                                  br(),
                                                                  plotOutput("violin_plot", width = "90%")))),
                                                fluidRow(
                                                  splitLayout(cellWidths = c("50%","50%"),
                                                              box(title = "Feature Plot", status = "success",solidHeader = T, width = "45%",
                                                                  p("Enter Gene Symbol:"),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      textInput("featinput",label = NULL, width = 100)),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      actionButton("click15","Plot the data")),
                                                                  br(),
                                                                  plotOutput("featuring", width = "90%")),
                                                              box(title = "Heatmap", status = "success",solidHeader = T, width = "45%",
                                                                  
                                                                  p("Enter Gene Symbol:"),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      selectizeInput("heatinput",width = 200, choices=NULL, label = NULL,
                                                                                     multiple = TRUE, options = list(create = TRUE))),
                                                                  div(style="display:inline-block;vertical-align:top;",
                                                                      actionButton("click16","Plot the data")),
                                                                  br(),
                                                                  plotOutput("heatmap", width = "90%"))))),
                                       tabPanel("Co-expression",
                                                br(),
                                                p("Enter more than a gene symbol and display graphs regarding its combined expression."),
                                                p("NOTE: In the case of the 'Co-expression UMAP' section two genes can be specified while in the 'Dotplot' section a list of more than two can be entered."),
                                                fluidRow(
                                                  splitLayout(
                                                    cellWidths = c("50%","50%"),
                                                    box(title = "Co-expression UMAP", status = "success", solidHeader = T, width = "45%",
                                                        p("Enter two Gene Symbols:"),
                                                        div(style="display:inline-block;vertical-align:top;",
                                                            textInput("coexp1",label = NULL, width = 100)),
                                                        div(style="display:inline-block;vertical-align:top;",
                                                            textInput("coexp2",label = NULL, width = 100)),
                                                        div(style="display:inline-block;vertical-align:top;",
                                                            actionButton("click17","Plot the data")),
                                                        br(),
                                                        plotOutput("coexpressionplot1")),
                                                    box(title = "Dotplot", status = "success", solidHeader = T, width = "45%",
                                                        p("Enter a list of Gene Symbols:"),
                                                        div(style="display:inline-block;vertical-align:top;",
                                                            selectizeInput("coexp3",width = 200, choices=NULL, label = NULL,
                                                                           multiple = TRUE, options = list(create = T))),
                                                        div(style="display:inline-block;vertical-align:top;",
                                                            actionButton("click18", "Plot the data")),
                                                        br(),
                                                        plotOutput("coexpressionplot2"))))))),
                  tabPanel("Trajectory analysis",
                                                p("Use dyno to compute cell trajectory inference. Select one of the methods to infer the trajectory that dyno has available, all methods are listed in the link below."),
                                                tags$p(HTML("All methods available for dyno can be found at: <a href=\"https://dynverse.org/reference/dynmethods/method/\">https://dynverse.org/reference/dynmethods/method/</a>.")),
                                                sidebarPanel(
                                                  h4("Run dyno for trajectory analysis"),
                                                  p("Preprocess and plot the data"),
                                                  div(style="display:inline-block;vertical-align:top;", p("Infer trajectory using:")),
                                                  div(style="display:inline-block;vertical-align:top;",textInput("select_dyno_method", label = NULL, value = "slingshot")),
                                                  actionButton("click23","Execute dyno"),
                                                  p(),
                                                  p(),
                                                  div(style="display:inline-block;vertical-align:top;",p("Specify a gene to plot")),
                                                  div(style="display:inline-block;vertical-align:top;",textInput("select_dyno_gene", label = NULL, value = NULL)),
                                                  actionButton("click24","Plot the data"),
                                                  hr(),
                                                  p("Plot dyno heatmap"),
                                                  div(style="display:inline-block;vertical-align:top;", p("Number of genes")),
                                                  div(style="display:inline-block;vertical-align:top;", numericInput("features_dyno_heat", label = NULL, value = 20)),
                                                  actionButton("click25","Plot the heatmap")
                                                ),
                                                mainPanel(
                                                  splitLayout(
                                                    cellWidths = c("50%","50%"),
                                                    box(title = "Dyno trajectory", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("dyno")),
                                                    box(title = "Dyno trajectory", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("dyno_feat"))),
                                                  box(title = "Dyno heatmap", status = "success", solidHeader = T,
                                                      plotOutput("dyno_heat")))),
                  
                  tabPanel("Cell communication",
                           tabsetPanel(type = "tabs",
                                       tabPanel("CellChat",
                                                br(),
                                                p("Use CellChat to perform cell-comunication inference, once the minimum number of cells that are taken into account when considering whether a cellular communication network is active."),
                                                sidebarPanel(
                                                  h4("Run CellChat"),
                                                  div(style="display:inline-block",p("Minimum number of cells")),
                                                  div(style="display:inline-block",numericInput("min_cell_comunication",value = 10, label = NULL, width = 100)),
                                                  actionButton("click26","Run CellChat")),
                                                mainPanel(
                                                  splitLayout(
                                                    cellWidths = c("50%","50%"),
                                                    box(title = "Number of interactions", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("twoplotscells_1")),
                                                    box(title = "Interaction weights/strength", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("twoplotscells_2"))))),
                                       tabPanel("Signal pathway analysis",
                                                br(),
                                                p("Display graphs regarding any of the encountered cellular communication networks that CellChat has found to be expressed."),
                                                sidebarPanel(
                                                  h4("Display Network"),
                                                  p("Specify network to display:"),
                                                  div(style="display:inline-block",textInput("networkdisplay",label = NULL,value = NULL)),
                                                  div(style="display:inline-block",actionButton("click27","Plot the data")),
                                                ),
                                                mainPanel(
                                                  h3("Detected networks"),
                                                  p("The following networks have been detected and can be displayed:"),
                                                  textOutput("networks"),
                                                  splitLayout(
                                                    cellWidths = c("50%","50%"),
                                                    box(title = "Circular", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("circular")),
                                                    box(title = "Chord", status = "success", solidHeader = T, width = "45%",
                                                        plotOutput("chord")))))))
                  
                )#navbarPage
)#FluidPage
#----

# Server ----
server <- function(input, output, session){
  
  options(shiny.maxRequestSize=30000000000000000*1024^2)
  
  #Observing----
  observe({
    updateSelectInput(session, "mitocondrial_pattern",
                      choices = list(`Mouse` = ("^mt-"),
                                     `Human` = ("^MT-")))
  })
  observe({
    updateNumericInput(session,"ngene.low")
  })
  observe({
    updateNumericInput(session,"ngene.high")
  })
  observe({
    updateNumericInput(session,"mito.high")
  })
  observe({
    updateNumericInput(session,"numi.low")
  })
  observe({
    updateNumericInput(session,"numi.high")
  })
  
  ##Data quality after filter
  observe({
    updateNumericInput(session,"minimum.genes")
  })
  observe({
    updateNumericInput(session,"maximum.umis")
  })
  observe({
    updateNumericInput(session,"maximum.mito")
  })
  observe({
    updateSelectInput(session, "normalization_method",
                      choices = list(`LogNormalize` = ("LogNormalize"),
                                     `CLR` = ("CLR"),
                                     `RC` = ("RC")))
  })
  observe({
    updateSelectInput(session, "margin",
                      choices = list(`Features` = (1),
                                     `Cell` = (2)))
  })
  observe({
    updateSelectInput(session, "scale.factor",
                      choices = list(`Default (10.000)` = (10000),
                                     `1000`= 1000,
                                     `100.000` = (100000),
                                     `1.000.000` = (1000000)))
  })
  observe({
    updateSelectInput(session, "high_method",
                      choices = list(`vst` = ("vst"),
                                     `mean.var.plot` = ("mvp"),
                                     `dispersion` = ("disp")))
  })
  observe({
    updateNumericInput(session,"number_of_features")
  })
  observe({
    updateNumericInput(session,"number_of_bins")
  })
  observe({
    updateNumericInput(session,"pcstouse")
  })
  observe({
    updateNumericInput(session,"expected_doublets")
  })
  observe({
    updateNumericInput(session, "pk_optim")
  })
  observe({
    updateNumericInput(session,"expect_optim")
  })
  observe({
    updateSelectInput(session, "resolution_to_use",
                      choices = list(`0.3` = "snn_res.0.3",
                                     `0.5` = "snn_res.0.5",
                                     `0.7` = "snn_res.0.7",
                                     `1` = "snn_res.1",
                                     `1.5` = "snn_res.1.5",
                                     `2` = "snn_res.2"))
  })
  observe({
    updateTextInput(session, "ridgevalue")
  })
  observe({
    updateTextInput(session, "violinvalue")
  })
  observe({
    updateTextInput(session, "featinput")
  })
  observe({
    updateSelectizeInput(session, "heatinput")
  })
  observe({
    updateSelectizeInput(session, "coexp3")
  })
  observe({
    updateTextInput(session, "coexp1")
  })
  observe({
    updateTextInput(session, "coexp2")
  })
  observe({
    updateNumericInput(session,"logf_thres")
  })
  observe({
    updateNumericInput(session,"pct_thres")
  })
  observe({
    updateNumericInput(session,"max_genes")
  })
  observe({
    updateTextInput(session, "select_dyno_method")
  })
  observe({
    updateTextInput(session, "select_dyno_gene")
  })
  observe({
    updateNumericInput(session,"features_dyno_heat")
  })
  observe({
    updateRadioButtons(session, "radio",
                       choices = list(`Yes` = "Yes", `No` = "No"))
  })
  observe({
    updateNumericInput(session, "min_cell_comunication")
  })
  observe({
    updateTextInput(session, "networkdisplay")
  })
  observe({
    updateRadioButtons(session, "radiorecode",
                       choices = list(`No` = "No", `Yes` = "Yes"))
  })
  #----
  
  #Data loading----
  seuobjjj <- eventReactive(input$seurat_object, {
    readRDS(input$seurat_object$datapath)
  })
  #----
  
  #Data denoise----
  seuobjj <- eventReactive(c(input$radiorecode, input$clickrecode),ignoreInit = T,{
    
    if ((input$radiorecode == "No") && (is.null(input$clickrecode )== FALSE)){
      return(seuobjjj())}
    
     if (input$radiorecode == "Yes" && is.null(input$clickrecode)== FALSE){
       seuobjjj <- seuobjjj()
       data <- as.matrix(seuobjjj[["RNA"]]@counts)
       data_RECODED <- RECODE(data)
       seuobjjj[["RECODED"]] <- CreateAssayObject(Matrix(data_RECODED, sparse = T))
       DefaultAssay(seuobjjj) <- "RECODED"
       return(seuobjjj)}
  })
  
  output$recodetext <- renderText({
    if (dim(seuobjj()) > 1 && input$radiorecode == "Yes"){ 
      return("RECODE has run successfuly")}
    if (dim(seuobjj()) > 1 && input$radiorecode == "No"){
      return("Data has not been denoised")
    }
  })
  #----
  
  
  #Recomended thresholds----
  output$recommended_maximum_genes <- renderText(round(mean(seuobjj()@meta.data$nFeature_RNA)*4))
  output$recommended_maximum_umis <- renderText(max(seuobjj()@meta.data$nCount_RNA)-10000)
  #----
  
  
  #Raw data numbers----
  output$raw_cnumber <- renderText(dim(seuobjj())[2])
  output$raw_ngenes <- renderText(median(seuobjj()@meta.data$nFeature_RNA))
  output$raw_nUMI <- renderText(median(seuobjj()@meta.data$nCount_RNA))
  output$raw_percentmito <- renderText({
    seuobjj <- seuobjj()
    seuobjj[["percent.mt"]] <- PercentageFeatureSet(seuobjj, pattern = input$mitocondrial_pattern)
    mito_mean <- mean(seuobjj@meta.data$percent.mt)
  })
  
  #Display Graph1 RAW
  output$quality_graph_1 <- renderPlot({
    seuobjj <- seuobjj()
    seuobjj[["percent.mt"]] <- PercentageFeatureSet(seuobjj, pattern = input$mitocondrial_pattern)
    ngene_low <- input$ngene.low
    ngene_high<- input$ngene.high
    mito_high <- input$mito.high
    numi_low <- input$numi.low
    numi_high <- input$numi.high
    meta_data <- seuobjj@meta.data
    meta_data$Quality <- rep("Good", nrow(meta_data))
    meta_data <- mutate(meta_data, Quality =
                          ifelse(nFeature_RNA < ngene_low | nFeature_RNA > ngene_high,"Bad",
                                 ifelse(nCount_RNA < numi_low | nCount_RNA > numi_high,"Bad",
                                        ifelse(percent.mt > mito_high,"Bad","Good"))))
    
    meta_data$Quality <- factor(meta_data$Quality)
    ggplot(meta_data, aes(x = nFeature_RNA, y = percent.mt, color = Quality))+
      geom_point(size = 1, 
                 color = case_when(
                   meta_data$Quality == "Good" ~ "#00B159",
                   meta_data$Quality == "Bad" ~ "#ff2d00",
                 ))+
      xlab("Number of genes")+
      ylab("Percentage of mitochondrial genes")+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks = element_line(colour = "black",size = 0.5),
            axis.text = element_text(size = 12.5),
            axis.title = element_text(size = 12.5))+
      guides(color = guide_legend(override.aes = list(size = 5)))  
    })
  
  #Display Graph2 RAW
  output$quality_graph_2 <- renderPlot({
    seuobjj <- seuobjj()
    seuobjj[["percent.mt"]] <- PercentageFeatureSet(seuobjj, pattern = input$mitocondrial_pattern)
    ngene_low <- input$ngene.low
    ngene_high<- input$ngene.high
    mito_high <- input$mito.high
    numi_low <- input$numi.low
    numi_high <- input$numi.high
    meta_data <- seuobjj@meta.data
    meta_data$Quality <- rep("Good", nrow(meta_data))
    meta_data <- mutate(meta_data, Quality =
                          ifelse(nFeature_RNA < ngene_low | nFeature_RNA > ngene_high,"Bad",
                                 ifelse(nCount_RNA < numi_low | nCount_RNA > numi_high,"Bad",
                                        ifelse(percent.mt > mito_high,"Bad","Good"))))
    ggplot(meta_data, aes(x = nCount_RNA, y = nFeature_RNA, color = Quality))+
      geom_point(size = 1, 
                 color = case_when(
                   meta_data$Quality == "Good" ~ "#00B159",
                   meta_data$Quality == "Bad" ~ "#ff2d00",
                 ))+
      xlab("Number of UMIs")+
      ylab("Number of genes")+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks = element_line(colour = "black",size = 0.5),
            axis.text = element_text(size = 12.5),
            axis.title = element_text(size = 12.5))+
      guides(color = guide_legend(override.aes = list(size = 5)))
  })
  
  #Display Graph3 RAW
  output$quality_graph_3 <- renderPlot({
    seuobjj <- seuobjj()
    seuobjj[["percent.mt"]] <- PercentageFeatureSet(seuobjj, pattern = input$mitocondrial_pattern)
    VlnPlot(seuobjj, features = c("percent.mt","nCount_RNA","nFeature_RNA" ), pt.size = 0) 
  })
  output$quality_graph_4 <- renderPlot({
    seuobjj <- seuobjj()
    seuobjj[["percent.mt"]] <- PercentageFeatureSet(seuobjj, pattern = input$mitocondrial_pattern)
    VlnPlot(seuobjj, features = c("percent.mt","nCount_RNA","nFeature_RNA" ), pt.size = 1) 
  })
  #----
  
  #Quality filter----
  #Subset
  seuobjfilt <- reactive({
    seuobjj <- seuobjj()
    seuobjj[["percent.mt"]] <- PercentageFeatureSet(seuobjj, pattern = input$mitocondrial_pattern)
    seuobj <- subset(seuobjj, subset = nFeature_RNA > input$minimum.genes &
                       nFeature_RNA < input$maximum.umis &
                       percent.mt < input$maximum.mito)
    return(seuobj)
    
  })
  
  #Display filtered data numbers
  output$ncells_filtered <- renderText({input$click
    isolate(dim(seuobjfilt())[2])})
  output$ngenes_filtered <- renderText({input$click
    isolate(median(seuobjfilt()@meta.data$nCount_RNA))})
  output$numis_filtered <- renderText({input$click
    isolate(median(seuobjfilt()@meta.data$nCount_RNA))})
  output$mitocondrial_genes_filtered <- renderText({input$click
    isolate(mean(seuobjfilt()@meta.data$percent.mt))})
  
  #Display
  #Display Graph1 filtered
  output$quality_graph_1_filtered <- renderPlot({
    ngene_low <- input$ngene.low
    ngene_high<- input$ngene.high
    mito_high <- input$mito.high
    numi_low <- input$numi.low
    numi_high <- input$numi.high
    meta_data <- seuobjfilt()@meta.data
    meta_data$Quality <- rep("Good", nrow(meta_data))
    meta_data <- mutate(meta_data, Quality =
                          ifelse(nFeature_RNA < ngene_low | nFeature_RNA > ngene_high,"Bad",
                                 ifelse(nCount_RNA < numi_low | nCount_RNA > numi_high,"Bad",
                                        ifelse(percent.mt > mito_high,"Bad","Good"))))
    
    meta_data$Quality <- factor(meta_data$Quality)
    ggplot(meta_data, aes(x = nFeature_RNA, y = percent.mt, color = Quality))+
      geom_point(size = 1, 
                 color = case_when(
                   meta_data$Quality == "Good" ~ "#00B159",
                   meta_data$Quality == "Bad" ~ "#ff2d00",
                 ))+
      xlab("Number of genes")+
      ylab("Percentage of mitochondrial genes")+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks = element_line(colour = "black",size = 0.5),
            axis.text = element_text(size = 12.5),
            axis.title = element_text(size = 12.5))+
      guides(color = guide_legend(override.aes = list(size = 5)))
  })
  
  #Display Graph2 filtered
  output$quality_graph_2_filtered <- renderPlot({
    ngene_low <- input$ngene.low
    ngene_high<- input$ngene.high
    mito_high <- input$mito.high
    numi_low <- input$numi.low
    numi_high <- input$numi.high
    meta_data <- seuobjfilt()@meta.data
    meta_data$Quality <- rep("Good", nrow(meta_data))
    meta_data <- mutate(meta_data, Quality =
                          ifelse(nFeature_RNA < ngene_low | nFeature_RNA > ngene_high,"Bad",
                                 ifelse(nCount_RNA < numi_low | nCount_RNA > numi_high,"Bad",
                                        ifelse(percent.mt > mito_high,"Bad","Good"))))
    
    ggplot(meta_data, aes(x = nCount_RNA, y = nFeature_RNA, color = Quality))+
      geom_point(size = 1, 
                 color = case_when(
                   meta_data$Quality == "Good" ~ "#00B159",
                   meta_data$Quality == "Bad" ~ "#ff2d00",
                 ))+
      xlab("Number of UMIs")+
      ylab("Number of genes")+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks = element_line(colour = "black",size = 0.5),
            axis.text = element_text(size = 12.5),
            axis.title = element_text(size = 12.5))+
      guides(color = guide_legend(override.aes = list(size = 5)))
  })
  
  #Display Graph3 filtered
  output$quality_graph_3_filtered <- renderPlot({
    VlnPlot(seuobjfilt(), features = c("percent.mt","nCount_RNA","nFeature_RNA"), pt.size = 0) 
  })
  output$quality_graph_4_filtered <- renderPlot({
    VlnPlot(seuobjfilt(), features = c("percent.mt","nCount_RNA","nFeature_RNA" ), pt.size = 1) 
  })
  
  #----
  
  #Preprocess----
  #Normalization
  seuobjfilt_norm <- eventReactive(input$click2,{
    seuobj <-NormalizeData(seuobjfilt(),
                           normalization.method = input$normalization_method,
                           margin = as.numeric(as.character(input$margin)),
                           scale.factor = as.numeric(as.character(input$scale.factor))) 
    return(seuobj)
  })
  output$soc.normal <- renderText({
    if (dim(seuobjfilt_norm()) > 1) 
      "The data has been normalized"
  })
  
  #Find highly variable features
  seuobjfilt_norm_high <- eventReactive(input$click3, {
    seuobj <-FindVariableFeatures(seuobjfilt_norm(),
                                  selection.method = input$high_method,
                                  nfeatures = input$number_of_features,
                                  num.bin = input$number_of_bins) 
    return(seuobj)
  })
  output$soc.normal2 <- renderText({
    if (dim(seuobjfilt_norm_high()) > 1) 
      "Higly variable features are displayed below"
  })
  
  #Plot variable features
  output$plot_find_variable_features <-renderPlot({
    top10 <- head(VariableFeatures(seuobjfilt_norm_high()),10)
    plot <- VariableFeaturePlot(seuobjfilt_norm_high())
    LabelPoints(plot = plot, points = top10, repel = T, xnudge = 0, ynudge = 0)
  })
  
  #Scale the data
  seuobjfilt_norm_high_scale <- eventReactive(input$click4, {
    all_the_genes <- rownames(seuobjfilt_norm_high())
    seuobj <- ScaleData(seuobjfilt_norm_high(), 
                        features = all_the_genes) 
    return(seuobj)
  })
  
  output$soc.normal3 <- renderText({
    if (dim(seuobjfilt_norm_high_scale()) > 1) 
      "The data has been scaled"
  })
  
  #Dimensonallity reduction
  seuobjfilt_norm_high_scale_dimred <- eventReactive(input$click5, {
    seuobj <- RunPCA(seuobjfilt_norm_high_scale(), 
                     features = VariableFeatures(seuobjfilt_norm_high_scale()))
    return(seuobj)
  })
  output$soc.normal4 <- renderText({
    if (dim(seuobjfilt_norm_high_scale_dimred()) > 0) 
      "PCs information is displayed below:" 
  })
  
  #Elbow plot
  output$elbowplot <- renderPlot({
    pct <- seuobjfilt_norm_high_scale_dimred()[["pca"]]@stdev /sum(seuobjfilt_norm_high_scale_dimred()[["pca"]]@stdev) *100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct <5)[1]
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    pcs <- min(co1, co2)
    pcschar <- as.character(as.numeric(pcs +1))
    una_frase <- c("Recomended PCs to use")
    plot_df <- data.frame(pct = pct, 
                          cumu = cumu, 
                          rank = 1:length(pct))
    ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
      geom_text() + 
      scale_color_manual(values = c("dodgerblue4", "#ff2d00"))+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks = element_line(colour = "black",size = 0.5),
            axis.text = element_text(size = 12.5),
            axis.title = element_text(size = 12.5))+
      guides(color = guide_legend(override.aes = list(size = 5)))+
      annotate("text", x = 50, y = max(pct), 
               label = paste(una_frase, pcschar, sep = ":"),
               parse = F)
  })
  #----
  
  ##DoubletFinder----
  #Preprocessing
  seuobjfilt_norm_high_scale_dimred_prepross <- eventReactive(input$click6, {
    Neighbours <- FindNeighbors(object = seuobjfilt_norm_high_scale_dimred(),
                                dims = 1:input$pcstouse)
    Clusters <- FindClusters(object = Neighbours, 
                             resolution = c(0.3,0.5,0.7,1,1.5,2))
    seuobj <- RunUMAP(object = Clusters,
                      dims = 1:input$pcstouse) 
    return(seuobj)
  })
  output$soc.normal5 <- renderText({
    if(dim(seuobjfilt_norm_high_scale_dimred_prepross()) > 1) 
      "Data processing has been succesful"
  })
  
  #Pk optimization
  optimus_pk_to_plot <- eventReactive(input$click7, {
    sweep_list <- paramSweep_v3(seu = seuobjfilt_norm_high_scale_dimred_prepross(),
                                PCs = 1:input$pcstouse,
                                sct = TRUE)
    sweep_stats <- summarizeSweep(sweep.list = sweep_list,
                                  GT = FALSE)
    return(find.pK(sweep.stats = sweep_stats))
  })
  
  output$soc.normal6 <- renderText({
    if (dim(optimus_pk_to_plot()) > 1) 
      "The pK value that corresponds with the biggest BCmetric value should be used"
  })
  
  output$optimus_pk <- renderPlot({
    ggplot(optimus_pk_to_plot(), aes(pK,BCmetric, group =1))+
      geom_point(show.legend = F)+
      geom_line(show.legend = F)+
      #      geom_vline(xintercept = pK[BCmetric == max(BCmetric)], 
      #                color="red")+
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.5),
            axis.ticks = element_line(colour = "black",size = 0.5),
            axis.text = element_text(size = 12.5),
            axis.text.x = element_text(angle = 90),
            axis.title = element_text(size = 12.5))
    
  })
  
  #Expected doublets adjustment
  expected_adjusted <- eventReactive(input$click8, {
    annotations <- seuobjfilt_norm_high_scale_dimred_prepross()@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExpdoubl <- round(input$expected_doublets*nrow(seuobjfilt_norm_high_scale_dimred_prepross()@meta.data))
    return(round(nExpdoubl*(1-homotypic.prop)))
  })
  
  output$expectdoublets_adjusted <- renderText({
    una_frase <- c("The adjusted number of doublets expected is")
    paste(una_frase, expected_adjusted(), sep = ":")
  })
  
  #DoubletFinder
  seurat_doublet <- eventReactive(input$click9, {
    seuobj <- doubletFinder_v3(seuobjfilt_norm_high_scale_dimred_prepross(),
                               PCs = 1:input$pcstouse,
                               pN = 0.25,
                               pK = input$pk_optim,
                               nExp = input$expect_optim,
                               reuse.pANN = FALSE, sct = T)
    return(seuobj)
  })
  
  output$summary_the_doublet_action <- renderTable({
    table(seurat_doublet()@meta.data$"DF.classification")
  })
  
  seurat_doublet_removal <- eventReactive(input$click10, {
    DF.name = colnames(seurat_doublet()@meta.data)[grepl("DF.classification", colnames(seurat_doublet()@meta.data))]
    return(seurat_doublet()[,seurat_doublet()@meta.data[,DF.name] == "Singlet"])
  })
  
  output$dims_final <- renderText({
    if (dim(seurat_doublet_removal())[2] > 0 &&
        dim(seurat_doublet_removal())[2] < dim(seurat_doublet())[2]) 
      "Doublets have been removed succesfuly"
  })
  #----
  
  #Clusterization----
  #UMAP
  umaping <- eventReactive(input$click11, {
    if (input$radio == "Yes"){
      return(seurat_doublet_removal())}
    if (input$radio == "No"){
      Neighbours <- FindNeighbors(object = seuobjfilt_norm_high_scale_dimred(),
                                  dims = 1:input$pcstouse_noclist)
      Clusters <- FindClusters(object = Neighbours, 
                               resolution = c(0.3,0.5,0.7,1,1.5,2))
      seuobj <- RunUMAP(object = Clusters,
                        dims = 1:input$pcstouse_noclist)
      return(seuobj)}
  })
  
  
  output$plot_umap <- renderPlot({
    if (input$radio == "Yes"){
      assay <- DefaultAssay(umaping())
      res_used <- input$resolution_to_use
      group <- paste(assay,res_used, sep = "_")
      return(DimPlot(umaping(),
                     reduction = "umap",
                     group.by = group,
                     label = T))}
    if (input$radio == "No"){
      assay <- DefaultAssay(umaping())
      res_used <- input$resolution_to_use
      group <- paste(assay,res_used, sep = "_")
      return(DimPlot(umaping(),
                     reduction = "umap",
                     group.by = group,
                     label = T))}
  })
  
  #TSNE
  seurat_doublet_removat_tsne <- eventReactive(input$click12, {
    if (input$radio == "Yes"){
      seuobj <- RunTSNE(object = umaping(),
                        dims.use = 1:input$pcstouse,
                        do.fast = T)
      return(seuobj)}
    if (input$radio == "No"){
      seuobj <- RunTSNE(object = umaping(),
                        dims.use = 1:input$pcstouse_noclist,
                        do.fast = T)
      return(seuobj)}
  })
  
  tsneing <- eventReactive(input$click12, {
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    return(DimPlot(seurat_doublet_removat_tsne(),
                   reduction = "tsne",
                   group.by = group,
                   label = T))
  })
  output$plot_tsne <- renderPlot({
    tsneing()
  })
  
  #Save the data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Seuratobject.rds")
    },
    content = function(file) {
      saveRDS(object = seurat_doublet_removat_tsne(), file = file)
    })
  #----
  
  #Annotate cluster----
  marking <- eventReactive(input$click19, {
    
    seurat_doublet_removat_tsne <- seurat_doublet_removat_tsne()
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    
    seurat_doublet_removat_tsne <- SetIdent(seurat_doublet_removat_tsne(),
                                            value = group)
    
    markers <- FindAllMarkers(seurat_doublet_removat_tsne,
                              only.pos = T,
                              min.pct = input$pct_thres,
                              logfc.threshold = input$logf_thres)
    return(markers)
  })
  
  output$soc.normal7 <- renderText({
    if (dim(marking()) > 1) 
      "The markers computation has finished"
  })
  
  #Render little table
  marking_resumed <- eventReactive(input$click20, {
    markers <- marking()%>%
      group_by(cluster) %>%
      slice_max(n = input$max_genes , order_by = avg_log2FC)
    return(markers)
  })
  
  output$table_marked <- renderDataTable({
    marking_resumed()
  })
  
  #Render graph
  umap_feat <- eventReactive(input$click21, {
    markers <- marking()
    
    markers2 <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC)
    
    markers_to_plot <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 1, order_by = avg_log2FC)
    
    return(FeaturePlot(seurat_doublet_removat_tsne(),
                       features = markers_to_plot$gene))
  })
  
  output$umap_marked <- renderPlot({
    umap_feat()
  })
  #----
  
  #Single gene----
  #Ridge
  ridging <- eventReactive(input$click13, {
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    return(RidgePlot(seurat_doublet_removat_tsne(),
                     features = input$ridgevalue,
                     group.by = group))
  })
  output$ridge <- renderPlot({
    ridging()
  })
  
  #Violin
  violin <- eventReactive(input$click14, {
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    return(VlnPlot(seurat_doublet_removat_tsne(),
                   features = input$violinvalue,
                   group.by = group,
                   log = TRUE))
  })
  
  output$violin_plot <- renderPlot({
    violin()
  })
  
  #Feature
  feat <- eventReactive(input$click15, {
    return(FeaturePlot(seurat_doublet_removat_tsne(),
                       features = input$featinput))
  })
  output$featuring <- renderPlot({
    feat()
  })
  
  #Heatmap
  heatm <- eventReactive(input$click16,{
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    return(DoHeatmap(seurat_doublet_removat_tsne(),
                     features = input$heatinput,
                     group.by = group))
  })
  output$heatmap <- renderPlot({
    heatm()
  })
  #----
  
  #Coexpression----
  #Primer plot
  coexpress1 <- eventReactive(input$click17, {
    feature1 <- input$coexp1
    feature2 <- input$coexp2
    return(FeaturePlot(seurat_doublet_removat_tsne(),
                       features = c(feature1,feature2),
                       blend = T,
                       blend.threshold = 0.5))
  })
  output$coexpressionplot1 <- renderPlot({
    coexpress1()
  })
  
  #Segon plot
  coexpress2 <- eventReactive(input$click18, {
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    return(DotPlot(seurat_doublet_removat_tsne(),
                   features = input$coexp3,
                   group.by = group))
  })
  output$coexpressionplot2 <- renderPlot({
    coexpress2()
  })
  #----
  
  #Trajectory----
  #Dyno
  dinosaur_counts <- eventReactive(input$click23,{
    return(Matrix::t(as(as.matrix(seurat_doublet_removat_tsne()@assays$RNA@counts), 'sparseMatrix')))
  })
  dinosaur_expression <- eventReactive(input$click23, {
    return(Matrix::t(as(as.matrix(seurat_doublet_removat_tsne()@assays$RNA@data), 'sparseMatrix')))
  })
  dinosaur_cellinfo <- eventReactive(input$click23, {
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    seurat_doublet_removat_tsne <- SetIdent(seurat_doublet_removat_tsne(),
                                            value = group)
    
    return(object_cellinfo <- seurat_doublet_removat_tsne@active.ident)
  })
  
  dyno_object_create <- eventReactive(input$click23, {
    dyno_object <- wrap_expression(counts = dinosaur_counts(),
                                   expression = dinosaur_expression())
    return(dyno_object)
  })
  
  #  model_to_use <- eventReactive(input$click23,  {
  #    guidelines <- guidelines(dyno_object_create())
  #    guidelines <- guidelines[[1]][[1]]
  #    guidelines <- guidelines[1]
  #    return(guidelines)
  #    })
  
  #  output$select_dyno_method <- renderText({
  #    model_to_use()
  #  })
  
  dinosaur_model <- eventReactive(input$click23, {
    model <- infer_trajectory(dyno_object_create(), method = input$select_dyno_method)
    model <- model %>%
      add_root()
    model <- model %>%
      add_dimred(dimred = as.matrix(seurat_doublet_removat_tsne()@reductions$umap@cell.embeddings),
                 expression_source = dyno_object_create()$expression)
    return(model)
  })
  
  dinosaur <- eventReactive(input$click23, {
    return(plot_dimred(dinosaur_model(), 
                       expression_source = dyno_object_create()$expression, 
                       grouping = dinosaur_cellinfo(),
                       color_density = "grouping"))
  })
  
  output$dyno <- renderPlot({
    dinosaur()
  })
  
  #Dyno feat
  dinosaur_featuring <- eventReactive(input$click24, {
    return(plot_dimred(dinosaur_model(),
                       expression_source = dyno_object_create()$expression,
                       grouping = dinosaur_cellinfo(),
                       feature_oi = input$select_dyno_gene,
                       color_cells = "feature"))
  })
  
  output$dyno_feat <- renderPlot({
    dinosaur_featuring()
  })
  
  #Dyno heatmap
  dinosaur_heat <- eventReactive(input$click25, {
    return(plot_heatmap(
      dinosaur_model(),
      expression_source = dyno_object_create()$expression,
      grouping = dinosaur_cellinfo(),
      features_oi = input$features_dyno_heat
    ))
  })
  
  output$dyno_heat <- renderPlot({
    dinosaur_heat()
  })
  
  #----
  
  #CellChat----
  cellchatt <- eventReactive(input$click26, {
    assay <- DefaultAssay(seurat_doublet_removat_tsne())
    res_used <- input$resolution_to_use
    group <- paste(assay,res_used, sep = "_")
    
    seuobj_norm <- SetIdent(seurat_doublet_removat_tsne(), 
                            value = group)
    
    seuobj_norm <- RenameIdents(object = seuobj_norm,
                                `0` = "cluster0")
    cellchat <- createCellChat(object = seuobj_norm, group.by = "ident", 
                               assay = DefaultAssay(seuobj_norm))
    
    
    if (input$mitocondrial_pattern == "^mt-"){
      CellChatDB <- CellChatDB.mouse}
    if (input$mitocondrial_pattern == "^MT-"){
      CellChatDB <- CellChatDB.human}
    
    #CellChatDB <- CellChatDB.mouse
    CellChatDB.use <- CellChatDB
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.mouse)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    
    cellchat <- filterCommunication(cellchat, min.cells = input$min_cell_comunication)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    return(cellchat)
  })
  
  netw1 <- eventReactive(input$click26, {
    groupSize <- as.numeric(table(cellchatt()@idents))
    return(netVisual_circle(cellchatt()@net$count, 
                            vertex.weight = groupSize, 
                            weight.scale = T, label.edge= F))
  })
  
  output$twoplotscells_1 <- renderPlot({
    netw1()
  })
  
  netw2 <- eventReactive(input$click26, {
    groupSize <- as.numeric(table(cellchatt()@idents))
    return(netVisual_circle(cellchatt()@net$weight, 
                            vertex.weight = groupSize, 
                            weight.scale = T, label.edge= F))
  })
  
  output$twoplotscells_2 <- renderPlot({
    netw2()
  })
  
  output$networks <- renderText({
    cellchatt()@netP$pathways
  })
  
  netwcircular <- eventReactive(input$click27, {
    return(netVisual_aggregate(cellchatt(), signaling = input$networkdisplay, layout = "circle"))
  })
  netwchord <- eventReactive(input$click27, {
    return(netVisual_aggregate(cellchatt(), signaling = input$networkdisplay, layout = "chord"))
  })
  
  output$circular <- renderPlot({
    netwcircular()
  })
  output$chord <- renderPlot({
    netwchord()
  })
  #----
  
  #----
}


#----

shinyApp(ui = ui, server = server)