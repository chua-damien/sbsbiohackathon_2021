#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(RColorBrewer)
library(EnhancedVolcano)
library(ggplot2)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(VennDiagram)
library(shinythemes)

# Define UI for application
ui <-shinyUI(
    navbarPage(title = "Data Exploration",
        
        tabPanel(title = "Home",
        fluidPage(theme = shinytheme("cosmo"),
                  
            titlePanel(
                             
            tags$h1(HTML(paste0("Megacephalus - ", tags$i("Campylobacter jejuni"))), align = "center")),
        includeHTML("C:/Users/Damien/Documents/NTU/_Biohackathon2021/app/homepagev3.html"))
        ),
        
        tabPanel(title = "Gene Counts",
                 sidebarLayout(
                     sidebarPanel(
                         selectizeInput(
                             inputId = "geneExpID",
                             label = "Input geneID",
                             choices = NULL,
                             selected=NULL,
                             multiple = F, # allow for multiple inputs
                             options = list(create = FALSE, placeholder = "E.g OIT94223") # if TRUE, allows newly created inputs
                         ),
                         radioButtons(
                             inputId = "geneComparison",
                             label = "Choose Comparison:",
                             choiceNames = c("Ss-Se", "Ss-Fs"),
                             choiceValues = c("SsSe", "SsFs"),
                             selected = "SsSe"),
                         actionButton("geneCountButton", "Get Data")
                     ),
                     mainPanel(plotOutput("GeneExpPlot"))
                 )
        ),
        
        
        tabPanel(title ="DGE analysis",
            sidebarLayout(
                sidebarPanel(
                    selectizeInput(
                        inputId="pairwise",
                        label = "Input desired pairwise comparison.",
                        choices = c("Spiral Stationary VS Spiral Exponential",
                                    "Spiral Stationary VS Filamentous Exponential"),
                        selected=NULL,
                        multiple = F,
                        options = list(create=F, placeholder= "EG: SsSe")
                    ),
                    
                    sliderInput(
                        inputId = "log2fc_input",
                        label = "Log2 fold change threshold",
                        value = 1.5,
                        min = 0,
                        max = 5,
                        step = 0.01
                    ),
                    
                    sliderInput(
                        inputId = "deg_pval",
                        label = "p-value threshold",
                        value = 0.05,
                        min = 0,
                        max = 0.1,
                        step = 0.001
                    ),
                    
                    actionButton("Pairwise2", "Get Data")
                ),
                
                mainPanel(
                    tabsetPanel(
                        tabPanel("Filtered Dataframe", dataTableOutput("DFdeg")),
                        tabPanel("Volcano plot", plotOutput("Volcanoplot")), 
                        tabPanel("Heatmap", InteractiveComplexHeatmapOutput("DEGHeatmap", 
                                                                            output_ui_float = TRUE, action = "hover")
                        )
                        
                        
                        )
                        
                        
                        
                        
                        
                        
                    ) 
                
                
            
        )),
        
        tabPanel(title = "Venn Diagram",
            verticalLayout(
                fluidRow(
                    column(12, align = "center",plotOutput("drawVenn",width = "50%",))
                ),
                sidebarLayout(
                     sidebarPanel(
                        radioButtons("vennRadio", "Choose condition(s):",
                            choiceNames = c("SsSe only", "SsFs only" , "SsSs&SsFe"),
                            choiceValues = c("SsSe", "SsFs", "pair"),
                            selected = "SsSe"),
                        actionButton("vennButton", "Get Data!")
                    ), 
                    mainPanel(dataTableOutput("vennTable"))
                )
            )
        ),
        
        tabPanel(title = "PCA Plot",
            sidebarLayout(
                sidebarPanel(
                    radioButtons(
                        inputId = "PCAComparison",
                        label = "Choose Comparison:",
                        choiceNames = c("Ss-Se", "Ss-Fs"),
                        choiceValues = c("SsSe", "SsFs"),
                        selected = "SsSe"),
                    actionButton("PCAButton", "Plot")
                ),
                mainPanel(plotOutput("PCAPlot"))
            )
        )
    ))



