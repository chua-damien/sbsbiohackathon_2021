#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


# Define UI for application that draws a histogram
shinyUI(
    navbarPage(title = "Data Exploration",
               
        tabPanel(title = "Gene Expression",
            
            sidebarPanel(
                selectizeInput(
                    inputId = "geneID",
                    label = "Input geneID",
                    choices = NULL,
                    selected=NULL,
                    multiple = F, # allow for multiple inputs
                    options = list(create = FALSE, placeholder = "E.g AH39280") # if TRUE, allows newly created inputs
                ),
                selectizeInput(
                    inputId = "treatment",
                    label = "Select treatment type",
                    choices = NULL,
                    selected=NULL,
                    multiple = T, # allow for multiple inputs
                    options = list(create = FALSE) # if TRUE, allows newly created inputs
                ),
                actionButton("gene.exp", "Plot")
            ),
            
            mainPanel(
                plotOutput("GeneExpPlot")
            )
        )
    )
)