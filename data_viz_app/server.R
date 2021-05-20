#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggfortify)

#read data
metadata <- read.csv("../data/Campylobacter jejuni_exp_anno.txt", sep = "\t")
exp.mat <- read.csv("../data/Campylobacter jejuni.txt", sep = "\t", row.names = 1)
gene.anno <- read.csv("../data/Campylobacter jejuni_gene_anno.txt", sep = "\t", header = F,
                      col.names = c("geneID", "Description"))


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    updateSelectizeInput(session, "geneID", choices = rownames(exp.mat), selected=character(0), server = TRUE)
    
    updateSelectizeInput(session, "treatment", choices = unique(metadata$Treatment)[unique(metadata$Treatment) != ""], 
                         selected=character(0), server = TRUE)
    
    updateSelectizeInput(session, "pca.treatment", choices = unique(metadata$Treatment)[unique(metadata$Treatment) != ""], 
                         selected=character(0), server = TRUE)
    
    output$GeneExpPlot <- renderPlot({
        if (input$plot.gene.exp == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$plot.gene.exp
        #getting inputs
        geneID <- isolate(input$geneID)
        expt.treatment <- isolate(input$treatment)
        #making df for plotting
        meta.plot <- metadata[metadata$Treatment %in% expt.treatment,]
        exp.mat.plot <- t(exp.mat[geneID, colnames(exp.mat) %in% meta.plot$Run])
        exp.mat.plot <- cbind(meta.plot, exp.mat.plot)
        colnames(exp.mat.plot)[ncol(exp.mat.plot)] <- "TPM"
        #plot 
        ggplot(data = exp.mat.plot) + 
            facet_wrap(~ Treatment, scales = "free_x") +
            geom_point(mapping = aes(x=Run, y = TPM, col = Treatment)) + 
            theme_light() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  legend.text=element_text(size=8),
                  strip.background = element_rect(colour="grey80", fill="grey80"),
                  strip.text = element_text(colour = 'black'))
    })
    
    output$PCAPlot <- renderPlot({
        if (input$plot.pca == 0)
            return()
        
        #freeze the input till confirmation is given to plot
        input$plot.pca
        #getting inputs
        pca.expt.treatment <- isolate(input$pca.treatment)
        pca.meta.plot <- metadata[metadata$Treatment %in% pca.expt.treatment,]
        pca.exp.mat.plot <- exp.mat[, colnames(exp.mat) %in% pca.meta.plot$Run]
        pca.exp.mat.plot <- pca.exp.mat.plot[rowSums(pca.exp.mat.plot) > 0,]
        pca <- prcomp(t(pca.exp.mat.plot), scale. = T, center = T)
        ggplot() + 
            geom_point(mapping = aes(pca$x[,1], pca$x[,2], color = pca.meta.plot$Treatment)) + 
            labs(col = "Treatment") + 
            xlab(paste0("PC1, ", round(pca$sdev[1] / sum(test$sdev), digits = 2), "% variation")) +
            ylab(paste0("PC2, ", round(pca$sdev[2] / sum(test$sdev), digits = 2), "% variation")) 
    })
})
