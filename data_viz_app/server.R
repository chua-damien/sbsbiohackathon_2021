#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
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
library(d3heatmap)
library(VennDiagram)

#read data
load("data/all_data.RData")

# Define server logic
shinyServer(function(input, output, session) {
    
    updateSelectizeInput(session, "pairwise", choices= unique(dfCombine$pairwise),
                         selected=character(0), server = T)
    updateSelectizeInput(session, "geneExpID", choices = rownames(cjejuniEM), selected=character(0), server = TRUE)
    
    output$DFdeg <- renderDataTable({
        if (input$Pairwise2 == 0 )
            return()
        
        input$Pairwise2
        #Get inputs
        pairWiseGrp <- isolate(input$pairwise)
        log2fc <- as.numeric(isolate(input$log2fc_input))
        pval <- as.numeric(isolate(input$deg_pval))
        
        #Plotdf
        if (pairWiseGrp == "SsSe"){
            dfFilter <- left_join(resmfSsSe, combineTerms, by = "geneID") %>%
                dplyr::select(-c("baseMean", "stat", "pvalue"))
            dfFilter <- dfFilter[, c(4, 1, 2, 3, 5, 6)]
            colnames(dfFilter) <- c("GeneID", "log2FC", "logfc.SE", "P-adjusted", "Gene Description", "GOterms")
            
            dfFilter <- dfFilter[which(abs(dfFilter$log2FC) > log2fc) ,]
            dfFilter <- dfFilter[which(dfFilter[,4] < pval) ,]
            
            dfFilter
        }
        else {
            tmp <- as.data.frame(resmfSsFs)
            tmp$geneID <- rownames(tmp)
            dfFilter <- left_join(tmp, combineTerms, by = "geneID") %>%
                dplyr::select(-c("baseMean", "stat", "pvalue"))
            dfFilter <- dfFilter[, c(4, 1, 2, 3, 5, 6)]
            colnames(dfFilter) <- c("GeneID", "log2FC", "logfc.SE", "P-adjusted", "Gene Description", "GOterms")
            
            dfFilter <- dfFilter[which(abs(dfFilter$log2FC) > log2fc) ,]
            dfFilter <- dfFilter[which(dfFilter[,4] < pval) ,]
            
            dfFilter
        }
        })
    
    output$Volcanoplot <- renderPlot({
        if (input$Pairwise2 == 0 ) 
            return()
        
        input$Pairwise2
        #Get inputs
        pairWiseGrp <- as.character(isolate(input$pairwise))
        
        #Plotvolcano
        if (pairWiseGrp == "SsSe") {
            
            resmfSsSe <- as.data.frame(resmfSsSe) #convert to df
            resmfSsSe$geneID <- rownames(resmfSsSe) #convert col1 to genes
            
            log2fc <- as.numeric(isolate(input$log2fc_input))
            pval <- as.numeric(isolate(input$deg_pval))
            
            #Custom colors
            keyvals.colour2 <- ifelse(
                resmfSsSe$log2FoldChange < -log2fc, 'royalblue',
                ifelse(resmfSsSe$log2FoldChange > log2fc, 'gold',
                       'black')) 
            names(keyvals.colour2)[keyvals.colour2 == 'gold'] <- 'Upregulated'
            names(keyvals.colour2)[keyvals.colour2 == 'black'] <- 'No biological sig'
            names(keyvals.colour2)[keyvals.colour2 == 'royalblue'] <- 'Downregulated'
            
            # #Resize boundaries
            # dev.new(width=20, height=20)
            
            #Plotting
            EnhancedVolcano(resmfSsSe,
                            lab = as.character(resmfSsSe$geneID),
                            x = 'log2FoldChange',
                            y = 'padj',
                            selectLab = as.character(top_n(resmfSsSe, -10, resmfSsSe$padj)[,7]),
                            xlim = c(-8,8),
                            xlab = expression(Log[2]~Fold~Change),
                            ylim = c(0, 100),
                            title = "DESeq2",
                            pCutoff = pval,
                            #pLabellingCutoff = 'p-value cutoff',
                            FCcutoff = log2fc,
                            pointSize = 2.0,
                            labSize = 3.0,
                            labCol = 'black',
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colCustom = keyvals.colour2,
                            legend=c('NS',expression(Log[2]~FC),'P-value',
                                     expression(p-value~and~log[2]~FC)),
                            legendLabels = TRUE,
                            legendPosition = 'right',
                            legendLabSize = 14,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 1.0,
                            colConnectors = 'black')
        } else {
            resmfSsFs <- as.data.frame(resmfSsFs) #convert to df
            resmfSsFs$geneID <- rownames(resmfSsFs) #convert col1 to genes
            log2fc <- as.numeric(isolate(input$log2fc_input))
            pval <- as.numeric(isolate(input$deg_pval))
            #Custom colors
            keyvals.colour2 <- ifelse(
                resmfSsFs$log2FoldChange < -1.5, 'royalblue',
                ifelse(resmfSsFs$log2FoldChange > 1.5, 'gold',
                       'black')) 
            names(keyvals.colour2)[keyvals.colour2 == 'gold'] <- 'Upregulated'
            names(keyvals.colour2)[keyvals.colour2 == 'black'] <- 'No biological sig'
            names(keyvals.colour2)[keyvals.colour2 == 'royalblue'] <- 'Downregulated'
            
            # #Resize boundaries
            # dev.new(width=20, height=20)
            # 
            #Plotting
            EnhancedVolcano(resmfSsFs,
                            lab = as.character(resmfSsFs$geneID),
                            x = 'log2FoldChange',
                            y = 'padj',
                            selectLab = as.character(top_n(resmfSsFs, -10, resmfSsFs$padj)[,7]),
                            xlim = c(-8,8),
                            xlab = expression(Log[2]~Fold~Change),
                            ylim = c(0, 100),
                            title = "DESeq2",
                            pCutoff = pval,
                            #pLabellingCutoff = 'p-value cutoff',
                            FCcutoff = log2fc,
                            pointSize = 2.0,
                            labSize = 3.0,
                            labCol = 'black',
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colCustom = keyvals.colour2,
                            legend=c('NS',expression(Log[2]~FC),'P-value',
                                     expression(p-value~and~log[2]~FC)),
                            legendLabels = TRUE,
                            legendPosition = 'right',
                            legendLabSize = 14,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 1.0,
                            colConnectors = 'black')
        }
    })  
    
    #plotting heatmap
    
    output$DEGHeatmap <- renderD3heatmap({
        if (input$Pairwise2 == 0 ) 
            return()
        
        input$Pairwise2
        #Get inputs
        pairWiseGrp <- as.character(isolate(input$pairwise))
        log2fc <- as.numeric(isolate(input$log2fc_input))
        pval <- as.numeric(isolate(input$deg_pval))
        
        # Filtering with user input
        resmfSsSe <- as.data.frame(resmfSsSe) #convert to df
        resmfSsSe$geneID <- rownames(resmfSsSe) #convert col1 to genes
        
        dfFilter <- dfCombine %>% filter(pairwise == pairWiseGrp) %>%
            dplyr::select(-c("baseMean", "stat", "pvalue", "pairwise"))
        colnames(dfFilter) <- c("GeneID", "log2FC", "logfc.SE", "P-adjusted", "Gene Description", "GOterms")
        
        dfFilter.pos <- dfFilter[which(dfFilter$log2FC > 0 ) ,]
        dfFilter.neg <- dfFilter[which(dfFilter$log2FC < 0 ) ,]
        dfFilter <- dfFilter[c(which(dfFilter.pos$log2FC > log2fc ), which(dfFilter.neg$log2FC < -log2fc )) ,]
        
        
        dfFilter_tpm <- exp.mat.campy.tpm[which(rownames(exp.mat.campy.tpm) == dfFilter$GeneID),]
        
        d3heatmap(exp.mat.campy.tpm)
        
        # #Resize boundaries
        # dev.new(width=100, height=200)
        # 
        # #PlottingHeatmap
        # 
        # renderD3heatmap({d3heatmap(exp.mat.campy.tpm)})
        
        #Heatmap(as.matrix(dfFilter_tpm))
        
        #Heatmap(as.matrix(upreg_all_GO_v4),col = pvalue_col_fun, name = "-log(p-value)", 
        #        width = unit(4, "cm"), 
        #        row_names_side = "left", row_names_gp = gpar(fontsize = 9), 
        #        row_dend_side = "right", 
        #        column_names_centered = TRUE, column_names_rot = 45, column_names_gp = gpar(fontsize = 10.5),
            
        #        left_annotation = ha )
        
        
    })
    
    #Venn tab
    output$drawVenn <- renderPlot({
        draw.pairwise.venn(nrow(dfSsSe),
                           nrow(dfSsFs),
                           sum(dfSsFs$geneID %in% dfSsSe$geneID),
                           category = c("SsFs", "SsSe"), scaled= F, cat.pos = c(1,1),
                           fill = c("#3b5a9d", "#ffa633"))
    })
    
    output$vennTable <- renderDataTable({
        if (input$vennButton == 0 )
            return()
        #Get inputs
        grpChoice <- isolate(input$vennRadio)
        base::print(grpChoice)
        #Plotdf
        venndf <- pairWisedata %>% filter(pairwise == grpChoice) %>%
            dplyr::select(-c("pairwise", "descr", "GO"))
        
        notableGenes <- goTermsall[(goTermsall$geneID %in% venndf$geneID), 1:3]
        venndf <- left_join(venndf, notableGenes, by="geneID")
        
        venndf})
    
    #Gene tab
    output$GeneExpPlot <- renderPlot({
        if (input$geneCountButton == 0)
            return()
        
        input$geneCountButton
        geneID <- isolate(input$geneExpID)
        compare <- isolate(input$geneComparison)
        if (compare == "SsSe"){
            d <- DESeq2::plotCounts(ddsmfSsSe, gene=geneID, intgroup="Treatment", returnData = T)
            ggplot(d, aes(x=Treatment, y=count)) + 
                geom_violin(mapping=aes(fill=Treatment), alpha = 0.6) + 
                geom_jitter(height = 0, width = 0.1) + 
                theme_light() + theme(legend.position = "none") + 
                ylab("Normalized Counts")
        }
        else {
            d <- DESeq2::plotCounts(ddsmfSsFs, gene=geneID, intgroup="Treatment", returnData = T)
            ggplot(d, aes(x=Treatment, y=count)) + 
                geom_violin(mapping=aes(fill=Treatment), alpha = 0.6) + 
                geom_jitter(height = 0, width = 0.1) + 
                theme_light() + theme(legend.position = "none") +
                ylab("Normalized Counts")
        }
    })
    
    #PCA tab
    output$PCAPlot <- renderPlot({
        if (input$PCAButton == 0)
            return()
        
        input$PCAButton
        comparePCA <- isolate(input$PCAComparison)
        
        if (comparePCA == "SsSe"){
            p <- DESeq2::plotPCA(vsdmfSsSe,intgroup=c("Strain", "Treatment"))
            p + theme_light()
        }
        else {
            p <- DESeq2::plotPCA(vsdmfSsFs,intgroup=c("Strain", "Treatment"))
            p + theme_light()
        }
        
    })
})
