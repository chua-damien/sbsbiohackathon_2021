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
library(InteractiveComplexHeatmap)
library(VennDiagram)

#read data
#file.dir <- file.path("~/NTU/_Biohackathon2021/data")
#load("~/NTU/_Biohackathon2021/app/all_data.RData")

#exp.mat.campy.tpm <- read.csv(paste0(file.dir,"/expmat(c.jejuni_tpm).csv"))
#rownames(exp.mat.campy.tpm) <- exp.mat.campy.tpm[,1]
#exp.mat.campy.tpm<-  exp.mat.campy.tpm[,2:ncol(exp.mat.campy.tpm)]

#gene.anno.campy <- read.csv(paste0(file.dir,"/c_jejuni_genedescr.txt"), sep = "\t", header = F,
#                            col.names = c("geneID", "Description"))
#gene.anno.campy<- gene.anno.campy[match(rownames(exp.mat.campy.tpm),gene.anno.campy$geneID),]
#gene.anno.campy[which(is.na(gene.anno.campy$Description)),2] <- "NA"
#gene.anno.campy$geneID <- rownames(exp.mat.campy.tpm)
#metadata.campy <- read.csv(paste0(file.dir,"/exp_annot_c.jejuni_new.txt"), sep = "\t")

# Define server logic

server<- shinyServer(function(input, output, session) {
    
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
                resmfSsFs$log2FoldChange < -log2fc, 'royalblue',
                ifelse(resmfSsFs$log2FoldChange > log2fc, 'gold',
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
    
    #plotting heatmap
    observeEvent(input$Pairwise2, {
        
        #Get inputs
        pairWiseGrp <- as.character(isolate(input$pairwise))
        log2fc <- as.numeric(isolate(input$log2fc_input))
        pval <- as.numeric(isolate(input$deg_pval))
        
        if (pairWiseGrp == "SsSe"){
            tmp <- as.data.frame(resmfSsSe)
        }
        else {
            tmp <- as.data.frame(resmfSsFs)
        }
        # Filtering with user input
        resmfSsSe <- as.data.frame(resmfSsSe) #convert to df
        resmfSsSe$geneID <- rownames(resmfSsSe) #convert col1 to genes
        
        
        tmp$geneID <- rownames(tmp)
        dfFilter <- left_join(tmp, combineTerms, by = "geneID") %>%
            dplyr::select(-c("baseMean", "stat", "pvalue"))
        dfFilter <- dfFilter[, c(4, 1, 2, 3, 5, 6)]
        colnames(dfFilter) <- c("GeneID", "log2FC", "logfc.SE", "P-adjusted", "Gene Description", "GOterms")
        
        dfFilter.pos <- dfFilter[which(dfFilter$log2FC > 0 ) ,]
        dfFilter.neg <- dfFilter[which(dfFilter$log2FC < 0 ) ,]
        dfFilter <- dfFilter[c(which(dfFilter.pos$log2FC > log2fc ), which(dfFilter.neg$log2FC < -log2fc )) ,]
        
        dfFilter_tpm <- exp.mat.campy.tpm[which(rownames(exp.mat.campy.tpm) %in% dfFilter$GeneID),]
        
        
        ha <- HeatmapAnnotation(df = data.frame(Strain = metadata.campy$Strain,
                                                Morphotype = metadata.campy$Treatment) ,
                                col = list(Strain = c("NCTC12661" = "orange",
                                                      "PT14" = "cyan"),
                                           Morphotype = c("Spiral morphotype exponential phase" = "dark green", 
                                                          "Spiral morphotype stationary phase" = "blue",
                                                          "Filamentous morphotype stationary phase" = "red")),
                                gp = gpar(col = "black"),
                                show_annotation_name = TRUE)
        
        gene_ha <- rowAnnotation(df = data.frame(gene_description = gene.anno.campy[which(gene.anno.campy$geneID %in% rownames(dfFilter_tpm)),]$Description),
                                 show_annotation_name = FALSE,
                                 simple_anno_size = unit(0.1, "cm"),
                                 annotation_width = unit(0.00000001, "cm"),
                                 show_legend = FALSE)
        
        #col_fun = colorRamp2(c(0, 2000, 8000), c("blue", "white", "red"))
        #lgd = Legend(col_fun = col_fun, title = "foo")
        
        col_fun = circlize::colorRamp2(c(0,300,8000), c("blue", "white", "red"))
        ht1 = Heatmap(as.matrix(dfFilter_tpm), name = "TPM_values", col = col_fun, 
                      top_annotation = ha,
                      left_annotation = gene_ha,
                      show_row_names = FALSE, show_column_names = FALSE
                      
                      )
        
        InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, ht1, "DEGHeatmap"
        )
        
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




#shinyApp(ui = ui, server = server)
