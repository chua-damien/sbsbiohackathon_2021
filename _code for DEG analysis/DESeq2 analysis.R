#Package required 
library(DESeq2)

#Load files 
cjejuniEM <- read.csv('expmat_cjejuni.csv', header=T)
cjejuniExp <- read.delim('exp_annot_cjejuni_new.txt', header=T, sep="\t")

#Experimental annotation 
coldata <- cjejuniExp[,c(1,2,7)]
coldata$Treatment <- factor(coldata$Treatment, 
                            levels = c("Spiral morphotype stationary phase",
                                       "Spiral morphotype exponential phase",
                                       "Filamentous morphotype stationary phase"))
rownames(coldata) <- coldata$Run
cjejuniEM <- as.matrix(cjejuniEM[,2:51])
cjejuniEM[,] <- as.integer(cjejuniEM[,])

#Filtering coldata / cjejuniEM for Spiral / Filamentatous stationary 
filtcoldataSsFs <- coldata[coldata$Treatment %in% c("Spiral morphotype stationary phase",
                                                    "Filamentous morphotype stationary phase"),]
rownames(filtcoldataSsFs) <- filtcoldataSsFs[,1]
filtEMdataSsFs <- cjejuniEM[,colnames(cjejuniEM) %in% rownames(filtcoldataSsFs)]

#Deseq2 object - SsFs
ddsSsFs <- DESeqDataSetFromMatrix(countData = filtEMdataSsFs,
                                  colData = filtcoldataSsFs,
                                  design= ~ Strain + Treatment + Strain:Treatment)
ddsSsFs <- ddsSsFs[rowSums(counts(ddsSsFs)) > 1,]
#Re-analyzing with multi-factorial design 
levels(ddsmfSsFs$Strain)
design(ddsmfSsFs) <- formula(~ Strain + Treatment)
ddsmfSsFs <- DESeq(ddsmfSsFs)
resmfSsFs <- results(ddsmfSsFs)



#Filtering coldata / cjejuniEM for Spiral stationary / Spiral exponential 
filtcoldataSsSe <- coldata[coldata$Treatment %in% c("Spiral morphotype stationary phase",
                                                    "Spiral morphotype exponential phase"),]
rownames(filtcoldataSsSe) <- filtcoldataSsSe[,1]
filtEMdataSsSe <- cjejuniEM[,colnames(cjejuniEM) %in% rownames(filtcoldataSsSe)]
#Sanity check
colnames(filtEMdataSsSe) == rownames(filtcoldataSsSe)

#Deseq2 object - SsSe
ddsSsSe <- DESeqDataSetFromMatrix(countData = filtEMdataSsSe,
                                  colData = filtcoldataSsSe,
                                  design= ~ Strain + Treatment + Strain:Treatment)
ddsSsSe <- ddsSsSe[rowSums(counts(ddsSsSe)) > 1,]

#Multi-factorial analysis
levels(ddsmfSsSe$Strain)
design(ddsmfSsSe) <- formula(~ Strain + Treatment)
ddsmfSsSe <- DESeq(ddsmfSsSe)
resmfSsSe <- results(ddsmfSsSe)
