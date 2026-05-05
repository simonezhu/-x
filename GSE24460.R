
library(ggplot2)
library(ggrepel)
library(limma)
library(affy)
library(arrayQualityMetrics)
library (EnhancedVolcano)
library(hgu133a2.db)

mydata<-ReadAffy(celfile.path="/home/simone/microarrays" 
)
mydata
eset <- rma(mydata)
exp_matrix <- exprs(eset)


group <- factor(c(rep("Control", 2), rep("Treat", 2)))
design <- model.matrix(~group)
fit <- lmFit(eset, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, number=Inf)
res

res <- as.data.frame(res)
res$logFC <- as.numeric(as.character(res$logFC))
res$P.Value <- as.numeric(as.character(res$P.Value))
res$adj.P.Val <- as.numeric(as.character(res$adj.P.Val))

columns(hgu133a2.db)
res$Symbol <- mapIds(hgu133a2.db,
                     keys = rownames(res),      
                     column = "TP53",        
                     keytype = "PROBEID",       
                     multiVals = "first")   

up_top10 <- head(res[order(res$P.Value), ], 100) 
up_top10 <- head(up_top10[up_top10$logFC > 0, ], 10)$Symbol
down_top10 <- head(res[order(res$P.Value), ], 100)
down_top10 <- head(down_top10[down_top10$logFC < 0, ], 10)$Symbol
top_genes <- c(up_top10, down_top10)

EnhancedVolcano(res,
                lab = res$Symbol,
                x = 'logFC', 
                y = 'P.Value',
                title = '差异表达分析',
                selectLab = top_genes,
                drawConnectors = TRUE,
                pCutoff = 0.05,
                FCcutoff = 1.0)
