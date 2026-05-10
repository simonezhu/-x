# 基础包
library(ggplot2)
library(ggrepel)
library(affy)
library(hgu133a2.db)
library(org.Hs.eg.db)
# 数据质量检测包
library(arrayQualityMetrics)
# 数据处理包
library(limma)
# 火山图包
library(EnhancedVolcano)
# 热图包
library(pheatmap)
# 代谢分析包
library (clusterProfiler)
library(enrichplot)

# 数据读取与预处理
mydata <- ReadAffy(celfile.path="/home/simone/microarrays/GSE24460") #把下载下来的.cel.gz文件转化为affybatch文件
eset <- rma(mydata) #标准化
exp_matrix <- exprs(eset) # 原始探针水平矩阵

# 数据质量检验
arrayQualityMetrics (mydata,outdir="quality_assement") #在/home/simone生成quality_assement

# 获取所有探针对应的 Symbol
all_symbols <- mapIds(hgu133a2.db,  # 芯片型号
                      keys = rownames(exp_matrix), #取exp_matrix的列名称为对象
                      keytype = "PROBEID", # 把探针ID
                      column = "SYMBOL")# 转换为基因名称（或者可以用GENENAME）
                      #multiVals = "first") #这行是糊弄事的，如果想要糊弄甲方就可以用，然后就不用针对非特异性结合再数据处理了

# 数据处理，针对潜在的探针非特异性结合造成的误差
exp_matrix_gene <- avereps(exp_matrix, ID = all_symbols) # 使用 avereps 对表达矩阵按 Symbol 求平均
exp_matrix_gene <- exp_matrix_gene[!is.na(rownames(exp_matrix_gene)), ] # 剔除掉没有 Symbol 的探针，并合并重复基因

# 差异表达分析 (基于去重后的基因矩阵)
group <- factor(c(rep("Control", 2), rep("Treat", 2))) # 2个对照组，2个实验组
design <- model.matrix(~group) #生成“disign”这个矩阵

# 使用limma对去重后的矩阵进行线性拟合
fit <- lmFit(exp_matrix_gene, design) # 输入“design”，将其赋值给“fit”
fit <- eBayes(fit) # 用经验贝叶斯模型对“fit”进行处理
res <- topTable(fit, coef=2, number=Inf) #生成数据框“res”，设定提取“fit”第2列对应的比较结果（组间差异），提取所有基因的计算结果（inf）
res$Symbol <- rownames(res) # 在“res”中生成一个叫“res$Symbol”的新列，其显示基因名简称

# 给volcano map和pheat map筛选用于展示的 Top 基因 (从去重后的结果中选)
res <- as.data.frame(res) # 将生成的对象强制转换为标准的数据框
up_top10 <- head(res[order(res$P.Value), ], 100) # 按 P 值从小到大（最显著到最不显著）排序，取排序后的前 100 个基因作为候选
up_top10 <- head(up_top10[up_top10$logFC > 0, ], 10)$Symbol # 在候选池中筛选出 logFC 大于 0 的基因，从这些上调基因中取前 10 个，只提取它们的基因名
down_top10 <- head(res[order(res$P.Value), ], 100) # 同上
down_top10 <- head(down_top10[down_top10$logFC < 0, ], 10)$Symbol #在候选池筛选出 logFC 小于 0 的基因，从这些上调基因中取前 10 个，只提取它们的基因名
top_genes <- c(up_top10, down_top10) # 将两个字符向量（10 个上调 + 10 个下调）拼接在一起，生成向量“top_genes”

# volcanno map
EnhancedVolcano(res, #取“res”为对象
                lab = res$Symbol, # 指定图中每一个点对应的标签名
                x = 'logFC', # 设置x轴为“logFC”
                y = 'P.Value', # 设置y轴为“P.Value”
                title = '基因水平差异分析 (去重后)', # 标题
                selectLab = top_genes, # 显示“top_genes”
                drawConnectors = TRUE, # 当标签离点较远时，画一条直线连着，防止混淆
                pCutoff = 0.05, # p=0.05
                FCcutoff = 1.0) # 显示是否有差异的分界线为log2FC=1

# pheatmap
plot_matrix <- exp_matrix_gene[top_genes, ] # 创建矩阵“plot_matrix”，其为“exp_matrix_gene”中“top_genes”所对应的那些
annotation_col <- data.frame(Group = group) # 创建一个注释条
rownames(annotation_col) <- colnames(plot_matrix) # 注释条会在热图的最上方显示一排色块，标明哪些列是对照组，哪些是实验组
bk <- c(seq(-2, -0.1, length.out=50), seq(0, 2, length.out=50))
pheatmap(plot_matrix, # 取“plot_matrix”为对象
         scale = "row", # 归一化              
         clustering_distance_rows = "correlation", # 根据基因表达的相关性进行聚类，把表达模式相似的基因排在一起
         annotation_col = annotation_col, 
         breaks = bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 设置热图颜色
         main = "Top DEGs Heatmap (Gene Level)")

# 使用bitr将基因名转换为ID以便clusterProfiler后续处理（e.g.GO）
geneID <- bitr(res$Symbol, # 创建变量“geneID”，以“res$Symbol”为其赋值
      fromType="SYMBOL", # 将基因名简称
      toType="ENTREZID", # 转换为ID
       OrgDb = org.Hs.eg.db) # 样本来自智人

# GO
ego <- enrichGO(gene          = geneID$ENTREZID, # 输入“geneID”
                OrgDb         = org.Hs.eg.db, # 样本来自智人
                ont           = "BP", # 在Biological Process（BP）层面分析
                pAdjustMethod = "BH", # 使用Benjamini-Hochberg（BH）控制假阳性结果
                pvalueCutoff  = 0.05) # 设定p=0.05
dotplot(ego, showCategory = 20) + ggtitle("GO Pathway Enrichment")# 生成名为“GO Pathway Enrichment”的GO的气泡图
write.csv(as.data.frame(ego), "my_go_results.csv") # 生成一份名为ego的GO的csv文件

# KEEG
kk <- enrichKEGG(gene         = geneID$ENTREZID, # 输入“geneID”
                 organism     = 'hsa',   # 样本来自智人
                 pvalueCutoff = 0.05, # 设定p=0.05
                 qvalueCutoff = 0.2) # 多重假设检验矫正后的阈值=0.02
dotplot(kk, showCategory = 20) + ggtitle("KEGG Pathway Enrichment") # 生成名为"KEGG Pathway Enrichment"的KEEG气泡图
write.csv(as.data.frame(kk), "my_kegg_results.csv") #生成一份名为kk的KEGG的csv文件
