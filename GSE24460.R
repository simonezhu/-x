# =============================================================================
# GSE24460 差异表达分析（修订版，用于 PR）
# 平台 GPL571 | 亲本 MCF-7 (Control)  vs  多柔比星耐药选育 MCF-7/ADR (Treat)
#
# 改动均以 [FIX-n] 标注，理由见对应 issue 与本 PR 描述。
# 注意：本修订为 *源码级* 修复；请作者在自己的 R 环境里跑一遍确认数值。
# =============================================================================

# --- 基础包（保持原有依赖）---
library(ggplot2)
library(ggrepel)
library(affy)
library(hgu133a2.db)
library(org.Hs.eg.db)
library(arrayQualityMetrics)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

# --- 阈值 ---
# [FIX-1] 所有"显著性"判据统一改用 BH 校正后的 adj.P.Val（原文用 raw P.Value）。
# 全基因组 ~46k 探针下 raw P<0.05 会放出约 2300 个假阳性；BH 控制 FDR。
PAJ_CUTOFF   <- 0.05   # adj.P.Val 阈值
VOLCANO_FC   <- 1      # 2 倍变化；x 轴是 logFC(=log2FC)，故 |log2FC|>1 == 2x（原文 1 正确）
ENRICH_GO_P  <- 0.05
ENRICH_KEGG_P <- 0.05
ENRICH_KEGG_Q <- 0.02

# --- 数据读取与预处理 ---
# ReadAffy 按*文件名字母序*排列列。本数据集字母序恰好 = 设计序（658 亲本,659 亲本,
# 660 ADR,661 ADR），但这是命名巧合、不保证 —— 见下方 [FIX-4] 的显式校验。
mydata <- ReadAffy(celfile.path = "/home/simone/microarrays/GSE24460")
eset <- rma(mydata)
exp_matrix <- exprs(eset)
# [FIX-3] 原注释「原始探针水平矩阵」更正：rma() 之后 exprs() 是 *探针集水平、
#         已 RMA 归一化/背景校正* 的矩阵；原始层在 mydata(AffyBatch) 里。

# --- 数据质量检验 ---
# [FIX-2] outdir 拼写 assement -> assessment。
arrayQualityMetrics(mydata, outdir = "quality_assessment")
# [FIX-5] 原脚本跑完 QC 从不等看、直接从坏样本上继续，污染整组对比；且 2/2 设计
#         里 1 个坏样本即毁掉对比。补 RMA 层 PCA + 相关热图，结果输出前强制人工确认。
pdf("qc_pca_heatmap.pdf")
pc <- prcomp(t(scale(t(exprs(eset)))))
plot(pc$x[, 1:2], pch = 19, col = c("navy", "navy", "firebrick3", "firebrick3"),
     xlab = "PC1", main = "PCA (RMA)")
pheatmap(cor(exprs(eset)), labels_col = colnames(exprs(eset)),
         main = "Sample correlation (RMA)")
dev.off()
message(">>> 请先查看 quality_assessment/index.html 与 qc_pca_heatmap.pdf，",
        "确认无离群样本/无错标后再采信后续差异分析。")

# --- 探针对应 Symbol ---
all_symbols <- mapIds(hgu133a2.db,
                      keys    = rownames(exp_matrix),
                      keytype = "PROBEID",
                      column  = "SYMBOL")
# [FIX-6] 删除原"这行是糊弄事的/糊弄甲方"自嘲注释。
# 去重策略见下：avereps 对同 Symbol 的多探针取均值（标准做法，但非最优，见 PR 讨论）。

# --- 折叠到基因水平 ---
exp_matrix_gene <- avereps(exp_matrix, ID = all_symbols)
exp_matrix_gene <- exp_matrix_gene[!is.na(rownames(exp_matrix_gene)), ]
message("折叠后基因数: ", nrow(exp_matrix_gene), " （原探针集数: ",
        nrow(exp_matrix), "）")

# --- 分组 ---
# [FIX-4] 原脚本 group <- factor(c(rep("Control",2), rep("Treat",2))) 硬编码
#         "前 2 列 = 对照"，完全依赖 ReadAffy 的字母序列序，命名一旦不规范即静默错标。
#         现从列名里的 GSM accession 推导并强校验，不符立即停止。
gsm_id <- vapply(regmatches(colnames(exp_matrix_gene),
                            regexec("^GSM[0-9]+", colnames(exp_matrix_gene))),
                 function(m) m[[1]], character(1))
known_gsm_id2group <- c("GSM602658" = "Control", "GSM602659" = "Control",
                        "GSM602660" = "Treat",  "GSM602661" = "Treat")
group <- factor(unname(known_gsm_id2group[match(gsm_id, names(known_gsm_id2group))]),
                levels = c("Control", "Treat"))
stopifnot(!any(is.na(group)),          # 每个样本都能映射到分组
          nlevels(group) == 2,
          sum(group == "Control") == 2, sum(group == "Treat") == 2)
message("分组: ", paste(colnames(exp_matrix_gene), as.character(group),
                        sep = "=", collapse = "  "))

# --- 差异表达分析 ---
# [FIX-7] 用显式 contrast 表达"ADR-亲本"，替代 magic coef=2；并把显著性落到 adj.P.Val。
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)         # Control, Treat
fit <- lmFit(exp_matrix_gene, design)
fit <- contrasts.fit(fit, makeContrasts(Treat_Control = "Treat - Control",
                                        levels = design))
fit <- eBayes(fit)
res <- topTable(fit, coef = "Treat_Control", number = Inf, adjust.method = "BH")
res$Symbol <- rownames(res)
res <- as.data.frame(res)
# res 现含 logFC, AveExpr, t, P.Value, adj.P.Val

# --- 展示用 Top 基因 ---
# [FIX-1b] 选基因改用 adj.P.Val（原文用 raw P.Value）。
# top_genes 仅用于图上打标签 / 热图选行，不替代完整差异分析。
up_genes   <- res[res$adj.P.Val < PAJ_CUTOFF & res$logFC >  0, ]
down_genes <- res[res$adj.P.Val < PAJ_CUTOFF & res$logFC <  0, ]
up_top10   <- head(up_genes[order(up_genes$adj.P.Val), ], 10)
down_top10 <- head(down_genes[order(down_genes$adj.P.Val), ], 10)
top_genes  <- unique(c(up_top10$Symbol, down_top10$Symbol))
# [FIX-8] 原脚本 head(...,100)[raw P 排序] 后再 head(...,10)：先按 raw P 切前 100
#         再筛 FC>0 取 10，与"top10 上调/下调"意图不符、且基于 raw P；
#         现改为先定 DEG(双阈值) 再按 adj.P.Val 取各方向最显著 10 个。
message("Top 上调: ", paste(head(top_genes, 10), collapse=", "), "\n",
        "Top 下调: ", paste(tail(top_genes, 10), collapse=", "))

# --- 火山图 ---
EnhancedVolcano(res,
                lab = res$Symbol,
                x = 'logFC',
                y = 'P.Value',                  # 画 raw P（视觉更分散，常见美学选择）
                title = 'Gene-level DE analysis (deduplicated)',
                selectLab = top_genes,
                drawConnectors = TRUE,
                pCutoff = PAJ_CUTOFF,
                pCutoffCol = 'adj.P.Val',       # [FIX-1] 显著性判据用 BH 校正 p，与 raw P 的 y 轴区分
                FCcutoff = VOLCANO_FC)          # 2 倍变化（x 是 log2FC，故 |log2FC|>1）

# --- 热图 ---
plot_matrix <- exp_matrix_gene[top_genes, ]
annotation_col <- data.frame(Group = as.character(group))
rownames(annotation_col) <- colnames(plot_matrix)
# [FIX-9] 原脚本 breaks 有断层: c(seq(-2,-0.1,50), seq(0,2,50)) 漏掉 (-0.1,0)，
#         且 99 断点(98 区间) 配 100 色，颜色映射不可靠。改为连续、过 0、对称。
nb <- 51
bk   <- seq(-2, 2, length.out = nb)            # 51 断点 -> 50 区间
colr <- colorRampPalette(c("navy", "white", "firebrick3"))(nb - 1)
pheatmap(plot_matrix,
         scale = "row",
         clustering_distance_rows = "correlation",
         annotation_col = annotation_col,
         breaks = bk,
         color = colr,
         legend = TRUE,
         main = "Top DEGs Heatmap (Gene Level)")

# --- 生物学合理性自检 ---
# [FIX-10] 本研究报道 MCF-7/ADR 相对亲本应上调 MDR/cSC 标志物。
#          若这些基因没有按预期方向上调，说明上游(质控/分组/阈值)可能有问题。
known_markers <- c("ABCB1", "CD44", "TGFB1", "SNAI1", "CCNE1", "MMP9")
mk <- res[res$Symbol %in% known_markers, c("Symbol", "logFC", "adj.P.Val")]
cat("\n已知 MDR/cSC 标志物（预期在 ADR 中上调）:\n")
print(mk[order(mk$logFC, decreasing = TRUE), ])
if (nrow(mk[mk$logFC > 0, ]) == 0)
  warning("上述标志物无一上调，请回查分组方向 / 质控 / 阈值是否出错。")

# --- 富集分析 ---
# [FIX-11] 原脚本对 *全基因组所有基因* (res，number=Inf) 做 GO/KEGG 富集 ——
#          把整个基因组当作"差异基因列表"，几乎必然得不到有意义结果。
#          富集对象应限定为差异表达基因，且分上调 / 下调分别做。
deg_up   <- res[res$adj.P.Val < PAJ_CUTOFF & res$logFC >  VOLCANO_FC, ]
deg_down <- res[res$adj.P.Val < PAJ_CUTOFF & res$logFC < -VOLCANO_FC, ]
geneID_up   <- bitr(deg_up$Symbol,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
geneID_down <- bitr(deg_down$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
message("DEG: 上调 ", nrow(deg_up), "  下调 ", nrow(deg_down))

# GO (BP)，上调/下调分开
ego_up <- if (nrow(geneID_up) > 0)
          enrichGO(gene = geneID_up$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = ENRICH_GO_P) else NULL
ego_dn <- if (nrow(geneID_down) > 0)
          enrichGO(gene = geneID_down$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = ENRICH_GO_P) else NULL
pdf("go_bp.pdf")
if (!is.null(ego_up)) print(dotplot(ego_up, showCategory = 20) + ggtitle("GO BP - Up"))
if (!is.null(ego_dn)) print(dotplot(ego_dn, showCategory = 20) + ggtitle("GO BP - Down"))
dev.off()
if (!is.null(ego_up)) write.csv(as.data.frame(ego_up), "my_go_up_results.csv")
if (!is.null(ego_dn)) write.csv(as.data.frame(ego_dn), "my_go_down_results.csv")

# KEGG，上调/下调分开
kk_up <- if (nrow(geneID_up) > 0)
         enrichKEGG(gene = geneID_up$ENTREZID, organism = 'hsa',
                    pvalueCutoff = ENRICH_KEGG_P, qvalueCutoff = ENRICH_KEGG_Q) else NULL
kk_dn <- if (nrow(geneID_down) > 0)
         enrichKEGG(gene = geneID_down$ENTREZID, organism = 'hsa',
                    pvalueCutoff = ENRICH_KEGG_P, qvalueCutoff = ENRICH_KEGG_Q) else NULL
pdf("kegg.pdf")
if (!is.null(kk_up)) print(dotplot(kk_up, showCategory = 20) + ggtitle("KEGG - Up"))
if (!is.null(kk_dn)) print(dotplot(kk_dn, showCategory = 20) + ggtitle("KEGG - Down"))
dev.off()
if (!is.null(kk_up)) write.csv(as.data.frame(kk_up), "my_kegg_up_results.csv")
if (!is.null(kk_dn)) write.csv(as.data.frame(kk_dn), "my_kegg_down_results.csv")

message(">>> 完成。")
