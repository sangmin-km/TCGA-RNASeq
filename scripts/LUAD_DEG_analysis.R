# ============================================================
# TCGA-LUAD RNA-seq DEG Analysis
# Tumor vs Normal
# ============================================================

library(TCGAbiolinks)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(AnnotationDbi)

# 1. 데이터 다운로드
query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# 2. 전처리
data_sub <- data[, data$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
count_matrix <- assay(data_sub, "unstranded")
coldata <- data.frame(
    sample_type = factor(
        ifelse(data_sub$sample_type == "Primary Tumor", "Tumor", "Normal"),
        levels = c("Normal", "Tumor")
    ),
    row.names = colnames(count_matrix)
)

# 3. DESeq2
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = coldata,
    design    = ~ sample_type
)
keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
dds  <- DESeq(dds)

res    <- results(dds, contrast = c("sample_type", "Tumor", "Normal"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# DEG 라벨링
res_df$status <- "Not Significant"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
res_df <- res_df[!is.na(res_df$padj), ]

# 4. Gene Symbol / Entrez ID 변환
res_df$ensembl_clean <- gsub("\\..*", "", res_df$gene)
res_df$symbol <- mapIds(org.Hs.eg.db, keys=res_df$ensembl_clean,
                        column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_df$entrez <- mapIds(org.Hs.eg.db, keys=res_df$ensembl_clean,
                        column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# 5. Volcano plot
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=status)) +
    geom_point(alpha=0.4, size=0.8) +
    scale_color_manual(values=c("Up"="#E64B35","Down"="#4DBBD5","Not Significant"="grey70")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", linewidth=0.3) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.3) +
    labs(title="TCGA-LUAD: Tumor vs Normal", x="log2 Fold Change",
         y="-log10 adjusted p-value", color="Status") +
    theme_bw(base_size=13)
ggsave("results/volcano_LUAD.png", width=8, height=6, dpi=300)

# 6. GO Enrichment
up_genes   <- res_df$entrez[res_df$status=="Up"   & !is.na(res_df$entrez)]
down_genes <- res_df$entrez[res_df$status=="Down" & !is.na(res_df$entrez)]

go_up <- enrichGO(gene=up_genes, OrgDb=org.Hs.eg.db, ont="BP",
                  pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
go_down <- enrichGO(gene=down_genes, OrgDb=org.Hs.eg.db, ont="BP",
                    pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)

dotplot(go_up, showCategory=20, title="GO Enrichment - Up") +
    theme(axis.text.y=element_text(size=8))
ggsave("results/go_up.png", width=8, height=8, dpi=300)

dotplot(go_down, showCategory=20, title="GO Enrichment - Down") +
    theme(axis.text.y=element_text(size=8))
ggsave("results/go_down.png", width=8, height=8, dpi=300)

# 7. 결과 저장
saveRDS(res_df, file="data/res_df.rds")
write.csv(res_df, file="data/DEG_results.csv", row.names=FALSE)
cat("분석 완료!\n")
