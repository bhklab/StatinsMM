
library(SummarizedExperiment)

library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(rtracklayer)
options(warn=-1)


if(!file.exists("../data/gencodev26.RData")) {
  gencodev26.full <- import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz")
  gencodev26 <- gencodev26.full[gencodev26.full$type=="gene",]
  save(gencodev26, file="../data/gencodev26.RData")
} else {
  load("../data/gencodev26.RData")
}

mcols(gencodev26) <- mcols(gencodev26)[,c("gene_id", "gene_type", "gene_name", "level")]


sample.table <- read.delim("../data/samplesheet.tsv", row.names = NULL)
rownames(sample.table) <- sample.table$Sample.Name

sample.table$File.Path <- as.character(sample.table$File.Path)
files <- sample.table$File.Path

names(files) <- rownames(sample.table)

files <- paste0("../data/", files)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)


## Fix for unexpressed/missing genes:
txi.rsem$length[txi.rsem$length==0] <- 1
dds <- DESeqDataSetFromTximport(txi.rsem, sample.table, ~Group)


gencodev26 <- gencodev26[gencodev26$gene_id %in% rownames(dds),]

gencodev26 <- gencodev26[match(rownames(dds), gencodev26$gene_id),]

rowData(dds) <- gencodev26

dds <- dds[rowData(dds)$gene_type == "protein_coding",]

dds <- DESeq(dds, minReplicatesForReplace=Inf)


MM.Cells <- c("KMS11", "OPM2", "KHM11", "NCIH929", "OCIMY7", "RPMI8226", "KMM1", "EJM", "LP1", "OCIMY5", "SKMM1", "ANBL62", "U266", "JJN3")
MM.Cells.cat <- as.factor(c("sens", "sens", "sens", "sens", "sens", "sens", "sens", "res", "res", "res", "res", "res", "res", "res"))

Four_Fourteen_Status <- c(1,1,1,1,2,2,2,2,1,2,2,2,2,2)
Four_Fourteen_Status <- c("Yes", "No")[as.integer(Four_Fourteen_Status)]
tbl <- data.frame("Cells" = MM.Cells, "Category" = MM.Cells.cat, "Four_Fourteen_Status" = Four_Fourteen_Status)
print(tbl)

## Top 30 genes less expressed in Sensitive vs Resistant cells

# rule of thumb for cooks distance cutoff: 4/(n-k-1) = 0.6667
res <- DESeq2::results(dds, alpha = 0.05, cooksCutoff= 4/11, pAdjustMethod="bonferroni")

res$Symbol <- rowData(dds)$gene_name[match(rownames(res) ,rowData(dds)$gene_id)]
resSig <- subset(res, res$padj < 0.05)

resSig <- resSig[order(resSig$log2FoldChange),]
#rownames(resSig) <- NULL
write.csv(as.data.frame(resSig), file="../results/rsem_bonferroni_diff_exp_genes.csv")
res <- res[order(res$log2FoldChange),]

write.csv(as.data.frame(res), file="../results/rsem_all_diff_exp_genes.csv")

resSig <- resSig[order(resSig$log2FoldChange),]
print(as.data.frame(resSig[,c(1,2,6,7)]))

library(EnhancedVolcano)

pdf("../results/RSEM_volcano_plot.pdf", height = 10, width = 10)
EnhancedVolcano(res, res$Symbol, x = "log2FoldChange", y="pvalue", title=NULL, xlim = c(-12,12), pCutoff = 4.67e-6, FCcutoff = 2, cutoffLineType = "blank", transcriptLabSize = 10, transcriptLabvjust = -0.3, subtitle = NULL, transcriptPointSize = 2.5, axisLabSize = 26, legendPosition = "bottom", legendLabSize = 26, captionLabSize = 20, ylim=c(0,10))
dev.off() 




