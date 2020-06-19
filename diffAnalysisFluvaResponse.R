library(DESeq2)
library(tximport)
library(rtracklayer)

sample.table <- read.csv("../data/sampletable.csv", row.names = 1)

sample.table$files <- as.character(sample.table$files)

files <- sample.table$files

files <- paste0("../data/", files)

names(files) <- rownames(sample.table)


txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## Fix for unexpressed/missing genes:
txi.rsem$length[txi.rsem$length==0] <- 1

dds.sens <- DESeqDataSetFromTximport(txi.rsem, sample.table, ~cell.line+treatment.status)

dds.sens <- dds.sens[,dds.sens$cell.line %in% c("NCI-H929", "KMS-11")]

dds.sens$cell.line <- as.factor(as.character(dds.sens$cell.line))


if(!file.exists("../data/ensembl88.RData")){
  ensembl <- import("ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz")
  save(ensembl, file="../data/ensembl88.RData")
} else {
  load("../data/ensembl88.RData")
}
genes <- ensembl[ensembl$type=="gene"]
names(genes) <- genes$gene_id
mcols(genes) <- mcols(genes)[,c("source", "gene_id", "gene_name", "gene_biotype")]


rownames(dds.sens) <- gsub(rownames(dds.sens), pattern = "\\.[1-9]+", rep="")

genes.to.keep <- intersect(names(genes), rownames(dds.sens))

dds.sens <- dds.sens[genes.to.keep,]

genes <- genes[genes.to.keep,]


rowRanges(dds.sens) <- genes

genes.to.keep.2 <- rownames(dds.sens)[rowData(dds.sens)$gene_biotype == "protein_coding"]

dds.sens <- dds.sens[genes.to.keep.2,]

dds.sens$treatment.status <- relevel(x = dds.sens$treatment.status, ref="untrt")


dds.sens <- DESeq(dds.sens)

res.sens <- results(dds.sens,name="treatment.status_trt_vs_untrt" , alpha = 0.05, cooksCutoff= 4/6)

res.sens <- res.sens[order(res.sens$log2FoldChange),]

res.sens$Symbol <- rowData(dds.sens)$gene_name[match(rownames(res.sens) ,rowData(dds.sens)$gene_id)]


res.sens <- res.sens[complete.cases(res.sens),]

write.csv(res.sens, file="../results/differential_expression_after_treatment_all_genes.csv")


resSig <- subset(res.sens, res.sens$padj < 0.05)
resSig <- resSig[order(resSig$log2FoldChange),]
resSig$Symbol <- rowData(dds.sens)$gene_name[match(rownames(resSig) ,rowData(dds.sens)$gene_id)]
write.csv(resSig, file="../results/differential_expression_after_treatment_sig_genes.csv")
