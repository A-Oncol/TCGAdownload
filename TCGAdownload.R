######目的：TCGAbiolinks下载TCGA泛癌多组学数据
######作者：申奥
######日期：2023-09-14
######版本：R 4.3.1


# 0.环境设置----
rm(list = ls())
options(stringAsFactors = F)

library(TCGAbiolinks)
library(openxlsx)
library(getopt)

spec <- matrix(
  # 每行五个，第五个可选，也就是说第五列可以不写
  # byrow 按行填充矩阵的元素
  # ncol  每行填充五个元素
  c("cancerType",  "c", 2, "character", "Input the cancer type"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)


# 1.设置癌种名称----
projects <- getGDCprojects()$project_id
projects
tcga_projects <- projects[grepl("^TCGA", projects, perl = T)]
tcga_projects
cancer_type <- opt$cancerType


# 2.下载并整理表达谱数据----
## 2.1 下载----
query <- GDCquery(project = cancer_type,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query = query, directory = "GDCdata", files.per.chunk = 5)
GDCprepare(query = query, save = T, save.filename = paste0("data/RNA/", cancer_type, "_RNA-seq.RData"))

## 2.2 提取mRNA和lncRNA的count和TPM数据----
library(SummarizedExperiment)
load(paste0("data/RNA/", cancer_type, "_RNA-seq.RData"))
se <- data
se

names(assays(se))
rowdata <- rowData(se)
names(rowdata)
table(rowdata$gene_type)
head(rowdata$gene_name)
length(rowdata$gene_name)

se_mRNA <- se[rowdata$gene_type == "protein_coding",]
se_lnc <- se[rowdata$gene_type == "lncRNA",]
se_mRNA
se_lnc

mRNA_count <- assay(se_mRNA, "unstranded")
mRNA_tpm <- assay(se_mRNA, "tpm_unstrand")
lnc_count <- assay(se_lnc, "unstranded")
lnc_tpm <- assay(se_lnc, "tpm_unstrand")

mRNA_symbol <- rowData(se_mRNA)$gene_name
head(mRNA_symbol)
lnc_symbol <- rowData(se_lnc)$gene_name
head(lnc_symbol)

mRNA_count_annoted <- cbind(data.frame(mRNA_symbol),
                            as.data.frame(mRNA_count))
mRNA_tpm_annoted <- cbind(data.frame(mRNA_symbol),
                          as.data.frame(mRNA_tpm))
lnc_count_annoted <- cbind(data.frame(lnc_symbol),
                           as.data.frame(lnc_count))
lnc_tpm_annoted <- cbind(data.frame(lnc_symbol),
                         as.data.frame(lnc_tpm))
save(mRNA_count_annoted, file = paste0("data/RNA/", cancer_type, "_mRNA_count.RData"))
save(mRNA_tpm_annoted, file = paste0("data/RNA/", cancer_type, "_mRNA_tpm.RData"))
save(lnc_count_annoted, file = paste0("data/RNA/", cancer_type, "_lncRNA_count.RData"))
save(lnc_tpm_annoted, file = paste0("data/RNA/", cancer_type, "_lncRNA_tpm.RData"))


# 3.下载CNV数据----
query <- GDCquery(project = cancer_type,
                  data.category = "Copy Number Variation",
                  data.type = "Masked Copy Number Segment",
                  access = "open")
GDCdownload(query = query, directory = "GDCdata", files.per.chunk = 2)
GDCprepare(query = query, save = T, save.filename = paste0("data/CNV/", cancer_type, "_CNV.RData"))


# 4.下载突变数据----
query <- GDCquery(project = cancer_type,
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation",
                  access = "open")
GDCdownload(query = query, directory = "GDCdata", files.per.chunk = 5)
GDCprepare(query = query, save = T, save.filename = paste0("data/Mutation/", cancer_type, "_Mutation.RData"))


# 5.下载miRNA数据----
query <- GDCquery(project = cancer_type,
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")
GDCdownload(query = query, directory = "GDCdata", files.per.chunk = 2)
GDCprepare(query = query, save = T, save.filename = paste0("data/RNA/", cancer_type, "_miRNA.RData"))


# 6.下载临床数据----
clinic_query <- GDCquery(project = cancer_type,
                         data.category = "Clinical", 
                         data.format = "bcr xml")
GDCdownload(clinic_query)
clinical <- GDCprepare_clinic(clinic_query, clinical.info = "patient")
save(clinical, file = paste0("data/Clinical/", cancer_type, "_clinical.RData"))
write.xlsx(clinical, file = paste0("data/Clinical/", cancer_type, "_clinical.xlsx"))


# 7.下载甲基化数据----
query <- GDCquery(project = cancer_type,
                   data.category = "DNA Methylation",
                   data.type = "Methylation Beta Value",
                   platform = "Illumina Human Methylation 450")
GDCdownload(query = query, directory = "GDCdata", files.per.chunk = 2)
GDCprepare(query = query, save = T, save.filename = paste0("data/Methylation/", cancer_type, "_Meth.RData"))