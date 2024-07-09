rm(list = ls())
library("DESeq2")
library("IHW")
library("pheatmap")

############################
# Create sample tables#
############################

setwd("XXXXX/sample_tables/")
system("cat sample_table_Analysis2ap_Day1_SNP_vs_Day0_SNP.csv | sed -e 's/day1/day7/' -e 's/Day1/Day7/' > sample_table_Analysis2bp_Day7_SNP_vs_Day0_SNP.csv")
system("cat sample_table_Analysis2ap_Day1_SNP_vs_Day0_SNP.csv | sed -e 's/day1/month4/' -e 's/Day1/5/' > sample_table_Analysis2cp_Month4_SNP_vs_Day0_SNP.csv")

# Note: for day7, one file is missing. So after copy, delete sample 10.
# system("cat sample_table_Analysis2dp_Day1_C_vs_Day0_C.csv | sed -e 's/day1/day7/' -e 's/Day1/Day7/' > sample_table_Analysis2ep_Day7_C_vs_Day0_C.csv")
system("cat sample_table_Analysis2dp_Day1_C_vs_Day0_C.csv | sed -e 's/day1/month4/' -e 's/Day1/5/' > sample_table_Analysis2fp_Month4_C_vs_Day0_C.csv")

############################
#Analysis2ap_Day1_SNP_vs_Day0_SNP####
############################

title="Analysis2ap_Day1_SNP_vs_Day0_SNP"
directory = "XXXXX"
adjusted_p_value_cutoff=0.05

setwd(directory)
dir.create(title)
setwd(title)
# Already prepared sample table.
file1 = paste0("../../sample_tables/sample_table_", title, ".csv")
file2 = "sample_table.csv"
file.copy(file1, file2)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(directory,title), design = ~ ID+condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C") #Compare compareA against the compareB as a reference
dds <- DESeq(dds)

res <- results(dds, filterFun=ihw, alpha=adjusted_p_value_cutoff)
resOrdered <- res[order(res$padj),]
significant_genes <- rownames(res[!is.na(res$padj) & res$padj < adjusted_p_value_cutoff, ])

summary(resOrdered, alpha = adjusted_p_value_cutoff)

#1: Fold_change_and_p_values_all_genes.csv
my_filename = paste0("RNA_seq_results_Fold_change_and_p_values_", title, "_all_genes.csv")
write.csv(as.data.frame(resOrdered), file=my_filename)

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts.csv")


####################################
#Analysis2bp_Day7_SNP_vs_Day0_SNP####
####################################

title="Analysis2bp_Day7_SNP_vs_Day0_SNP"
directory = "XXXXXX"
adjusted_p_value_cutoff=0.05

setwd(directory)
dir.create(title)
setwd(title)
# Already prepared sample table.
file1 = paste0("../../sample_tables/sample_table_", title, ".csv")
file2 = "sample_table.csv"
file.copy(file1, file2)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(directory,title), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C") #Compare compareA against the compareB as a reference
dds <- DESeq(dds)

res <- results(dds, filterFun=ihw, alpha=adjusted_p_value_cutoff)
resOrdered <- res[order(res$padj),]
significant_genes <- rownames(res[!is.na(res$padj) & res$padj < adjusted_p_value_cutoff, ])

summary(resOrdered, alpha = adjusted_p_value_cutoff)

#1: Fold_change_and_p_values_all_genes.csv
my_filename = paste0("RNA_seq_results_Fold_change_and_p_values_", title, "_all_genes.csv")
write.csv(as.data.frame(resOrdered), file=my_filename)

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts.csv")
