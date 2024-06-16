rm(list = ls())
library("DESeq2")
library("IHW")
library("pheatmap")
library(dplyr)
library(ggplot2)
library(ggrepel)

title="Homologous_Vaccination_Convalescent_cross_sectional"
directory = "/Users/josephyun/Desktop/AI/Data_Analysis/LeeHyeKyung/"
subdirectory = paste0(directory,title)
setwd(directory)
dir.create(title)
setwd(title)

adjusted_p_value_cutoff=0.05

createSampleTable <- function(baseNames, day, max_number, filePath = "../../data/", outputFileName = "sample_table.csv") {
  # Create an empty data frame to store the table entries
  sampleTable <- data.frame(sampleName = character(),
                            fileName = character(),
                            condition = character(),
                            stringsAsFactors = FALSE)
  
  # Extract IDs from the base names
  baseIds <- as.numeric(sub(".*_(\\d+)$", "\\1", baseNames))
  
  # Create entries for the specified base names
  for (baseName in baseNames) {
    id <- sub(".*_(\\d+)$", "\\1", baseName)  # Extract ID from base name
    entry <- data.frame(
      sampleName = paste("sample", id, paste0("day", day), sep = "_"),
      fileName = paste(filePath, "clean-", baseName, "_Day", day, "_R1_001.fastq-paired.fastq_sorted.txt", sep = ""),
      condition = "E",
      stringsAsFactors = FALSE
    )
    # Combine the new entry with the existing table
    sampleTable <- rbind(sampleTable, entry)
  }
  
  # Create control entries for the remaining samples
  for (i in 1:max_number) {
    if (!(i %in% baseIds)) {
      entry <- data.frame(
        sampleName = paste("sample", i, paste0("day", day), sep = "_"),
        fileName = paste(filePath, "clean-BNT_Convalescent_", i, "_Day", day, "_R1_001.fastq-paired.fastq_sorted.txt", sep = ""),
        condition = "C",
        stringsAsFactors = FALSE
      )
      # Combine the new entry with the existing table
      sampleTable <- rbind(sampleTable, entry)
    }
  }
  
  # Write the table to a CSV file
  write.csv(sampleTable, file = outputFileName, row.names = FALSE, quote = FALSE)
  
  # Return the created data frame (optional, for verification in interactive use)
  return(sampleTable)
}

####################
#Analysis List#
#1 STAT2:NM_001385110:exon19:c.G1749C:p.M583I lines: 73-300
#2 JAK2:NM_001322204:exon11:c.G1402T:p.V468F lines: 304-530
#3 JAK3:NM_000215:exon24:c.C3217T:p.L1073F lines: 531-769
#4 JAK3:NM_000215:exon16:c.G2164A:p.V722I lines: 770-1005
#5 TYK2:NM_001385199:exon22:c.C3124G:p.P1042A lines: 1006-1240
#6 TYK2:NM_001385199:exon14:c.T1865G:p.I622S lines: 1241-1477
#7 TYK2:NM_001385197:exon8:c.G1087A:p.G363S lines: 1478-1712
#8 TYK2:NM_001385197:exon8:c.G1084T:p.V362F lines: 1713-1935
#9 TYK2:NM_001385197:exon3:c.G157A:p.A53T lines: 1936-2169
####################


####################
# Analysis1a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_4", "BNT_Convalescent_15", "BNT_Convalescent_16")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis1"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis1a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis1b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis1b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("TLR7", "FAM26F", "AGER", "AHR")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}


####################
# Analysis2a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_9")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis2"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis2a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis2b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis2b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("PSMB10", "CD7", "BAX", "PRKDC", "OSBP", "LINC00664")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}


####################
# Analysis3a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_13")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis3"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis3a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis3b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis3b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("AGER", "AREG", "BTLA", "APPL1", "A2M", "CASP6")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}





####################
# Analysis4a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_2")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis4"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis4a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis4b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis4b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("A2M", "AATBC", "ACVRL1", "ANXA1", "ADM", "AKIRIN2")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}






####################
# Analysis5a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_2","BNT_Convalescent_4","BNT_Convalescent_5")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis5"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis5a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis5b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis5b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("FCHO1", "MAP3K15", "TNFAIP3", "RGL1", "RAB39A")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}







####################
# Analysis6a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_15","BNT_Convalescent_4","BNT_Convalescent_5","BNT_Convalescent_7")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis6"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis6a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis6b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis6b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("LOC101927243", "IL20RB", "RAB39A", "SOCS1", "KCNK7", "B3GNT8", "UNC13D", "G0S2")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}


####################
# Analysis7a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_10","BNT_Convalescent_11","BNT_Convalescent_14","BNT_Convalescent_16")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis7"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis7a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis7b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis7b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("HLA-DRB5", "RIMBP3", "HGS", "SFPQ", "CD38", "MAVS")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}



####################
# Analysis8a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_16","BNT_Convalescent_2","BNT_Convalescent_7","BNT_Convalescent_8","BNT_Convalescent_9")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis8"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis8a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis8b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis8b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("RMRP", "MRS2P2", "GALK1", "FUT7", "AREG", "PRLR", "RAG1")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}



####################
# Analysis9a_Day1_SNP_vs_Day1_C
####################

title="Homologous_Vaccination_Convalescent_cross_sectional"
# Num samples
max_number <- 16
baseNames <- c("BNT_Convalescent_2")
vacinationType="BNT_Convalescent"

setwd(subdirectory)
analysisTitle="Analysis9"
analysisDirectory=paste0(subdirectory,"/",analysisTitle)
dir.create(analysisDirectory)
setwd(analysisDirectory)


subTitle="Analysis9a_Day1_SNP_vs_Day1_C"
analysisDirectoryA=paste0(analysisDirectory,"/",subTitle)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=1, max_number)

# Analysis
this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
dds$condition <- relevel(dds$condition, ref="C")
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

# Get mean and median.
normalized_counts_df <- as.data.frame(normalized_counts)

samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Create a data frame to format the results
results_a <- data.frame(
  normalized_counts_df,
  Day1_Mean_E = mean_e,
  Day1_Mean_C = mean_c,
  Day1_Median_E = median_e,
  Day1_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_a_filtered <- results_a[rownames(results_a) %in% immune_genes, ]
write.csv(results_a_filtered, "filtered_counts.csv")

####################
# Analysis9b_Day1_SNP_vs_Day0_SNP
####################

subTitle="Analysis9b_Day0_SNP_vs_Day0_C"
analysisDirectoryB=paste0(analysisDirectory,"/",subTitle)

setwd(analysisDirectory)
dir.create(subTitle)
setwd(subTitle)

createSampleTable(baseNames, day=0, max_number)

this_sampleTable <- read.csv(file="sample_table.csv", header=TRUE)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = this_sampleTable, directory = paste0(analysisDirectory), design = ~ condition)
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


# Calculate mean and median
normalized_counts_df <- as.data.frame(normalized_counts)
samples_e <- this_sampleTable$sampleName[this_sampleTable$condition == "E"]
samples_c <- this_sampleTable$sampleName[this_sampleTable$condition == "C"]

data_e <- select(normalized_counts_df, samples_e)
data_c <- select(normalized_counts_df, samples_c)
mean_e <- apply(data_e, 1, mean)
median_e <- apply(data_e, 1, median)
mean_c <- apply(data_c, 1, mean)
median_c <- apply(data_c, 1, median)

# Optionally, create a data frame to format the results
results_b <- data.frame(
  normalized_counts_df,
  Day0_Mean_E = mean_e,
  Day0_Mean_C = mean_c,
  Day0_Median_E = median_e,
  Day0_Median_C = median_c,
  res
)


immune_gene_list <- read.csv("../../../data/immune_gene_list.csv")
immune_genes <- immune_gene_list[, 1]

results_b_filtered <- results_b[rownames(results_b) %in% immune_genes, ]
print(head(results_b_filtered))

write.csv(results_b_filtered, "filtered_counts.csv")


# Create final filtered output
results_b_filtered$gene <- rownames(results_b_filtered)
results_a_filtered$gene <- rownames(results_a_filtered)

final_df <- inner_join(results_b_filtered, results_a_filtered, by = "gene")
rownames(final_df) <- final_df$gene
final_df$gene <- NULL

print(final_df)

write.csv(final_df, paste0("../", vacinationType, "_Normalized_Counts.csv"))

####################
# Heatmap and Plots #
####################

subTitle="Plots"
setwd("../")
dir.create(subTitle)
setwd(subTitle)

final_df$composite_value <- final_df$padj.x - final_df$padj.y
df_sorted <- final_df[order(final_df$composite_value, decreasing = TRUE), ]

write.csv(df_sorted[1:30,], "BNT_Convalescent_top_30_genes_sorted_deregulated.csv")

df_subset <- df_sorted[1:50, c("Day0_Mean_C", "Day1_Mean_C", "Day0_Mean_E", "Day1_Mean_E")]
df_subset$MeanC_Difference <- df_subset$Day0_Mean_C - df_subset$Day1_Mean_C
df_subset$MeanE_Difference <- df_subset$Day0_Mean_E - df_subset$Day1_Mean_E

signed_log2 <- function(x) {
  sign(x) * log2(1 + abs(x))
}
df_matrix <- as.matrix(df_subset[, c("MeanC_Difference", "MeanE_Difference")])
df_transformed <- apply(df_matrix, 2, signed_log2)

my_filename = "heatmap_top50.png"
png(filename = my_filename, width=5, height=5, unit = "in", res=300)

pheatmap(df_transformed, cluster_row=F, show_rownames = T, show_colnames = T, cluster_cols = F, cluster_rows = T, fontsize_row = 4)

dev.off()

# Plots
this_sampleTableA <- read.csv(paste0(analysisDirectoryA,"/sample_table.csv"))
this_sampleTableB <- read.csv(paste0(analysisDirectoryB,"/sample_table.csv"))

day1_C <- this_sampleTableA$sampleName[this_sampleTableA$condition == "C"]
day1_E <- this_sampleTableA$sampleName[this_sampleTableA$condition == "E"]
day0_C <- this_sampleTableB$sampleName[this_sampleTableB$condition == "C"]
day0_E <- this_sampleTableB$sampleName[this_sampleTableB$condition == "E"]

split_samples <- strsplit(c(day1_C, day1_E, day0_C, day0_E), "_")
sample_list <- sapply(split_samples, function(x) {
  paste(x[1], x[2], sep = "_")  
})

# Manually entered
n_range <- c("AATBC", "ACVRL1", "A2M", "ANXA1", "ADM", "AKIRIN2", "ARID5A")
for (n in n_range) {
  df <- df_sorted[n, c(day1_C, day1_E, day0_C, day0_E)]
  value_list <- as.vector(t(df)) 
  
  data <- data.frame(
    Sample = sample_list,
    Sample_Num = as.numeric(sub("sample_", "", sample_list)),
    Timepoint = c(rep("Day1", length(day1_C)), rep("Day1", length(day1_E)), rep("Day0", length(day0_C)), rep("Day0", length(day0_E))),  
    Value = value_list,  
    Condition = c(rep("Control", length(day1_C)), rep("SNP", length(day1_E)), rep("Control", length(day0_C)), rep("SNP", length(day0_E)))
  )
  
  my_filename = paste0(n, "_plot.png")
  title = paste0(n, " plot")
  png(filename = my_filename, width=5, height=5, unit = "in", res=300)
  
  print(ggplot(data, aes(x = Timepoint, y = Value, group = Sample, color = Timepoint)) +
          geom_point(size = 3) + 
          geom_line() + 
          geom_text_repel(aes(label=Sample_Num), color="black", data=subset(data, Timepoint=="Day0"), max.overlaps=20, box.padding = 0.1, point.padding = 0.5, nudge_x = -0.5, segment.color = alpha("blue", 0)) + 
          scale_color_manual(values = c("Day0" = "black", "Day1" = "red")) + 
          facet_wrap(~ Condition) + 
          labs(title = title, y = "Value", x = "Timepoint") +
          theme_minimal())
  
  dev.off()
  if (file.exists(my_filename)) {
    print(paste("Successfully generated plot:", my_filename))
  } else {
    print(paste("Failed to generate plot:", my_filename))
  }
}

