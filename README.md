NGS RNA-seq Data Analysis Pipeline
UVB vs Control (8 × 8 Biological Replicates)
Overview
This repository provides a fully reproducible, end-to-end RNA-sequencing (RNA-seq) analysis pipeline for bulk transcriptomics. The workflow combines:
Bash-based preprocessing for raw sequencing data
R-based differential expression analysis using DESeq2
Quality control and visualization, including PCA and volcano plots
The pipeline is demonstrated on 16 samples: UVB-exposed (n=8) vs Control (n=8), producing publication-ready outputs.
Pipeline Workflow
Raw FASTQ (paired-end)
        ↓
FastQC — Quality assessment
        ↓
Trim Galore — Adapter & quality trimming
        ↓
HISAT2 — Genome alignment
        ↓
SAMtools — BAM sorting & indexing
        ↓
HTSeq — Gene-level read counting
        ↓
Count matrix assembly
        ↓
DESeq2 — Differential expression analysis
        ↓
VST PCA & Volcano Plots — Sample clustering and visualization
        ↓
MultiQC — Aggregate QC summary
Detailed Description of Tools
1. FastQC
Purpose: Assess the quality of raw sequencing reads
Checks: Per-base quality scores, GC content, adapter contamination, sequence duplication
Input: Raw FASTQ files
Output: HTML and TXT reports (results/fastqc/)
Command example:
fastqc data/raw_fastq/*.fastq.gz -o results/fastqc/
2. Trim Galore
Purpose: Remove adapter sequences and low-quality bases
Based on: Cutadapt and FastQC
Input: Raw FASTQ files
Output: Trimmed FASTQ files (data/trimmed_fastq/)
Command example:
trim_galore --paired data/raw_fastq/CTRL_Rep01_R1.fastq.gz data/raw_fastq/CTRL_Rep01_R2.fastq.gz -o data/trimmed_fastq/
3. HISAT2
Purpose: Splice-aware alignment of RNA-seq reads to reference genome
Input: Trimmed FASTQ files
Output: SAM files, later converted to BAM (results/alignment/)
Command example:
hisat2 -x reference/genome_index -1 data/trimmed_fastq/CTRL_Rep01_R1_val_1.fq.gz -2 data/trimmed_fastq/CTRL_Rep01_R2_val_2.fq.gz -S results/alignment/CTRL_Rep01.sam
4. SAMtools
Purpose: Convert SAM → BAM, sort and index for downstream analysis
Commands example:
samtools view -bS results/alignment/CTRL_Rep01.sam > results/alignment/CTRL_Rep01.bam
samtools sort results/alignment/CTRL_Rep01.bam -o results/alignment/CTRL_Rep01_sorted.bam
samtools index results/alignment/CTRL_Rep01_sorted.bam
5. HTSeq
Purpose: Count reads mapping to annotated genes
Input: Sorted BAM files + GTF annotation
Output: Gene-level count files (results/counts/)
Command example:
htseq-count -f bam -r name -s no -i gene_id results/alignment/CTRL_Rep01_sorted.bam reference/genes.gtf > results/counts/CTRL_Rep01_counts.txt
6. MultiQC
Purpose: Aggregate QC reports from FastQC, Trim Galore, HISAT2 mapping rates
Output: Single HTML report (results/multiqc/multiqc_report.html)
Command example:
multiqc results/ -o results/multiqc/
Experimental Design
Factor	Description
Study type	Bulk RNA-seq
Comparison	UVB vs Control
Replicates	Control (n=8), UVB (n=8)
Total samples	16
Sequencing	Illumina paired-end
DESeq2-Ready Metadata
Sample Metadata (metadata/samples.tsv)
sample_id	condition	batch
CTRL_Rep01	Control	B1
CTRL_Rep02	Control	B1
CTRL_Rep03	Control	B1
CTRL_Rep04	Control	B1
CTRL_Rep05	Control	B2
CTRL_Rep06	Control	B2
CTRL_Rep07	Control	B2
CTRL_Rep08	Control	B2
UVB_Rep01	UVB	B1
UVB_Rep02	UVB	B1
UVB_Rep03	UVB	B1
UVB_Rep04	UVB	B1
UVB_Rep05	UVB	B2
UVB_Rep06	UVB	B2
UVB_Rep07	UVB	B2
UVB_Rep08	UVB	B2
condition = experimental group
batch = optional covariate for batch-aware analysis
DESeq2 Design Formula
# Standard design
design = ~ condition

# Batch-aware design
design = ~ batch + condition
DESeq2 Analysis, PCA & Volcano Plot
library(DESeq2)
library(ggplot2)

# Load count matrix
counts <- read.table("results/counts/counts_matrix.tsv", header=TRUE, row.names=1, sep="\t")
coldata <- read.table("metadata/samples.tsv", header=TRUE, row.names=1, sep="\t")
counts <- counts[, rownames(coldata)]

# DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~ batch + condition)
dds$condition <- relevel(dds$condition, ref="Control")
dds <- DESeq(dds)

# Differential expression
res <- results(dds, contrast=c("condition","UVB","Control"))
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "results/DESeq2/DESeq2_UVB_vs_Control_results.csv")

# VST & PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color=significant), alpha=0.6) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  labs(title="Volcano Plot: UVB vs Control",
       x="Log2 Fold Change", y="-log10 Adjusted P-value")
GEO / ENA Submission Guidance
Files to include: counts_matrix.tsv, samples.tsv, DESeq2 results, multiqc_report.html, experimental description
Metadata example:
Study design: Bulk RNA-seq
Organism: Homo sapiens
Comparison: UVB vs Control
Replicates: 8 per group
Sequencing: Illumina paired-end
Analysis: HISAT2 + HTSeq + DESeq2
Best-Practice Thresholds
Adjusted p-value: padj < 0.05
Log2 fold change: |log2FC| ≥ 1
Optionally filter low-expression genes before DESeq2
Author
Maruf Hasan
Postdoctoral Researcher
Expertise: Transcriptomics | Bioinformatics | UV-mediated Signaling | Vitamin D Biology | Translational
