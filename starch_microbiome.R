#Downloading Packages from CRAN
install.packages("ggplot2")
install.packages("vegan")
install.packages("igraph")
install.packages("corrr")
install.packages("pheatmap")
install.packages("tidyverse")

#downloading packages from Bioconductor
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")

#load packages
library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(igraph)
library(corrr)
library(pheatmap)
library(tidyverse)

#setting current directory
getwd()
setwd("C:/Users/Asus/Desktop/starch_microbiome/")

#setup output folder
outdir = "results/"
dir.create(outdir, showWarnings = FALSE)
dirs = c("plots","tables","asv","taxa","phyloseq","network")
sapply(file.path(outdir, dirs), dir.create, showWarnings = FALSE, recursive = TRUE)

#list input files
path = "./fastq_files/"
r1 = sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
r2 = sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))
sample.names = sapply(strsplit(basename(r1), "_"), `[`,1)
print(r1)
print(r2)
print(sample.names)

#quality profiles
r1_qc = plotQualityProfile(r1[1:10])
r2_qc = plotQualityProfile(r2[1:10])
print(r1_qc)
print(r2_qc)
ggsave("results/plots/r1_quality_check.png", r1_qc, width = 10, height = 6)
ggsave("results/plots/r2_quality_check.png", r2_qc, width = 10, height = 6)

#filter and trimming
filt_r1 = file.path(path, "filtered", paste0(sample.names, "_r1_filt.fastq.gz"))
filt_r2 = file.path(path, "filtered", paste0(sample.names, "_r2_filt.fastq.gz"))
names(filt_r1) = sample.names
names(filt_r2) = sample.names
filter.out = filterAndTrim(r1, filt_r1, r2, filt_r2,truncLen = c(240,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
print(filter.out)
write.csv(filter.out, "results/tables/filter_stat.csv")

#learn Error rates
err_r1 = learnErrors(filt_r1, multithread = TRUE)
err_r2 = learnErrors(filt_r2, multithread = TRUE)
err_r1_plot = plotErrors(err_r1, nominalQ = TRUE)
err_r2_plot = plotErrors(err_r2, nominalQ = TRUE)
print(err_r1_plot)
print(err_r2_plot)
ggsave("results/plots/r1_error_rate.png", err_r1_plot, width = 10, height = 6)
ggsave("results/plots/r2_error_rate.png", err_r2_plot, width = 10, height = 6)

#dereplication
derep_r1 = derepFastq(filt_r1); names(derep_r1) <- sample.names
derep_r2 = derepFastq(filt_r2); names(derep_r2) <- sample.names
saveRDS(derep_r1, "results/derep_r1.rds")
saveRDS(derep_r2, "results/derep_r2.rds")

#sample inference
dada_r1 = dada(derep_r1, err = err_r1, multithread = TRUE)
dada_r2 = dada(derep_r2, err = err_r2, multithread = TRUE)
print(dada_r1[[1]])
print(dada_r2[[1]])

#merge paired reads
mergers = mergePairs(dada_r1, derep_r1, dada_r2, derep_r2)
head(mergers[[1]])

#make sequence table
seqtab = makeSequenceTable(mergers)
print(dim(seqtab))
table(nchar(getSequences(seqtab)))
saveRDS(seqtab, "results/sequence_table.rds")

#remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method = "pooled", multithread=FALSE)
print(dim(seqtab.nochim))
sum(seqtab.nochim)/sum(seqtab)

#track reads through the pipeline
getN = function(x) sum(getUniques(x))
track = cbind(filter.out, sapply(dada_r1, getN), sapply(dada_r2, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) = c("input", "filtered", "denoised_r1", "denoised_r2", "merged", "nonchim")
rownames(track) = sample.names  
print(track)
write.csv(track, "results/tables/read_track_record.csv")


#assign taxonomy
taxa_species = assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread = TRUE)
taxa_genus = assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread = TRUE)
print(head(taxa_species))
print(head(taxa_genus))
saveRDS(taxa_species, file = "results/species_taxo.rds")
saveRDS(taxa_genus, file = "results/genus_taxo.rds")

taxa.print_species = taxa_species
rownames(taxa.print_species) = NULL
head(taxa.print_species)
write.csv(taxa.print_species, "species_taxo.csv")

taxa.print = taxa_genus
rownames(taxa.print) = NULL
head(taxa.print)
write.csv(taxa.print, "genus_taxo.csv")

taxa_species_names = cbind(ASV = rownames(taxa_species), taxa_species)
write.csv(taxa_species_names, "results/asv/species_taxonomy_with_ASV.csv", row.names = FALSE)

taxa_genus_names = cbind(ASV = rownames(taxa_genus), taxa_genus)
write.csv(taxa_genus_names, "results/asv/genus_taxonomy_with_ASV.csv", row.names = FALSE)

#evaluate accuracy
unqs.mock = seqtab.nochim["SRR13200563",]
unqs.mock = sort(unqs.mock[unqs.mock>0], decreasing = TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in mock communty. \n")

#load metadata
metadata = read.csv("metadata.csv", row.names = 1)
print(head(metadata))

#create a phyloseq object
ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa_species), sample_data(metadata))
print(ps)

#alpha diversity
alpha_plot = plot_richness(ps, measures = c("Shannon","Simpson"))
print(alpha_plot)
ggsave("results/plots/alpha_diversity_plot.png", alpha_plot, width = 8, height = 5)

#Beta Diversity
ordination = ordinate(ps, method = "PCoA", distance = "bray")
beta_plot = plot_ordination(ps, ordination, color = "Group")
print(beta_plot)
ggsave("results/plots/beta_diversity_plot.png", beta_plot, width = 8, height = 5)

#taxonomic barplots
ps.rel = transform_sample_counts(ps, function(x) x/ sum(x))
barplot_phylum = plot_bar(ps.rel, fill = "Phylum")
print(barplot_phylum)
ggsave("results/plots/taxonomic_plot.png", barplot_phylum, width = 10, height = 6)

#differential abundance (DESeq2)
dds = phyloseq_to_deseq2(ps, ~ Treatment)
dds = DESeq(dds)
res = results(dds)
print(head(res))
write.csv(as.data.frame(res), "results/tables/deseq_results.csv")

#SCFA correlation
top_taxa = names(sort(taxa_sums(ps), TRUE)[1:20])
cor_table = cor(as.data.frame(otu_table(ps)[, top_taxa]), metadata[, c("Acetate", "Butyrate", "Formate", "Lactate", "Succinate", "pH")], use = "complete.obs")
heatmap = pheatmap(cor_table)
ggsave("results/plots/heatmap.png", heatmap, width = 12, height = 8, dpi = 300)

#network analysis
otu.rel = transform_sample_counts(ps, function(x) x/sum(x))
cor_matrix = cor(otu_table(otu.rel), method = "spearman")
cor_matrix[abs(cor_matrix) < 0.6] <- 0
net <- graph.adjacency(cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
network = plot(net, vertex.size = 5, vertex_label=NA)
ggsave("results/plots/network_analysis_plot.png", network, width = 12, height = 8, dpi = 300)

rownames


#anxure
# #clear R-Env and plots
# rm(list = ls())
# graphics.off()
# cat("\014")
# r2[1:2]
