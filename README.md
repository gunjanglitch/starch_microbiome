# 🧬 Gut Microbiome Response to Dietary Starches

**One-line description:**  
Analyzed how human gut microbiota respond to different dietary starches using 16S rRNA sequencing and microbiome profiling.

---

## 📌 Project Overview

This project explores the impact of different dietary starches on human gut microbiota using publicly available 16S rRNA amplicon sequencing data from [PRJNA682640](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA682640).  
The goal is to observe how various fermentable substrates affect microbial diversity, taxonomy, and short-chain fatty acid (SCFA) production in vitro.

---

## 🧪 Dataset Summary

- **Source**: NCBI SRA Project [PRJNA682640](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA682640)
- **Samples I selected**: SRR13200563 | SRR13200564 | SRR13200565 | SRR13200568 | SRR13200569 | SRR13200570 | SRR13200573 | SRR13200574 | SRR13200578 | SRR13200580
- **Treatment**: Human fecal microbiota fermented with:
  - Banana
  - CornStarch
  - Tapioca
  - Amylopectin
  - HAM2, HAM4 (resistant starch)
  - Potato starches (PS_Ba, PS_Rb)
  - Water (control)
- **Measured Outputs**: SCFA levels (acetate, butyrate, lactate, formate, succinate) and pH

---

## 🔧 Analysis Pipeline

The analysis was conducted in **R** using the `DADA2` and `Phyloseq` packages.

1. **Download & organize raw data**
2. **Quality check** with `FastQC`
3. **Read filtering & trimming**
4. **Error rate learning**
5. **Dereplication and sample inference**
6. **Merging paired reads**
7. **Chimera removal**
8. **ASV table creation**
9. **Taxonomic classification** (SILVA database)
10. **Alpha diversity** analysis (Shannon, Simpson indices)
11. **Beta diversity** analysis (PCoA plots)
12. **Taxonomic visualization** (barplots, heatmaps)
13. **Differential abundance analysis**
14. **SCFA–microbiota correlation analysis**

---

## 📈 Key Findings

- **Tapioca** → Highest **butyrate** levels (21.08 mM), but **lowest microbial diversity**
- **Banana & CornStarch** → Supported more **diverse** and **balanced** microbial communities
- **Amylopectin** → Elevated **lactate**, minimal butyrate, suggesting possible microbial imbalance
- **Beta diversity**: Samples clustered distinctly by starch type — strong effect of diet on community structure

---

## 🛠️ Tools & Packages Used

| Tool         | Purpose                            |
|--------------|-------------------------------------|
| FastQC       | Quality control of raw reads        |
| DADA2        | Denoising, ASV generation           |
| Phyloseq     | Taxonomy, diversity, visualization  |
| ggplot2      | Plotting and customization          |
| R Base/Tidyverse | Data wrangling and analysis     |


---

## 📫 Contact

📧 Feel free to connect or reach out:  
🔗 🌐 [LinkedIn](https://www.linkedin.com/in/gunjan-sarode/) | 📫 gunjansarode.bioinfo@gmail.com
🐛 Found a bug or have a question? Open an issue on this repo!

---

> 🚀 This project helped me understand the connection between diet, microbiota, and function, and strengthened my skills in microbiome data analysis.




