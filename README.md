# **RT-qPCR analysis of virus accumulation and defense related gene expression after PGPR treatment**

A fully reproducible workflow for evaluating systemic virus infection, PGPR-mediated mitigation of CMV and PVY, and defense-related gene expression using RT-qPCR in *Arabidopsis thaliana* and *Nicotiana benthamiana*.
This repository provides all datasets, R code, and visualization files necessary to ensure transparency, reproducibility, and statistical rigor across all three research objectives.
This repository supports transparent, organized, and repeatable plant–microbe–virus interaction research. The study was conducted across three independent biological replicates, ensuring statistical robustness and reproducibility.

The study was carried out using three independent biological replicates, with technical replicates included for each RT-qPCR reaction.

All files including raw data, R code, outputs, and figures are available and organized to support repeatable analyses.
 
**Project overview**

This repository contains the complete workflows used to address three major objectives:

**Objective 1:** Systemic infection profiling

To determine whether Cucumber mosaic virus (CMV) and Potato virus Y (PVY) establish systemic infection in *Arabidopsis thaliana* or *Nicotiana benthamiana*. RT-qPCR was used to quantify viral RNA in both inoculated and systemic leaves across multiple timepoint.

**Objective 2:** PGPR-mediated mitigation of virus infection

To evaluate whether pre-treatment with plant-growth promoting rhizobacteria (PGPR): *Pseudomonas fluorescens*, *Serratia marcescens*, *Bacillus subtilis* reduces CMV or PVY accumulation.

**Note:**PGPR treatments were applied one week prior to virus inoculation to ensure time for induced systemic resistance (ISR) activation.

**Objective 3:**Defense gene expression after PGPR treatment and PGPR + Virus challenge

To characterize host immune activation, RT-qPCR was used to measure the expression of salicylic-acid, jasmonic-acid, and ethylene pathway genes under:

PGPR treatment alone

PGPR treatment followed by CMV infection

Genes include *PR1*, *PR4*, *PR5*, *NPR1*, *EDS1*, *PAD4*, *ICS1*, *PDF1.2*, *RdR1*, *OPR3*, and *Thi2.1*.

The full dataset used in this project is included in the repository under the folder named `Data/`.


### Folder and File Structure
- [Data](https://github.com/BibechanaPaudel/Reproducibility-project/tree/main/Data): Contains .csv files for each of the three biological replicates.

- [Nb_PGPR+PVY.Rmd](https://github.com/BibechanaPaudel/Reproducibility-project/blob/main/Nb_PGPR%2BPVY.Rmd): The complete R script with code, output, and narrative.

- [Nb_PGPR+PVY.html](https://github.com/BibechanaPaudel/Reproducibility-project/blob/main/Nb_PGPR%2BPVY.html): Rendered HTML version of the RMarkdown.

- [Nb_PGPR+PVY.md](https://github.com/BibechanaPaudel/Reproducibility-project/blob/main/Nb_PGPR%2BPVY.md): Markdown-exported version for GitHub-friendly viewing.

- [Nb_PGPR+PVY_files/figure-gfm](https://github.com/BibechanaPaudel/Reproducibility-project/tree/main/Nb_PGPR%2BPVY_files/figure-gfm): Output plots in high-quality vector format (.png). 


### Project Overview

**Objective:**  
To evaluate whether PGPR treatments mitigate PVY accumulation, based on qPCR-derived Cq values, analyzed over time (Dpi) in both inoculated and systemic leaf tissues of *N. benthamiana*.

**Treatments:**
- Control  
- *Pseudomonas fluorescens*  
- *Serratia marcescens*  
- *Bacillus amyloliquefaciens*  
- *Bacillus subtilis*

### Script workflow

#### **1. Data Import and Preparation**

- Three biological replicate datasets (`.csv`) are stored in the `Data/` folder.
- All data are processed in a single RMarkdown (`.Rmd`) script.
- Each file contains Cq values from RT-qPCR for PVY quantification, along with metadata (Treatment, Dpi, Technical Replicate).
- Viral load is calculated from Cq values using the regression equation **Y= −3.93X+49.153**, where Y= Cq value and X= copy number of virus from [Feng, J.L et al., 2006](https://academic.oup.com/abbs/article/38/10/669/217), then log-transformed for normality.

#### **2. Data Grouping and Statistical Analysis**

- Data are grouped by Treatment, Days post-inoculation (Dpi), and Technical Replicates.
- A linear model is fit to the log-transformed viral load to analyze interaction effects.
- ANOVA (type-II via car::Anova) is used to test for significant Treatment × Dpi interactions.
- Estimated marginal means are calculated using `emmeans`, followed by Tukey-adjusted pairwise comparisons.
- Significance letters (e.g., a, b, c) are generated using `multcompView`.
- Significance groupings (letters) are extracted for visualisation.

#### **3. Visualisation**

- Bar plots are generated for:
  - **Inoculated leaves** at 1, 4, 7, and 10 Dpi
  - **Systemic leaves** at 7 and 10 Dpi
- Visuals include:
  - Grouped bar charts with error bars (standard error)
  - Statistical significance labels over bars
  - Combined plots with shared legends
- Color-blind friendly palettes are used.

### Citation
If you use this code for your research, please cite it using the Zenodo DOI provided here:

[![DOI](https://zenodo.org/badge/966418246.svg)](https://doi.org/10.5281/zenodo.15258548)

### Disclaimer
The code provided herein has not been peer-reviewed and may contain errors. Users are encouraged to test the code thoroughly and verify its accuracy for their specific applications. The author is not responsible for any errors or inaccuracies in the results generated by this script.

### File Tree

```
├── Data
│   ├── Nb_PGPR+PVY_1st_Reproducibility.csv   ## Data from first biological replication
│   ├── Nb_PGPR+PVY_2nd_Reproducibility.csv   ## Data from second biological replication
│   └── Nb_PGPR+PVY_3rd_Reproducibility.csv   ## Data from third biological replication
├── Nb_PGPR+PVY.html   ##script in html format
├── Nb_PGPR+PVY.md     ##script in github favoured markdowm format
├── Nb_PGPR+PVY.Rmd    ##script in rmd format
├── Nb_PGPR+PVY_files
│   └── figure-gfm
│       ├── Combine fig 1st trial-1.png  ## combine fig from 1st biological replication
│       ├── Combine fig 2nd trial-1.png  ## combine fig from 2nd biological replication
│       ├── Combine fig 3rd trial-1.png  ## combine fig from 3rd biological replication
│       ├── PVY on IL_1st trial_NB-1.png 
│       ├── PVY on IL_2nd trial_NB-1.png
│       ├── PVY on IL_3rd trial_NB-1.png
│       ├── PVY on SL_1st trial_NB-1.png
│       ├── PVY on SL_2ndtrial_NB-1.png
│       └── PVY on SL_3rdtrial_NB-1.png
├── README.md 
└── Reproducibility-project.Rproj
```
