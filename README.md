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
- [Data](https://github.com/BibechanaPaudel/MS-thesis/tree/main/Data): Contains .csv files for each of the three biological replicates, and combined dataset of all three biological replicates that are divided based on objective.

- [At_Virus](https://github.com/BibechanaPaudel/MS-thesis/tree/main/At_Virus): Contains the R script, Rmd, HTML file, and a PDF with figures for all *A. thaliana* virus-accumulation experiments.

- [At_PGPR_Virus](https://github.com/BibechanaPaudel/MS-thesis/tree/main/At_PGPR_Virus): Contains the R script, Rmd, HTML file, and a PDF with figures for all *A. thaliana* virus-accumulation experiments after PGPR treatment.

- [At_Gene](https://github.com/BibechanaPaudel/MS-thesis/tree/main/At_Gene): Contains the R script, Rmd, HTML file, and a PDF with figures for all *A. thaliana* gene expression analysis experiments.

- [Nb_Virus](https://github.com/BibechanaPaudel/MS-thesis/tree/main/Nb_Virus): Contains the R script, Rmd, HTML file, and a PDF with figures for all *N. benthamiana* virus-accumulation experiments.

- [Nb_PGPR_Virus](https://github.com/BibechanaPaudel/MS-thesis/tree/main/Nb_PGPR_Virus): Contains the R script, Rmd, HTML file, and a PDF with figures for all *N. benthamiana* virus-accumulation experiments after PGPR treatment.

- [Nb_Gene](https://github.com/BibechanaPaudel/MS-thesis/tree/main/Nb_Gene): Contains the R script, Rmd, HTML file, and a PDF with figures for all *N. benthamiana* gene expression analysis experiments. 



### Script workflow

#### **1. Data Import and Preparation**

- The `Data/` folder contains three biological replicate (`.csv`) files, organized as Obj1, Obj2, and Obj3 for each objective. A separate combined file is also included, where each replicate represents the technical replicates within the three biological replicates.
- All data are processed in a different RMarkdown (`.Rmd`) script based on the host plant and experiment objective.
- Each file contains Cq values from RT-qPCR for virus quantification or gene expression analysis, along with metadata (Treatment, Dpi, Technical Replicate). The dataset includes a column Gene, which identifies each gene being tested.
- Viral load is calculated from Cq values using the regression equation **Y= −3.93X+49.153** for PVY and **Y = –3.353X + 37.416** for CMV, where Y= Cq value and X= copy number of virus from [Feng, J.L et al., (2006)](https://academic.oup.com/abbs/article/38/10/669/217), then log-transformed for normality.
- Gene expression is calculated using 2^(-ΔΔCq)[Taylor et al., (2019)](https://pubmed.ncbi.nlm.nih.gov/30654913/)

#### **2. Data Grouping and Statistical Analysis**

- Data are grouped by Treatment, Days post-inoculation (Dpi), and Technical Replicates.
- For virus-accumulation after PGPR treatment, log-transformed viral load is analyzed separately at each Dpi:
 - For each Dpi, a one-way ANOVA is fit using aov(logViralLoad ~ Treatment).
 - Fisher’s LSD test is performed with agricolae::LSD.test (α = 0.05, p.adj = "none") to compare Treatment means.
 - Significance letters (e.g., a, b, c) are extracted from lsd$groups and used for visualization.
- For gene expression analysis, genes with |Log2(FC)∣>1 are considered biologically significant.


### Citation
If you use this code for your research, please cite it using the Zenodo DOI provided here:

[![DOI](https://zenodo.org/badge/966418246.svg)](https://doi.org/10.5281/zenodo.15258548)

### Disclaimer
The code provided herein has not been peer-reviewed and may contain errors. Users are encouraged to test the code thoroughly and verify its accuracy for their specific applications. The author is not responsible for any errors or inaccuracies in the results generated by this script.

### File Tree
```
├── At_Gene
│   ├── At_PGPR_CMV_Gene.html
│   ├── At_PGPR_CMV_Gene.Rmd
│   ├── At_PGPR_Combine.html
│   ├── At_PGPR_Combine.pdf
│   └── At_PGPR_Combine.Rmd
├── At_PGPR_Virus
│   ├── At_PGPR_CMV.html
│   ├── At_PGPR_CMV.pdf
│   ├── At_PGPR_CMV.Rmd
│   ├── At_PGPR_CMV_Combine.pdf
│   ├── At_PGPR_PVY.html
│   ├── At_PGPR_PVY.pdf
│   ├── At_PGPR_PVY.Rmd
│   └── At_PGPR_PVY_files
│       └── figure-html
├── At_Virus
│   ├── At_CMV.html
│   ├── At_CMV.Rmd
│   ├── At_CMV_CP.pdf
│   ├── At_CMV_RNA3.pdf
│   ├── At_PVY.html
│   ├── At_PVY.Rmd
│   ├── At_PVY_PVY1.pdf
│   └── At_PVY_PVY3.pdf
├── CMV_RNA3_combine_paper.pdf
├── Data
│   ├── Obj1
│   │   ├── At_CMV_1st_CP_1ng.csv
│   │   ├── At_CMV_1st_RNA3_1ng.csv
│   │   ├── At_CMV_2nd_CP_1ng.csv
│   │   ├── At_CMV_2nd_RNA3_1ng.csv
│   │   ├── At_PVY_1st_PVY1_10ng.csv
│   │   ├── At_PVY_1st_PVY3_10ng.csv
│   │   ├── At_PVY_2nd_PVY1_10ng.csv
│   │   ├── At_PVY_2nd_PVY3_10ng.csv
│   │   ├── Nb_CMV_1st_CP_10ng.csv
│   │   ├── Nb_CMV_1st_CP_1ng.csv
│   │   ├── Nb_CMV_1st_RNA3_10ng.csv
│   │   ├── Nb_CMV_1st_RNA3_1ng.csv
│   │   ├── Nb_CMV_2nd_CP_1ng.csv
│   │   ├── Nb_CMV_2nd_RNA3_1ng.csv
│   │   ├── Nb_PVY_1st_PVY1_10ng.csv
│   │   ├── Nb_PVY_1st_PVY1_1ng.csv
│   │   ├── Nb_PVY_1st_PVY1_1ng.xlsx
│   │   ├── Nb_PVY_1st_PVY3_10ng.csv
│   │   ├── Nb_PVY_1st_PVY3_1ng.csv
│   │   ├── Nb_PVY_2nd_PVY1_10ng.csv
│   │   └── Nb_PVY_2nd_PVY3_10ng.csv
│   ├── Obj2
│   │   ├── At_PGPR_1st_trial.csv
│   │   ├── At_PGPR_2nd_trial.csv
│   │   ├── At_PGPR_3rd_trial.csv
│   │   ├── At_PGPR_CMV_Gene_1st_trial.csv
│   │   ├── At_PGPR_CMV_Gene_2nd_trial.csv
│   │   ├── At_PGPR_CMV_Gene_3rd_trial.csv
│   │   ├── At_PGPR_CMV_Gene_Combine.csv
│   │   ├── At_PGPR_Combine.csv
│   │   ├── Nb_PGPR_1st_trial.csv
│   │   ├── Nb_PGPR_2nd_trial.csv
│   │   ├── Nb_PGPR_3rd_trial.csv
│   │   ├── Nb_PGPR_CMV_Gene_1st_trial.csv
│   │   ├── Nb_PGPR_CMV_Gene_2nd_trial.csv
│   │   ├── Nb_PGPR_CMV_Gene_3rd_trial.csv
│   │   ├── Nb_PGPR_CMV_Gene_Combine.csv
│   │   └── Nb_PGPR_Combine.csv
│   └── Obj3
│       ├── At_PGPR_CMV_1st_trial_New.csv
│       ├── At_PGPR_CMV_2nd_trial_New.csv
│       ├── At_PGPR_CMV_3rd_trial_New.csv
│       ├── At_PGPR_CMV_Combine.csv
│       ├── At_PGPR_PVY_1st_trial_New.csv
│       ├── At_PGPR_PVY_2nd_trial_New.csv
│       ├── At_PGPR_PVY_3rd_trial_New.csv
│       ├── At_PGPR_PVY_Combine.csv
│       ├── Nb_PGPR_CMV_1st_trial_New.csv
│       ├── Nb_PGPR_CMV_2nd_trial_New.csv
│       ├── Nb_PGPR_CMV_3rd_trial_New.csv
│       ├── Nb_PGPR_CMV_Combine.csv
│       ├── Nb_PGPR_PVY_1st_trial_New.csv
│       ├── Nb_PGPR_PVY_2nd_trial_New.csv
│       ├── Nb_PGPR_PVY_3rd_trial_New.csv
│       └── NB_PGPR_PVY_Combine.csv
├── MS-thesis.Rproj
├── Nb_Gene
│   ├── NB_PGPR_CMV_Combine.pdf
│   ├── Nb_PGPR_CMV_Gene.html
│   ├── Nb_PGPR_CMV_Gene.Rmd
│   ├── Nb_PGPR_Combine.html
│   ├── Nb_PGPR_Combine.md
│   ├── NB_PGPR_Combine.pdf
│   ├── Nb_PGPR_Combine.Rmd
│   └── Nb_PGPR_Combine_files
│       └── figure-gfm
│           ├── JA_Nb-1.png
│           ├── SA_NB-1.png
│           └── unnamed-chunk-8-1.png
├── Nb_PGPR_Virus
│   ├── NB_PGPR_CMV.html
│   ├── NB_PGPR_CMV.pdf
│   ├── NB_PGPR_CMV.Rmd
│   ├── Nb_PGPR_PVY.html
│   ├── NB_PGPR_PVY.pdf
│   ├── Nb_PGPR_PVY.Rmd
│   └── Nb_PGPR_PVY_files
│       └── figure-html
├── Nb_PVY.Rmd
├── Nb_Virus
│   ├── Nb_CMV.html
│   ├── Nb_CMV.Rmd
│   ├── NB_CMV_CP.pdf
│   ├── NB_CMV_RNA3.pdf
│   ├── NB_CMV_RNA3.png
│   ├── Nb_PVY.html
│   ├── Nb_PVY.Rmd
│   ├── NB_PVY_PVY1.pdf
│   └── Nb_PVY_PVY3.pdf
├── PVY_PVY3_Combine_Paper.pdf
├── README.html
└── README.md
```
