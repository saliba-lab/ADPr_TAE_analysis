# MBE analysis (ADPr-TAE project)
This repository contains scripts for analysis and visualization of mutational base editing (MBE) starting from CRISPResso outputs. They were used in the following paper:

[**Landscape and age dynamics of immune cells in the Egyptian rousette bat (2022)**](https://www.helmholtz-hiri.de/en/research/organisation/teams/team/rna-synthetic-biology/)

[Virginia Friedrichs, Christophe Toussaint, Alexander Schäfer, Melanie Rissmann, Oliver Dietrich, Thomas C. Mettenleiter, Gang Pei, Anne Balkema-Buschmann, Antoine-Emmanuel Saliba, Anca Dorhoi](https://www.helmholtz-hiri.de/en/research/organisation/teams/team/rna-synthetic-biology/)

[*https://doi.org/10.1016/j.celrep.2022.111305*](https://www.helmholtz-hiri.de/en/research/organisation/teams/team/rna-synthetic-biology/)

## Data Accessibility
Example data to run the scripts are available in the **data** folder of this repository.

The actual raw sequencing data used in the paper are available at [SRA](https://www.helmholtz-hiri.de/en/research/organisation/teams/team/rna-synthetic-biology/) and should be processed with CRISPResso beforehand.

Generally speaking, these scripts can be reused (after modifications) with any result table obtained from CRISPResso.

## Repository Structure
**data:** Directory containing example data.

**analysis:** Directory containing the R and python scripts.

**outputs:** This directory is created when running the scripts. It will contain the processed data and different tables.

## Code Execution
#### 1- Download repository
Option 1: Download manually the repository as a ZIP archive and extract it locally on your computer

Option 2: Clone the repository
```shell
git clone https://github.com/saliba-lab/MBE_analysis.git
cd MBE_analysis/analysis
```


#### 2- Install R dependencies 
See Dependencies section.


#### 3- Run sequentially the scripts in the **analysis** directory 
Make sure to set the **analysis** directory as the working directory when running the scripts.

[Must add a short description of each script]

## Dependencies
List of R packages necessary to run the scripts.

- R                 4.0.3
- dplyr             1.0.10
- openxlsx          4.2.3
- [Must be completed]

