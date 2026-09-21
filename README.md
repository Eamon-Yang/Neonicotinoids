# Systemic underestimation of aquatic neonicotinoid-associated ecological risks
<p align="left">
<img src="https://img.shields.io/badge/Python-3776AB.svg?style&logo=Python&logoColor=white" alt="Python" />
<img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT License" />
</p>
This repository provides a rational workflow for the identification and prioritization of pollutants in environmental samples by the combination of nontarget screening and MCDA prioritization approach.
<br>
<br>
<p align="center">
  <img src="Graphical_Abstract.png" alt="Graphical_Abstract" width="600"/>
</p>

---

## Overview

### 1. **Suspect screening database** ([1. Suspect screening database.xlsx](https://github.com/Eamon-Yang/Neonicotinoids/blob/main/1.Suspect_screening_database.xlsx))

The suspect screening database contains **1,349 compounds**, including parent neonicotinoids (**p-NEOs**), transformation products (**t-NEOs**), and structural analogues (**a-NEOs**).

The database was compiled and expanded using the following sources and approaches:

- **Industrial chemical inventories**
  - U.S. Toxic Substances Control Act (**TSCA**) Inventory
  - Canadian Domestic Substances List (**DSL**)
  - European Union **REACH** database
  - Inventory of Existing Chemical Substances in China (**IECSC**)

- **Chemical databases and suspect lists**
  - U.S. EPA **CompTox Chemicals Dashboard**
  - U.S. EPA **ToxCast** chemical library
  - **NORMAN Suspect List Exchange**

- **Literature-based screening**
  - Web of Science searches using combinations of *"neonicotinoid insecticides"* or *"neonicotinoids"* with environmental and pesticide-related terms.

- **Structure-based screening**
  - PubChem searches using characteristic NEO structural fragments, including **2-chloro-5-methylpyridine** and **2-chloro-5-methylthiazole**.

- **Transformation product prediction**
  - **BioTransformer 3.0**
  - **Chemical Transformation Simulator (CTS) v1.3.2.2**
  - **enviPath**

All data sources were accessed on **October 1, 2024**.

**Database:** Please see the separate Excel file for the complete list of **1,349 suspect compounds**.



### 2. Structural Similarity Analysis for a-NEOs ([2.Molecular_similarity.py](https://github.com/Eamon-Yang/Neonicotinoids/blob/main/2.Molecular_similarity.py))

Potential structural analogues of neonicotinoids (**a-NEOs**) were identified through structural similarity analysis followed by manual screening and structural verification.

The analysis was performed in **Python** using the **RDKit** cheminformatics toolkit.

- **Reference compound**
  - A representative parent neonicotinoid (**p-NEO**) was used as the reference structure for similarity calculation.

- **Molecular representation**
  - Morgan fingerprints
  - Radius = **2**
  - Bit vector length = **2048**
  - Pharmacophoric features enabled (`useFeatures=True`)

- **Similarity calculation**
  - Structural similarity between the reference p-NEO and each compound in the database was calculated using the **Tanimoto coefficient**.
  - Similarity scores were calculated for all compounds with valid SMILES and exported for subsequent screening.

- **Candidate selection**
  - Compounds with a Tanimoto similarity ≥ **0.5** were manually selected as candidate structural analogues.

- **Structural verification**
  - Candidate compounds were further manually curated based on characteristic structural features of NEOs.
  - Particular attention was given to the presence of the **2-chloro-5-methylpyridine moiety** (`ClC1=NC=C(C)C=C1`), a characteristic structural moiety occurring in many NEO-related compounds.

- **Input**
  - CAS Registry Number
  - SMILES

- **Output**
  - CAS Registry Number
  - SMILES
  - Tanimoto similarity score

Invalid or empty SMILES entries were excluded from the similarity calculation.

**Code:** Please see the accompanying Python script Molecular_similarity.py.






### 2. **t-NEOs and a-NEOs Extraction** ([2_Extract_OPCs.py](https://github.com/WestonSu/Organophosphorus/blob/main/2_Extract_OPCs.py))
The `2_Extract_OPCs.py` script extracts OPCs by:
- Filtering out compounds containing counterions (e.g., Na+/K+/Cl−/Br−) and metal/metalloid-containing compounds.
- Standardizing the resulting list into an "MS-ready" format for HRMS analysis.

### 5. **t-NEOs and a-NEOs similarity calculation** ([5_Tanimoto.py](https://github.com/WestonSu/Organophosphorus/blob/main/5_Tanimoto.py))
This script calculates the structural similarity between compounds using the **Tanimoto coefficient**, which is commonly used in cheminformatics for molecular similarity assessments.


## Repository Contents
- `1_Industrial_chemicals.csv`: The combined industrial chemical database.
- `2_Extract_NEO-related compounds.py`: Script to extract OPCs and generate an MS-ready suspect list.
- `5_Tanimoto similarity.py`: Python script for calculating structural similarity between chemicals.

### Python Packages:
- `pandas`
- `numpy`
- `scikit-learn`

## Citations 
We ask users to directly cite the following paper:

This project also builds on a number of other projects, algorithms and ideas. Please consider citing the following full list of papers when relevant: 

1. Su, W.; Li, P.; Zhong, L.; Liang, W.; Li, T.; Liu, J.; Ruan, T.; Jiang, G. Occurrence and distribution of antibacterial quaternary ammonium compounds in Chinese estuaries revealed by machine learning-assisted mass spectrometric analysis. Environ. Sci. Technol. 2024. 58 (26), 11707−11717. DOI: 10.1021/acs.est.4c02380
2. Helmus, R.; ter Laak, T. L.; van Wezel, A. P.; de Voogt, P.; Schymanski, E. L. patRoon: open source software platform for environmental mass spectrometry based non-target screening. J. Cheminform. 2021, 13 (1), 1. DOI: 10.1186/s13321-020-00477-w
3. Ruttkies, C.; Schymanski, E. L.; Wolf, S.; Hollender, J.; Neumann, S. MetFrag relaunched: incorporating strategies beyond in silico fragmentation. J. Cheminform 2016, 8 (1), 3. DOI: 10.1186/S13321-016-0115-9
4. Dührkop, K.; Fleischauer, M.; Ludwig, M.; Aksenov, A. A.; Melnik, A. V.; Meusel, M.; Dorrestein, P. C.; Rousu, J.; Böcker, S. SIRIUS 4: a rapid tool for turning tandem mass spectra into metabolite structure information. Nat. Methods 2019, 16 (4), 299−302. DOI: 10.1038/s41592-019-0344-8
5. Su, W. et al. Identification and Prioritization of Emerging Organophosphorus Compounds Beyond Organophosphate Esters in Chinese Estuarine Waters. Environ. Sci. Technol. 2025, 59, 8, 4080–4091. DOI: 10.1021/acs.est.4c09869

