# Systemic underestimation of aquatic neonicotinoid-associated ecological risks

![Suspect Screening](https://img.shields.io/badge/Suspect%20Screening-Neonicotinoids-2E8B57)
![Python](https://img.shields.io/badge/Python-3.14-3776AB?logo=python&logoColor=white)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-00A6D6)
![Monte Carlo](https://img.shields.io/badge/Monte%20Carlo-Simulation-6C3483)
![License](https://img.shields.io/badge/License-MIT-yellow)

This repository provides the data, code, and reproducible computational workflows supporting the suspect screening, structural characterization, and ecological risk assessment of parent neonicotinoids (p-NEOs), transformation products (t-NEOs), and structural analogues (a-NEOs).

<img width="850" height="265" alt="Figure1_GitHub" src="https://github.com/user-attachments/assets/841e096e-c297-4865-9d37-e7e15043c555" />

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



### 2. Structural Similarity Screening for a-NEOs ([2.Molecular_similarity.py](https://github.com/Eamon-Yang/Neonicotinoids/blob/main/2.Molecular_similarity.py))
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
**Code:** Please see the accompanying Python script `Molecular_similarity.py`.



### 3. Substructure-Based Screening of NEO-related compounds ([3.Substructure_extraction.py](https://github.com/Eamon-Yang/Neonicotinoids/blob/main/3.Substructure_extraction.py))
NEO-related compounds were systematically extracted from existing chemical databases through **substructure-based screening** using the **RDKit** cheminformatics toolkit in Python.
The screening was based on a set of predefined structural patterns associated with NEO-related compounds.
- **Structural patterns**
  - Six predefined NEO-related structural patterns were encoded as SMILES and converted into RDKit molecular query objects.
- **Substructure matching**
  - SMILES structures of compounds in the input database were parsed and validated using RDKit.
  - Each valid molecular structure was screened against the predefined structural patterns using RDKit substructure matching (`HasSubstructMatch`).
  - A compound was extracted if it contained at least one of the predefined structural patterns.
- **Input**
  - CAS Registry Number
  - SMILES
- **Output**
  - CAS Registry Number
  - SMILES
Invalid or empty SMILES entries were excluded prior to substructure matching.
**Code:** Please see the accompanying Python script `Substructure_extraction.py`.



### 4. Monte Carlo Simulation of Global NEO Risk Increments ([4.Monte_Carlo_simulation.py](https://github.com/Eamon-Yang/Neonicotinoids/blob/main/4.Monte_Carlo_simulation.py))
Monte Carlo simulation was used to estimate the potential increase in global ecological risks associated with transformation products (**t-NEOs**) and structural analogues (**a-NEOs**) relative to their corresponding parent neonicotinoids (**p-NEOs**).
The simulation was implemented in Python using observed risk increment factors derived from measured p-NEOs, t-NEOs, and a-NEOs, together with global parent-NEO risk quotients (RQs).
### Simulation Framework
- **Risk increment factors**
  - Observed risk increment factors greater than 1 (**M > 1**) were retained for distribution fitting.
  - The retained increment factors were log-transformed and fitted to log-normal distributions.
  - Distribution parameters were estimated separately for transformation products of individual p-NEOs and collectively for structural analogues.
- **Monte Carlo simulation**
  - Number of iterations: **10,000**
  - Random seed: **42**
  - Risk increment factors were randomly sampled from the fitted log-normal distributions.
  - Simulated increment factors were constrained to values ≥ 1.
  - Sampled factors were applied to parent-NEO RQs at global monitoring sites.
- **Transformation products**
  - Compound-specific simulations were conducted for parent NEOs
- **Structural analogues**
  - A collective analogue risk increment distribution was used to estimate the additional risk associated with a-NEOs relative to the summed RQs of p-NEOs.
- **Combined risk scenarios**
  - `p-NEOs`: parent NEOs only
  - `p+t-NEOs`: parent NEOs adjusted for transformation-product risk increments
  - `p+a-NEOs`: parent NEOs adjusted for structural-analogue risk increments
  - `p+t+a-NEOs`: transformation-product-adjusted RQs further adjusted using the structural-analogue increment factor
- **Simulation interval**
  - Risk increment ratios are summarized using the median and **95% simulation interval**, defined by the **2.5th and 97.5th percentiles** of the simulated ratio distribution.
- **Input**
  - The simulation requires two tab-delimited input files.
- **`RQ_Global.txt`**
  - Global parent-NEO risk quotients compiled from the literature-based global dataset.
  - The file also contains geographic information such as `Continent` for descriptive summaries.
- **`Magnification.txt`**
  - Observed risk increment factors derived from measured transformation products and structural analogues.
- **Output**
  - The script generates statistical summaries and cumulative probability plots.
**Code:** Please see the accompanying Python script `Monte_Carlo_simulation`.


### Repository Contents
- `1. Suspect screening database.xlsx`: The combined industrial chemical database.
- `2.Molecular_similarity.py`: Calculates molecular similarity between a reference NEO and compounds in an existing chemical database using Morgan fingerprints and the Tanimoto coefficient.
- `3.Substructure_extraction.py`: Screens existing chemical databases for compounds containing predefined NEO-related structural patterns using RDKit substructure matching.
- `4.Monte_Carlo_simulation.py`: Performs Monte Carlo simulations to estimate global ecological risk increments associated with t-NEOs and a-NEOs, including p+t-NEO, p+a-NEO, and p+t+a-NEO risk scenarios.


### Requirements
The scripts were implemented in Python and require the following packages:
- `Python`
- `NumPy`
- `pandas`
- `SciPy`
- `Matplotlib`
- `RDKit`


### References

If you use the data, code, or computational workflows provided in this repository, please cite:

> Yang, Z. et al. **Systemic underestimation of aquatic neonicotinoid-associated ecological risks.** *xxx* **202x**, *xx* (xx), xxx–xxx. DOI: 

> Rogers, D.; Hahn, M. **Extended-Connectivity Fingerprints.** *J. Chem. Inf. Model.* **2010**, *50* (5), 742–754. DOI: [10.1021/ci100050t](https://doi.org/10.1021/ci100050t)

> Bajusz, D. et al. **Why Is Tanimoto Index an Appropriate Choice for Fingerprint-Based Similarity Calculations?** *J. Cheminform.* **2015**, *7*, 20. DOI: [10.1186/s13321-015-0069-3](https://doi.org/10.1186/s13321-015-0069-3)

> Su, W. et al. **Identification and Prioritization of Emerging Organophosphorus Compounds Beyond Organophosphate Esters in Chinese Estuarine Waters.** *Environ. Sci. Technol.* **2026**. DOI: [10.1021/acs.est.4c09869](https://doi.org/10.1021/acs.est.4c09869)

> Li, P. et al. **Occurrence and Temporal Trends of Benzotriazole UV Stabilizers in Mollusks (2010–2018) from the Chinese Bohai Sea Revealed by Target, Suspect, and Nontarget Screening Analysis.** *Environ. Sci. Technol.* **2022**, *56* (23), 16759–16767. DOI: [10.1021/acs.est.2c04143](https://doi.org/10.1021/acs.est.2c04143)
