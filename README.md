<h1 align="center"><em>Tribolium castaneum</em> transcriptomics pipeline</h1>
<h3 align="center">ECOEVODYN GROUP VERSION (PRIVATE REPOSITORY)</h3>
<p align="center">
<img src="https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation/assets/25666459/0536e2b9-67b2-4e8a-9227-a8ca3354eb63">
<br/>
</p>

This repository contains the transcriptomics pipeline I developed during my postdoc in the [Eco-Evolutionary Dynamics group](https://www.helsinki.fi/en/researchgroups/eco-evolutionary-dynamics), from Sept. 2021 to Dec. 2022. It deals with various facets of the transcriptomic data obtained from the _Tribolium castaneum_ laboratory adaptive experiment published in [Koch & Guillaume (2020a)](https://doi.org/10.1371/journal.pgen.1008768), [Koch & Guillaume (2020b)](https://doi.org/10.1111/mec.15607) and [Koch et al. (2020)](https://doi.org/10.1111/evo.14119).

The pipeline is described here for publication and educational purposes.

Do not hesitate to <a href="mailto:charles DOT rocabert AT hhu DOT de">contact me</a> for any remarks/questions.

# Table of contents
- [Overview](#overview)
- [Authors](#authors)
- [Copyright](#copyright)
- [Publications](#publications)
- [Dependencies](#dependencies)
- [Description of the pipeline](#pipeline_description)
  - [Scripts](#scripts)
    - [Puhti scripts](#scripts_0)
    - [Reorganizing BAM files](#scripts_1)
    - [Variants detection](#scripts_2)
    - [Selecting populations](#scripts_3)
    - [Genotype imputation](#scripts_4)
    - [Allelic frequency changes (AFCs)](#scripts_5)
    - [Preparing read counts](#scripts_6)
    - [Preparing phenotypes](#scripts_7)
    - [eQTLs](#scripts_8)
    - [WGCNA modules](#scripts_9)
    - [Allele specific expression (ASE)](#scripts_10)
    - [LD map with Lep-MAP3](#scripts_11)
    - [Haplotype blocks](#scripts_12)
    - [FST clusters](#scripts_13)
  - [Analyses](#analyses)
  - [Data](#data)
- [License](#license)

# Overview

![final_pipeline](https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation/assets/25666459/61c9eded-050b-42f4-b29e-a79c62fb6438)

# Authors <a name="authors"></a>

Charles Rocabert, Eva L. Koch, Frédéric Guillaume.

# Copyright <a name="copyright"></a>

Copyright © 2021-2023 Charles Rocabert, Eva L. Koch, Frédéric Guillaume. All rights reserved.

# Publications <a name="publications"></a>

• Eva L. Koch, Charles Rocabert, Frédéric Guillaume, in prep.

• Charles Rocabert, Eval L. Koch, Frédéric Guillaume. Deciphering the genetic architecture of polygenic adaptation in Tribolium castaneum. <em>ESEB 2022: Congress of the European Society for Evolutionary Biology</em> (Aug. 2022, Prague, Czech Republic)

# Dependencies <a name="dependencies"></a>

### • Software:
- PLINK 2.0,
- GATK-4.2.3.0,
- samtools 1.14
- picard-toolkit
- bcftools-1.14
- snpEff-5.1
- Beagle-5.4
- Lep-MAP3
- subread-2.0.3
- GEMMA-0.98.5

### • Programming languages:
- Python3+
- R
- shell

### • R-packages
- tidyverse
- cowplot
- edgeR
- limma
- ggfortify
- tibble
- WGCNA

### • Python libraries
- swiftclient
- paramiko

# Description of the pipeline <a name="pipeline_description"></a>

      ├── scripts
      ├── analyses
      ├── data
      └── README.md

The pipeline is splitted in three categories:
- `scripts`: This folder contains all the specific omics tasks to build the data that will be used for further, higher-level analyses. This includes _e.g._ variants detection, eQTLs analysis, ASE, etc.
- `analyses`: This folder contains higher-level analyses. This include _e.g._ the detection of signals of selection by interesecting gene expression and allele frequency datasets.
- `data`: This folder contains all the data produced by the pipeline. This folder is **not included** in this repository. A download link will be provided in the future.

The pipeline is described in details below.

## Scripts <a name="scripts"></a>

      └── scripts
           ├── 0_Puhti_scripts
           ├── 1_BAM_files_reorganization
           ├── 2_variant_call
           ├── 3_select_population
           ├── 4_genotype_imputation
           ├── 5_AFCs
           ├── 6_read_counts
           ├── 7_phenotypes
           ├── 8_eQTLs
           ├── 9_WGCNA_modules
           ├── 10_ASE
           ├── 11_LD_map
           ├── 12_haplotype_blocks
           └── 13_Fst_clusters

Omics tasks are numbered and separated into folders for the purpose of clarity.
For each task, the scripts are globally numbered in the order of their execution, and are split between local (`local` folder) and HPC (`hpc` folder scripts).

Sometimes, a shell script is also available to run all **local scripts** in the right order (see below). For HPC scripts, the user must update and deploy Puhti scripts (`./scripts/0_puhti_scripts`, see below).
    
### 📂 Puhti scripts <a name="scripts_0"></a>

      └── scripts
           └── 0_Puhti_scripts

This folder contains various scripts and wrappers to deploy Python or R tasks on Puhti (single jobs or array jobs), and connect to Allas.

⚠️ **User and project names must be updated before usage**.

```mermaid
flowchart LR
A("Single job task") --> B[["Run wrapper"]]
```

```mermaid
flowchart LR
A("Parallel jobs task") --> B[["Array wrapper"]]
```

### 📂 Reorganizing BAM files <a name="scripts_1"></a>

      └── scripts
           └── 1_BAM_files_reorganization
                └── local
                     ├── 1_CreateBamMap.py
                     ├── 2_SplitBamMap.R
                     ├── 3_CheckBamReadgroups.py
                     └── 4_RelocateBamFiles.py

The goal of this task is to parse the original source of BAM files (here a hard-drive) to rationalize their organization, make some adjustements and extract information.

Associated data folder(s): `./data/tribolium_bam`.

#### ⚙️ `1_CreateBamMap.py (local)`:
> This script parses the original BAM files folder and extract the information in the form of a table `./data/tribolium_bam/bam_map.csv`.

#### ⚙️ `2_SplitBamMap.R (local)`:
> This script splits the file `bam_map.csv` into many sample subsets, depending on the reference genome, the environment, the line or the generation. All files are stored in `./data/tribolium_bam` folder.
**Sample files will be used all along the pipeline**.

#### ⚙️ `3_CheckBamReadgroups.py (local)`:
> This script parses every BAM files to check the absence of the "read group" entry ("RG" label), which is mandatory for further analyses (will be added later).

#### ⚙️ `4_RelocateBamFiles.py (local)`:
> This script relocates BAM files from the original hard-drive for further analysis (:warning: needs debugging).
> **Ultimately, BAM files are transfered to Allas with an independent script**.

```mermaid
flowchart LR
subgraph local
direction LR
A[("Source<br/>BAMs")] --> B("1_CreateBamMap.py<br/><em>(local)</em>")
B --> C[("BAM<br/>map")]
C --> D("2_SplitBamMap.R<br/><em>(local)</em>")
D --> E[("Samples")]
E --> F("3_CheckBamReadgroups.py<br/><em>(local)</em>")
A[("Source<br/>BAMs")] --> F
E --> G("4_RelocateBamFiles.py<br/><em>(local)</em>")
A[("Source<br/>BAMs")] --> G
G --> H[("Ready for<br/>analysis<br/>BAMs")]
end
```

### 📂 Variants detection <a name="scripts_2"></a>

      └── scripts
           └── 2_variant_call
                └── hpc
                     ├── 1_BamPreProcessing.py
                     ├── 2_HaplotypeCaller.py
                     └── 3_JointCaller.py
                └── local
                     └── 4_SelectFilterAnnotateVariants.py
                └── variant_detection_pipeline.sh

The goal of this task is to perform a variant call on the entire transcriptomic data (version Tcas3.30 of the genome). It uses various software and ultimately produces a filtered VCF file containing raw SNPs.
This raw SNPs VCF file is stored in `./data/tribolium_vcf/Tribolium_castaneum_ALL_Tcas3.30.vcf.gz`.

Associated data folder(s): `./data/tribolium_vcf`.

#### ⚙️ `1_BamPreProcessing.py` (HPC array wrapper):
> This script pre-processes BAM files by (in this order):
> - Importing the unedited BAM file from Allas,
> - Importing the reference genome and generate indices,
> - Adding the read group,
> - Uncompressing the BAM file to SAM,
> - Editing the SAM file to recalibrate MAPQ values (255 --> 60),
> - Compressing edited SAM file to BAM,
> - Copying a version of the BAM file to mark duplicates,
> - Generating BAI index file,
> - Handling splicing events,
> - Exporting edited BAM files to Allas.

#### ⚙️ `2_HaplotypeCaller.py` (HPC array wrapper):
> This script uses marked duplicates BAM files to run the complete pipeline for per-sample variant call by (in this order):
> - Importing the BAM file from Allas
> - Importing the reference genome and generate indexes,
> - Generating BAI index file,
> - Running GATK HaplotypeCaller. This task can take several hours,
> - Exporting GVCF files to Allas.

#### ⚙️ `3_JointCaller.py` (HPC run wrapper):
> This script runs the complete pipeline for the join call by (in this order):
> - Generating GenomicsDB database by:
>   - Importing the list of samples,
>   - Importing all GVCFs from Allas,
>   - Importing the reference genome from Allas and compute index files,
>   - Generating the sample map,
>   - Generating the interval list,
>   - Running GATK GenomicsDBImport,
>   - Exporting the database to the scratch.
> - Performing the joint-call cohort by:
>   - Importing the consolidated database from the scratch if needed,
>   - Running GATK GenotypeGVCFs,
>   - Exporting the joint-call file (VCF) to the scratch.

#### ⚙️ `4_SelectFilterAnnotateVariants.py` (local):
> This script selects bi-allelic SNP variants from the original VCF file by (in this order):
> 
> - Selecting bi-allelic SNP variants,
> - Tagging `DP=0` genotypes as missing (`./.`),
> - Filtering out low quality SNPs,
> - Annotating SNPs,
> - Adding SNP unique identifiers.

**➡️ Local script to run with the shell script `variant_detection_pipeline.sh`**.

```mermaid
flowchart TB
subgraph HPC
direction LR
A[("Unedited<br/>BAMs")] --> B("1_BamPreProcessing.py<br/><em>(HPC array)</em>")
B --> C[("Edited<br/>BAMs")]
C --> D("2_HaplotypeCaller.py<br/><em>(HPC array)</em>")
D --> E[("GVCFs")]
E --> F("3_JointCaller.py<br/><em>(HPC)</em>")
F --> G[("VCF")]
end
subgraph local
direction LR
H("variant_detection_pipeline.sh<br/><em>(local)</em>") --> I[("bi-allelic SNPs<br/>VCF")]
end
HPC --> |"Download<br/>VCF"| local
```

### 📂 Selecting populations <a name="scripts_3"></a>

      └── scripts
           └── 3_select_population
                └── local
                     ├── 1_SelectPopulation.py
                     ├── 2_FilterGenotypes.py
                     └── 3_VariantsToTable.py
                ├── VCF_ALL_Beagle_pipeline.sh
                ├── VCF_CT_G1_eQTL_imputed_pipeline.sh
                ├── VCF_HD_G1_eQTL_imputed_pipeline.sh
                ├── VCF_imputed_genotypes_line_separation_pipeline.sh
                └── VCF_CT_HD_G1_LepMAP3_pipeline.sh

This task is a general function which separate population-level files (VCF or read counts) into sub-population files on user request.
**This function is used in further tasks involving sub-population analyses**.

Associated data folder(s): `./data/tribolium_snp`, `./data/tribolium_counts`.

#### ⚙️ `1_SelectPopulation.py` (local):
> This script selects a sub-population of samples based on a sample list (see section [3.1.1](#scripts_1")) by selecting the subset of samples in a VCF file or a read counts file.

#### ⚙️ `2_FilterGenotypes.py` (local):
> This script filters SNPs based on genotype coverage by (in this order):
> - Filtering genotypes based on F_MISSING,
> - Removing non-variant SNPs,
> - Removing LOWQUAL SNPs.

#### ⚙️ `3_VariantsToTable.py` (local):
> This script extracts VCF statistics into a text table.

**➡️ 5 pipelines are available for the selection of VCF sub-populations:**
- `VCF_ALL_Beagle_pipeline.sh`: Select ALL genotypes for Beagle imputation pipeline,
- `VCF_CT_G1_eQTL_imputed_pipeline.sh`: Select CT-G1 imputed genotypes for the eQTLs analysis,
- `VCF_HD_G1_eQTL_imputed_pipeline.sh`: Select HD-G1 imputed genotypes for the eQTLs analysis,
- `VCF_imputed_genotypes_line_separation_pipeline.sh`: Separate imputed genotypes into lines and generations for the AFC pipeline,
- `VCF_CT_HD_G1_LepMAP3_pipeline.sh`: Select CT/HD-G1 genotypes for the lep-MAP3 pipeline.

```mermaid
flowchart LR
subgraph local
direction LR
A[("VCF")] --> B("1_SelectPopulation.py<br/><em>(local)</em>")
B --> C("2_FilterGenotypes.py<br/><em>(local)</em>")
C --> D("3_VariantsToTable.py<br/><em>(local)</em>")
D --> E[("Sub-population<br/>VCF")]
end
```

**OR**

```mermaid
flowchart LR
subgraph local
direction LR
A[("Read<br/>counts")] --> B("1_SelectPopulation.py<br/><em>(local)</em>")
B --> C[("Sub-population<br/>read<br/>counts")]
end
```

### 📂 Genotype imputation <a name="scripts_4"></a>

      └── scripts
           └── 4_genotype_imputation
                └── local
                     ├── 1_ExtractNoMissingMarkers.py
                     ├── 2_ImputationTests.py
                     └── 3_ImputeGenotypes.py
                └── VCF_ALL_imputation_pipeline.sh

This task imputes missing genotypes on the global VCF file (including all samples) using Beagle software.
It is also possible to run a test to evaluate the quality of genotype imputations.

Associated data folder(s): `./data/tribolium_snp`.

#### ⚙️ `1_ExtractNoMissingMarkers.py` (local):
> This script extracts from the raw VCF file all the markers with a 100% call rate (_i.e_ without missing genotypes).
> The dataset is saved in binary format in `./data/tribolium_snp/imputation_tests`.
> 
#### ⚙️ `2_ImputationTests.py` (local):
> This script tests imputation capabilities of Beagle with toy datasets.
> Test results are saved in `./data/tribolium_snp/imputation_tests`.

#### ⚙️ `3_ImputeGenotypes.py` (local):
> Imputes genotypes of the VCF file with Beagle.

**➡️ Local script to run with the shell script `VCF_ALL_imputation_pipeline.sh`**.

```mermaid
flowchart LR
subgraph "local (benchmark)"
direction LR
B("1_ExtractNoMissingMarkers.py<br/><em>(local)</em>") --> C("2_ImputationTests.py<br/><em>(local)</em>")
C --> D[("Imputation<br/>success rate")]
end
subgraph "local (main pipeline)"
direction LR
E("2_ImputeGenotypes.py<br/><em>(local)</em>") --> F[("Imputed<br/>VCF")]
end
A[("VCF")] --> B
A --> E
```

### 📂 Allelic frequencies changes (AFCs) <a name="scripts_5"></a>

      └── scripts
           └── 5_AFCs
                └── local
                     ├── 1_MergeAFs.R
                     ├── 2_ComputeAFCs.R
                     ├── 3_SplitAFCs.R
                     └── 4_PlotAFCDistributions.R

This task computes AFCs per line and environment and splits the resulting dataset by line for further analyses.

Associated data folder(s): `./data/tribolium_afc`.

#### ⚙️ `1_MergeAFs.R` (local):
> This script merges and saves allelic frequencies.

#### ⚙️ `2_ComputeAFCs.R` (local):
> This script computes AFCs and saves the result.

#### ⚙️ `3_SplitAFCs.R` (local):
> This script Split AFCs by line and environment and saves the result.

#### ⚙️ `4_PlotAFCDistributions.R` (local):
> Optional script to display AFCs distibutions.

```mermaid
flowchart LR
subgraph "local"
direction LR
A[("SNPs<br/>table")] --> B("1_MergeAFs.R<br/><em>(local)</em>")
B --> C("2_ComputeAFCs.R<br/><em>(local)</em>")
C --> D("3_SplitAFCs.R<br/><em>(local)</em>")
D --> E[("AFCs<br/>per line<br/>and environment")]
end
```

### 📂 Preparing read counts <a name="scripts_6"></a>

      └── scripts
           └── 6_read_counts
                └── hpc
                     ├── 1_FeatureCounts.py
                     └── 2_MergeFeatureCounts.py
                └── local
                     ├── 3_TransformReadCounts.R
                     ├── 4_DetectLowExpressedTranscripts.R
                     └── 5_StandardizeReadCounts.R
                └── read_counts_preparation_pipeline.sh

This task handles BAM files (without marked duplicates) to perform a reads count. The final output is a text file containing read counts for every called samples.

Associated data folder(s): `./data/tribolium_counts`.

#### ⚙️ `1_FeatureCounts.py` (HPC array wrapper):
> This script calculates feature counts for every individual samples by (in this order):
> - Importing `subread` package and compile it,
> - Importing reference genome annotation,
> - Importing the list of samples,
> - For each BAM file, if the read counts file does not exist:
>   - Running `subread FeatureCounts`,
>   - Exporting resulting files to Allas.

#### ⚙️ `2_MergeFeatureCounts.py` (HPC run wrapper):
> This script merges feature counts from every individual samples by (in this order):
> - Importing the list of samples and all read counts,
> - Mergeing read counts in a single file,
> - Exporting the resulting file to Allas.

#### ⚙️ `3_TransformReadCounts.R` (local):
> This script transforms a gene expression dataset by filtering it, calculating the TMM normalization, and removing run, batch and line effects.

#### ⚙️ `4_DetectLowExpressedTranscripts.R` (local):
> This script detects low expressed transcripts and save the list of expressed ones.

#### ⚙️ `5_StandardizeReadCounts.R` (local):
> This script standardizes a gene expression dataset by quantile normalization.

**⚙➡️ Local scripts to run with the shell script `read_counts_preparation_pipeline.sh`**.

```mermaid
flowchart TB
subgraph HPC
direction LR
A[("BAMs")] --> B("1_FeatureCounts.py<br/><em>(HPC array)</em>")
B --> C[("Per-sample<br/>read counts")]
C --> D("2_MergeFeatureCounts.py<br/><em>(HPC array)</em>")
D --> E[("Global<br/>read counts")]
end
subgraph local
direction LR
F("3_TransformReadCounts.R<br/><em>(local)</em>") --> H[("Transformed<br/>read counts")]
G("4_DetectLowExpressedTranscripts.R<br/><em>(local)</em>") --> I[("Expressed<br/>reads")]
H --> J("5_StandardizeReadCounts.R<br/><em>(local)</em>")
I --> J
J --> K[("Standardized<br/>read counts")]
end
HPC --> |"Download<br/>read counts"| local
```

### 📂 Preparing phenotypes <a name="scripts_7"></a>

      └── scripts
           └── 7_phenotypes
                └── local
                     ├── 1_CopyExpressionFiles.sh
                     ├── 2_ComputePlasticityResponse.R
                     ├── 3_ComputePhenotypicNoise.R
                     └── 4_ComputeRelativeFitness.R
                └── phenotype_preparation_pipeline.sh

This task uses ready-to-use gene expression data to compute various phenotypes, including plasticity and phenotypic noise.
The task also computes relative fitnesses.
**Calculated phenotypes will be used in further analyses, _e.g._ eQTLs or WGCNA modules**.

Associated data folder(s): `./data/tribolium_phenotypes`.

#### ⚙️ `1_CopyExpressionFiles.sh` (local):
> This script copies ready-to-use expression files to the folder `./data/tribolium_phenotypes/`.

#### ⚙️ `2_ComputePlasticityResponse.R` (local):
> This script computes plasticity response between CT and HD.

#### ⚙️ `3_ComputePhenotypicNoise.R` (local):
> This script estimates phenotypic noise (V<sub>e</sub>) in HD.

#### ⚙️ `4_ComputeRelativeFitness.R` (local):
> This script computes relative fitnesses per line for CT and HD individuals (G1).

```mermaid
flowchart LR
subgraph "local"
direction LR
A[("read<br/>counts")] --> B("1_CopyExpressionFiles.sh<br/><em>(local)</em>")
B --> C("2_ComputePlasticityResponse.R<br/><em>(local)</em>")
B --> D("3_ComputePhenotypicNoise.R<br/><em>(local)</em>")
C --> E[("Plastic response<br/>phenotypes")]
D --> F[("V<sub>e</sub><br/>phenotypes")]
G[("Fitness<br/>measurements")] --> H("4_ComputeRelativeFitness.R<br/><em>(local)</em>")
H --> I[("Relative<br/>fitnesses")]
end
```

### 📂 eQTLs <a name="scripts_8"></a>

      └── scripts
           └── 8_eQTLs
                └── hpc
                     ├── CheckGemmaFiles.py
                     ├── DeleteGemmaFiles.py
                     ├── 3_Gemma.py
                     └── 4_CollectSignificantEQTLs.R
                └── local
                     ├── EditFam.R
                     ├── ExtractGenePos.py
                     ├── MergeEQTLsDatasets.R
                     ├── ComputeCorrelationToFitness.R
                     ├── 1_eQTLs_data_preparation_pipeline.sh
                     ├── 2_UploadGemmaFiles.py
                     ├── 5_DownloadGemmaFiles.py
                     └── 6_merge_eQTLs_pipeline.sh
                └── phenotype_preparation_pipeline.sh

This task runs GWAAs on various phenotypes using GEMMA software to detect significant eQTLs at genome-scale.

Associated data folder(s): `./data/tribolium_eqtl`.

#### ⚙️ `1_eQTLs_data_preparation_pipeline.sh` (local):
> This script generates input files necessary to the GWAA with GEMMA software:
> - `.bed` and `.bim` files using `plink2`,
> - `.fam` and `.pheno` files using the script `EditFam.R`.
> Files are generated at once for all the phenotypes (expression, plasticity, noise and fitness). 

#### ⚙️ `2_UploadGemmaFiles.py` (local):
> This script exports all the input files (`.bed`, `.bim`, `.fam` and `.pheno`) to the distant storage server Allas.

#### ⚙️ `CheckGemmaFiles.py` (HPC maintenance script):
> This maintenance script checks if GEMMA output files are missing after an array job.

#### ⚙️ `DeleteGemmaFiles.py` (HPC maintenance script):
> This maintenance script delete GEMMA output files to prevent duplicates before a new run.

#### ⚙️ `3_Gemma.py` (HPC array wrapper):
> This script runs GWAAs on a given number of phenotypes (the max number of jobs being limited on Puhti) by (in this order):
> - Importing GEMMA software and datasets from Allas,
> - Calculating the kinship matrix,
> - Running GEMMA on the right set of phenotypes,
> - Converting the output to RDS format,
> - Exporting resulting files to Allas.

#### ⚙️ `4_CollectSignificantEQTLs.R` (HPC independent script):
> This script collect significant eQTLs in a dedicated file saved in the scratch. Two steps are applied to filter significant eQTLs:
> Collect all significant eQTL associations given a p-value threshold.
> - Genomic correction is applied and the FDR is calculated,
> - eQTLs with a p-value < 0.05 are selected.

#### ⚙️ `5_DownloadGemmaFiles.py` (local):
> This script downloads significant eQTL files from Allas, based on a list of phenotypes.

#### ⚙️ `6_merge_eQTLs_pipeline.sh (local):
> This pipeline collect and merge significant eQTLs with additional information such as gene position (`ExtractGenePos.py` script) and gene and phenotype annotations (`MergeEQTLsDatasets.R` script).
> All phenotypes (expression, plasticity, noise and fitness) are treated at once.
> This script also calculates the correlation to fitness of every phenotypes (`ComputeCorrelationToFitness.R` script).

```mermaid
flowchart TB
subgraph sg1["local (1)"]
direction LR
A[("Selected<br/>phenotypes")] --> B("1_eQTLs_data_preparation_pipeline.sh<br/><em>(local)</em>")
C[("Reference<br/>genome")] --> B
B --> D[(".bed")]
B --> E[(".bim")]
B --> F[(".fam")]
B --> G[(".pheno")]
end
subgraph sg2["HPC"]
direction LR
H("3_Gemma.py<br/><em>(HPC array wrapper)</em>") --> I("4_CollectSignificantEQTLs.R<br/><em>(distant script)</em>")
I --> J[("Significant<br/>eQTLs")]
end
subgraph sg3["local (2)"]
direction LR
K("6_merge_eQTLs_pipeline.sh<br/><em>(local)</em>") --> L[("Annotated<br/>significant<br/>eQTLs")]
end
sg1 --> |"Upload<br/>input files<br/>(2_UploadGemmaFiles.py)"| sg2
sg2 --> |"Download<br/>output files<br/>(5_DownloadGemmaFiles.py)"| sg3
```

### 📂 WGCNA modules <a name="scripts_9"></a>

      └── scripts
           └── 9_WGCNA_modules
                └── local
                     ├── 1_ExploreModuleFitnessCorrelations.R
                     ├── 2_PlotModuleExploration.R
                     ├── 3_ComputeModules.R
                     └── 4_PlotSelectedModules.R
                ├── 1_module_exploration_pipeline.sh
                ├── 2_plot_module_exploration_pipeline.sh
                └── 3_module_calculation_pipeline.sh

This task computes optimal gene expression functional modules using WGCNA R-package.
The procedure involves two main steps:
- Estimate the "soft power threshold" as the best compromise between a scale-free structure and correlation to fitness,
- Calculate and save final modules using the optimal threshold.
Running the pipeline is facililated by three shell scripts to **(i)** run the exploration, **(ii)** plot the result to make graphic reading decision, and **(iii)** run the proper calculation of modules.

Associated data folder(s): `./data/tribolium_modules`.

#### ⚙️ `1_ExploreModuleFitnessCorrelations.R` (local):
> This script explores a range of soft power thresholds from 1 to 30 (integer values only) to generate information on the scale-free structure and modules' correlation to fitness.

#### ⚙️ `2_PlotModuleExploration.R` (local):
> This script plots the result of the soft power threshold exploration in order to make graphic-reading decision on the threshold leading to the best compromise between a scale-free structure and modules' correlation to fitness.

#### ⚙️ `3_ComputeModules.R` (local):
> This script computes WGCNA functional modules given a phenotype and the soft power threshold.

#### ⚙️ `4_PlotSelectedModules.R` (local):
> This script generates various plots from the resulting modules.

```mermaid
flowchart TB
subgraph "local"
direction LR
A[("Phenotype")] --> B("1_ExploreModuleFitnessCorrelations.R<br/><em>(local)</em>")
B --> C("2_PlotModuleExploration.R<br/><em>(local)</em>")
C --> D[("Seleted<br/>threshold")]
D --> E("3_ComputeModules.R<br/><em>(local)</em>")
A --> E
E --> F[("WGCNA<br/>modules")]
end
```

### 📂 Allele specific expression (ASE) <a name="scripts_10"></a>

      └── scripts
           └── 10_ASE
                └── hpc
                     ├── 1_ASEReadCounter.py
                     └── 2_DetectSignificantASE.R
                └── local
                     └── 3_SummarizeASEData.R

This task uses GATK software to calculate ASE per sample, detects significant ASE-SNPs using an exact binomial test and FDR normalization, and merges the results.

Associated data folder(s): `./data/tribolium_ase`.

#### ⚙️ `1_ASEReadCounter.py` (local):
> This script uses `GATK ASEReadCounter` method to calculates ASE counts per sample.

#### ⚙️ `2_DetectSignificantASE.R` (local):
> This script detects significant ASE using an exact binomial test and FDR normalization.

#### ⚙️ `3_SummarizeASEData.R` (local):
> This script merges and summarizes the ASE dataset at the population level.

```mermaid
flowchart TB
subgraph sg1["HPC"]
direction LR
A[("BAMs")] --> B("1_ASEReadCounter.py<br/><em>(HPC array wrapper)</em>")
B --> C("2_DetectSignificantASE.R<br/><em>(HPC)</em>")
C --> D[("Significant<br/>ASEs")]
end
subgraph sg2["local"]
direction LR
E("3_SummarizeASEData.R<br/><em>(HPC array wrapper)</em>") --> F[("Summarized<br/>ASEs")]
end
sg1 --> |"Download<br/>ASEs"| sg2
```

### 📂 LD map with Lep-MAP3 <a name="scripts_11"></a>
This task calculates a genetic map (or LD map) using Lep-MAP3 software. The pipeline is complex and will not be described here.
**Please consult the PDF report available in the folder `./doc/Lep-MAP3 report` for a detailed description of the pipeline**.

### 📂 Haplotype blocks <a name="scripts_12"></a>

      └── scripts
           └── 12_haplotype_blocks
                └── local
                     └── 1_ExtractHaplotypeBlocks.R

This task extracts haplotype blocks from the genetic map calculated with Lep-MAP3 (see above).
By default, SNPs are considered to belong to the same haplotype block if their genetic distance is equal.

Associated data folder(s): `./data/tribolium_ld`.

#### ⚙️ `1_ExtractHaplotypeBlocks.R` (local):
> This script extracts haplotype blocks from the genetic map and export the result in a file.

```mermaid
flowchart TB
subgraph sg2["local"]
direction LR
A[("Genetic<br/>map")] --> B("1_ExtractHaplotypeBlocks.R<br/><em>(local)</em>")
B --> C[("Haplotype<br/>blocks")]
end
```

### 📂 FST clusters <a name="scripts_13"></a>
➡️ Work in progress ⬅️

Associated data folder(s): `./data/tribolium_diversity`.

## Analyses <a name="analyses"></a>
Analysis pipelines are related to manuscripts under preparation and will be described here later.

## Data <a name="data"></a>

      └── data
           ├── experiment_data: contains various gathered experimental data
           ├── tribolium_afc: contains the result of AFCs calculations
           ├── tribolium_ase: contains the result of ASE tests
           ├── tribolium_bam: contains BAM files and sample information
           ├── tribolium_counts: contains read counts datasets
           ├── tribolium_diversity: contains datasets related to genetic diversity calculations
           ├── tribolium_eqtl: contains significant eQTLs for all phenotypes
           ├── tribolium_filters: contains various VCF filters
           ├── tribolium_genome: contains annotated genomes
           ├── tribolium_haplotype_blocks: contains haplotype blocks extracted from the genetic map
           ├── tribolium_ld: contains the genetic map
           ├── tribolium_modules: contains WGCNA module files
           ├── tribolium_pedigree: contains pedigrees related to the family stucture of the population
           ├── tribolium_phenotypes: contains calculated phenotypes
           ├── tribolium_snp: contains VCF files containing bi-allelic variant SNPs
           └── tribolium_vcf: contains the raw VCF file obtained from the variant call

# License <a name="license"></a>
➡️ License agreement to be displayed here ⬅️

