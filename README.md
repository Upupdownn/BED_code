# BED: Bayesian End-Motif Decomposition

## Description

**Bayesian End-motif Decomposition (BED)** is a latent variable probabilistic framework designed to address signal obscuration and cross-cohort instability in cfDNA end-motif analysis. Based on the **cfDNA mixture hypothesis**, BED separates observed EDM profiles into dominant physiological background and tumor-associated aberrant components to derive purified aberrant features for robust cancer detection.

## Overview
![Pipeline Diagram](assets/workflow.png)

The BED workflow is an integrated pipeline designed for:
1. **Feature Extraction**: Efficiently extracting end-motif (4-mer) frequencies from cfDNA fragments.
2. **Bayesian Decomposition**: Separating observed EDM profiles into physiological background and aberrant components to derive purified aberrant features together with a training-free Aberrant Weight (AW) indicator for aberrant signal quantification.
3. **Cancer Detection**: Training and evaluating SVM classifiers using both raw EDM frequencies and purified aberrant features for robust cancer detection.
4. **Cross-Cohort Validation**: Assessing diagnostic consistency and robustness across independent validation cohorts.
5. **Downstream Analysis**: Supporting tissue-of-origin localization, prognostic analysis, and biological interpretation of aberrant EDM patterns.



## Project Structure

```text
.
├── BED_workflow.py           # Main orchestration script
├── bin/                      # Core functional scripts
│   ├── 01_bam_to_tsv.py
│   ├── 02_extract_edm_features.py
│   ├── 03_merge_edm_features.py
│   ├── 04_generate_ref_pilot.py
│   ├── 05_bayesian_decompo.py
│   └── 06_svm_train_val.py
├── scripts/                  # other functional scripts
├── examples/                 # Small example dataset for testing
│   ├── frag_file/
│   └── info.tsv
└── environment.yml           # Conda environment configuration
```

## Installation
We recommend using Conda to manage the environment; while the pipeline is fully CPU-compatible for accessibility, a CUDA-enabled PyTorch setup is recommended for significantly faster decomposition.

```bash
# Clone the repository
git clone https://github.com/Upupdownn/BED_code.git
cd BED_code

# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate bed_analysis
```

## Preparation

The BED pipeline supports two input formats for starting the analysis. You can either provide raw alignment files (BAM) or pre-processed fragment files (TSV).

**1. Input Fragment Data**

* **Option A: BAM Files**

  * Standard genomic alignment files are supported.

  * **Requirements:** Files must be sorted and indexed (e.g., sample.bam and sample.bam.bai).

* **Option B: Fragment TSV Files**

  * If you have already extracted fragment information, provide a TSV file with a header and the following columns: chr, start, end, mapq, and strand.

  * **Example Data:** For testing purposes, sample fragment TSV files are provided in the `examples/frag_file/` directory.

**2. Reference Genome (2bit format)**

A `.2bit` file of the reference genome is required for sequence extraction and end-motif frequency calculation. Please choose the version (e.g., hg19 or hg38) that matches your alignment.

To run the provided example dataset, you can download the hg19 reference genome using the following command:

```Bash
# Download hg19 reference genome from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
```

**Note:** Ensure the path to the .2bit file is correctly passed to the --tb_file argument in the workflow script.

## Usage for BED_workflow.py

### 1. Prepare Inputs
* **Input Directory**: A folder containing `.bam` or `.tsv` fragment files (including `.tsv.gz` files).

* **Reference File**: A `.2bit` reference genome file (e.g., hg19.2bit).

* **Info File**: A TSV file with at least `id` and `label` (0 for healthy, 1 for cancer) columns.

### 2. Run the Workflow
The BED_workflow.py script automates the entire process:
```bash
./BED_workflow.py \
    /path/to/input_dir \
    /path/to/output_dir \
    --tb_file /path/to/hg19.2bit \
    --train_info_file train_info.tsv \
    --processes 20
```

### 3. Key Arguments
* `input_dir`: Path to BAM or TSV files.

* `output_dir`: Directory where all results and intermediate files will be saved.

* `--tb_file`: Path to the reference genome in .2bit format.

* `--train_info_file`: TSV containing sample metadata for training.

* `--val_info_file`: (Optional) TSV containing sample metadata for validation.

* `-p`: Number of parallel processes (default: 10).

### Expected Outputs
After a successful run, the output_dir will contain:

* `features/`: Merged EDM and ACW (decomposed) feature matrices.

* `svm_edm/`: SVM models and performance scores based on raw motifs.

* `svm_acw/`: SVM models and performance scores based on BED-decomposed features.

* `pilot/`: Results from the learning rate search and pilot decomposition.

## Sub-scripts Usage
All scripts are located in the scripts and bin directories. Each script provides a command-line help message and can be inspected using the `-h` parameter.


## Contacts

If you have any questions or feedback, please contact us at:

**Email**: upupdownn@gmail.com

## Software License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
