# Z-RNA Project

This repository contains analysis pipelines and utility scripts for processing and analyzing Z-RNA sequencing data. The project includes raw data preprocessing, Z-RNA cluster detection, secondary structure prediction, structure feature extraction, and k-mer enrichment analysis.

## Directory Structure

```
.
├── APP_Index.py
├── GenoTransPosShuttle.pl
├── basic.sh
├── bp_paired.py
├── bp_paired2.py
├── kmer_speed4_0.5.py
├── loop_stem_sequence.py
└── rawdata_process.sh
```

---

## Script Descriptions

### 1. `rawdata_process.sh`

Pipeline for preprocessing paired-end Z-RNA sequencing data. It performs:

* Quality control (using `fastqc`)
* Adapter trimming (`fastp`)
* Genome alignment (`STAR`)
* Filtering reads mapped to plasmid sequences (identified as `chrPlasmid`)
* Strand specificity check (`infer_experiment.py`)

**Dependencies:**

* `fastqc`, `fastp`, `STAR`, `samtools`, `infer_experiment.py`
* Reference genome and plasmid genome (with plasmid chromosome named as `chrPlasmid`)

**Usage:**

* Set `rootpath` to your working directory.
* Set `Datapath` to the folder containing your input fastq.gz files (`*_1.fastq.gz` and `*_2.fastq.gz`).
* Set `refpath` to your reference genome and annotation files.

---

### 2. `basic.sh`

Performs core Z-RNA analysis:

* **Cluster Definition**: Identifies clusters based on RNA editing density (10 editing sites within 100 bp in transcript coordinates ).
* **Secondary Structure Prediction**: Uses tools like `RNAfold`, `RNAstructure`, and SHAPE-seq data to predict RNA secondary structures.

**Dependencies:**

* `RNAfold`
* `RNAstructure`
* `bedtools`
* `GenoTransPosShuttle.pl`
* SHAPE-seq data for structure-assisted folding

---

### 3. `GenoTransPosShuttle.pl`

A Perl script to convert genomic coordinates to transcript coordinates based on GTF/GFF annotation. Useful for integrating structural and expression data at the transcript level.

---

### 4. `bp_paired.py`

Extracts RNA stem (base-paired) regions from `.ct` format secondary structure prediction files.

* Only supports fixed stem lengths of 6–10 bp.
* Outputs paired base positions suitable for further sequence analysis.

---

### 5. `bp_paired2.py`

An enhanced version of `bp_paired.py`:

* Allows users to define a **minimum stem length** threshold.
* More flexible for exploring variable stem structures in RNA.

---

### 6. `kmer_speed4_0.5.py`

Performs **k-mer enrichment analysis** and calculates the **APP Index**:

* Compares two input fasta files (`fastaA`, `fastaB`)
* For each k-mer (user-defined length `k`), calculates:

  * Frequency in both files
  * Fold change (fastaA vs. fastaB)
  * APP index (enrichment score)
* k-mers not found in one file are assigned a pseudo count of 0.5.

---

### 7. `loop_stem_sequence.py`

Extracts **stem and loop sequences** from `.ct` files:

* Parses base-pairing information
* Outputs fasta-formatted sequences for stem and loop regions
* Useful for sequence motif and structure-based analysis

---

## How to Cite

If you use this pipeline or any component in your work, please cite the corresponding publication (under submission, please check back for DOI).

---

## Contact

For questions or suggestions, please contact me at LegolasdeAda@163.com or open an issue.
