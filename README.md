# Meta-Storms 2

![Version](https://img.shields.io/badge/Version-%202.3%20beta-brightgreen.svg)
![Release date](https://img.shields.io/badge/Release%20date-Feb.%204%2C%202020-brightgreen.svg)



## Contents

- [Introduction](#introduction)
- [System Requirement and dependency](#system-requirement-and-dependency)
- [Installation guide](#installation-guide)
- [Pre-computing](#pre-computing)
- [Example dataset](#example-dataset)
- [Build a MSE database](#build-a-mse-database)
- [Search the MSE database](#search-the-mse-database)
- [Multiple classification based on search result](#multiple-classification-based-on-search-result)
- [Microbiome Novelty Score (MNS) based on search results](#microbiome-novelty-score-mns-based-on-search-results)
- [Input and output file format](#input-and-output-file-format)
- [Supplementary](#supplementary)
- [Citation](#citation)
- [Contact](#contact)

## Introduction

Meta-Storms 2 is the standalone implementation of the Microbiome Search Engine (MSE; http://mse.ac.cn). MSE is a search engine designed to efficiently search a database of microbiome samples and identify similar samples based on phylogenetic or functional relatedness. Meta-Storms 2 consists of the following steps: (i) creating a database composed of reference microbiome samples, and (ii) searching for similar samples in the database with given query microbiome sample(s) via phylogenetic similarities. Meta-Storms 2 relies on an advanced indexing algorithm, providing a fast, and constant, search speed in very large databases. Now Meta-Storms 2 supports OTU-based (for 16S rRNA), Species-based (for WGS) and KEGG-Ontology-based (for function) search, which is compatible with profiling tools of Parallel-META 3, QIIME/QIIME2 and MetaPhlAn2. 


## System Requirement and dependency

### Hardware Requirements

Meta-Storms 2 only requires a standard computer with sufficient RAM to support the operations defined by a user. For typical users, this would be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

  RAM: 8+ GB  
  CPU: 4+ cores, 3.3+ GHz/core

### Software Requirements

OpenMP library is the C/C++ parallel computing library. Most Linux releases have OpenMP already been installed in the system. In Mac OS X, to install the compiler that supports OpenMP, we recommend using the Homebrew package manager:
```
brew install gcc --without-multilib
```

## Installation guide

### Automatic Installation (recommended)

At present, Meta-Storms 2 provides a fully automatic installer for easy installation.

a. Download the package:
```
git clone https://github.com/qibebt-bioinfo/meta-storms.git
```

b. Install by installer:
```
cd meta-storms
source install.sh
```

The package should take less than 1 minute to install on a computer with the specifications recommended above.

### Tips for Automatic Installation

1. The automatic installer configures the environment variables to the default configuration specified in the file of "\~/.bashrc" or "\~/.bash_profile". If you prefer to configure the environment variables to other configuration file, please choose the option of manual installation below.

2. If the environment variables are not activated automatically, please enable them manually by running the command "source \~/.bashrc".

3. If the automatic installer fails, Meta-Storms 2 can still be installed manually by the following options.

### Manual Installation

If the automatic installer fails, Meta-Storms 2 can still be installed manually.

a. Download the package:
```
git clone https://github.com/qibebt-bioinfo/meta-storms.git
```

b. Configure the environment variables (the default environment variable configuration file is "\~/.bashrc"):

```
export MetaStorms=Path to Meta-Storms 2
export PATH="$PATH:$MetaStorms/bin/"
source ~/.bashrc
```

c. Compile the source code (this is required **only** when installing the source code package):
```
cd meta-storms
make
```

## Pre-computing

Meta-Storms 2 accepts microbiome samples profiled into OTUs (for 16S) or species (for shotgun) or KEGG Orthologies (KO, for both 16S and shotgun).

### OTU (for 16S sequences)

16S rRNA amplicon sequences can be picked into OTUs against GreenGenes 13-8 (97% level) reference by [Parallel-META 3](https://github.com/qibebt-bioinfo/parallel-meta) (recommended) or [QIIME](http://qiime.org/). For a give sequence file (FASTA or FASTQ format, eg. sample1.fa), OTUs can be profiled from 16S sequences in the three alternative methods:

a. by [Parallel-META 3](https://github.com/qibebt-bioinfo/parallel-meta) (recommended):

```
PM-parallel-meta -f F -r sample1.fa -o sample1.out
```

Then the output file "sample1.out/classification.txt" is qualified as the input for Meta-Storms 2 (refer to [Single sample](#single-sample-file-and-sample-list)). For multiple samples as input, samples should be listed in the sample list (refer to [Sample list](#single-sample-file-and-sample-list)).

b. by QIIME


```
pick_otus.py -m uclust_ref --suppress_new_clusters -i sample1.fa -o sample1.out
MetaDB-parse-qiime-otu -i sample1.out/sample1_otus.txt -o sample1.out/classification.txt
```

Then the output file "sample1.out/classification.txt" is qualified as the input for Meta-Storms 2 (refer to [Single sample](#single-sample-file-and-sample-list)). For multiple samples as input, samples should be listed in the sample list (refer to [Sample list](#single-sample-file-and-sample-list)).

c. by QIIME2: integrate Meta-Storms2 as [QIIME2-plug-in](https://github.com/qibebt-bioinfo/q2-metastorms.git).

### Species (for shotgun sequences)
Metagenomic shotgun sequences can be annotated into species by [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) (recommended). With a input sequence file “sample1.fa”, to get the species using [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) by:
```
metaphlan2.py sample_1.fa --input_type fasta --tax_lev s --ignore_viruses --ignore_eukaryotes --ignore_archaea > profiled_sample_1.sp.txt
```
Then the output file "profiled_sample_1.sp.txt" is qualified as the input for Meta-Storms 2 (refer to [Single sample](#single-sample-file-and-sample-list)). For multiple samples as input, samples should be listed in the sample list (refer to [Sample list](#single-sample-file-and-sample-list)).

### Functions (KEGG Orthology, for both 16S and shotgun)
Both 16S and shotgun sequences can be annotated into KEGG Orthologies (KO) by [Parallel-META 3](https://github.com/qibebt-bioinfo/parallel-meta) (recommended for 16S) or [Humann2](http://www.huttenhower.org/humann) (recommended for shotgun). With a input sequence file "sample1.fa":

a. 16S sequences by Parallel-META 3
```
PM-parallel-meta -r sample1.fa -o sample1.out
```
in the output directory "sample1.out", the result file "functions.txt" is qualified as the input for Meta-Storms 2 (refer to [Single sample](#single-sample-file-and-sample-list)). For multiple samples as input, samples should be listed in the sample list (refer to [Sample list](#single-sample-file-and-sample-list)).

b. Shotgun sequences by Humann2
```
humann2 --input sample.fa --output sample1.out
```
## Example dataset

Here we provide a demo dataset with 20 human oral microbiome samples in two different healthy statuses from *Huang, et al., 2014*. The pre-computing result (in the [OTU table](#otu-table) format and derived from Parallel-META 3 and the meta-data are in the "[example](#example-dataset)" folder in the installation package. We use this dataset to demonstrate all the following example commands.

Please change your work directory to the "example" folder by

```
cd example
sh Readme
```

\* Huang, S., et al., *Predictive modeling of gingivitis severity and susceptibility via oral microbiota*. ISME J, 2014. 8(9): p. 1768-80.

## Build a MSE database

The command "MetaDB-make-otu" builds a new MSE database for Meta-Storms 2 based search from the given samples. Samples are listed in either *(i)* [single sample list](#single-sample-file-and-sample-list) (for Parallel- META 3 format, by -i or -l with optional –p,), or *(ii)* [OTU table](#otu-table) (OTU table format, by -T). It outputs a database file (*.mdb).

**Usage:**
```
	MetaDB-make-otu [-option] value
	
	[Input options]
		-i or -l Input filename list
		-p List file path prefix for '-i' or '-l' [Optional for -i and -l]
	or
		-T (upper) Input OTU table (*.OTU.Count)
	or
		-d (*.mdb) Make the HDD mode data files for a database
	
	[Output options]
		-o Output database name, default is "database.mdb"
		-H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F
	
	[Other options]
		-h Help
```

Example (make sure you are in "[example](#example-dataset)" path):
```
MetaDB-make-otu -T taxa.OTU.Count -o database
```

You can also build a MSE database by species (MetaDB-make-sp) or function (MetaDB-make-func) .

### HDD mode

The HDD (Hard Drive Disk) mode uses the re-encoding technique to minimize the RAM usage for database search (although the mode is slower). When the HDD mode is enabled via "–H t", [MetaDB-make-otu/func/sp](#build-a-mse-database) will generate accessory data named as \*.mdb.hdd under the same directory of the output database (\*.mdb). For extremely large databases (e.g., sample number > 10,000), we strongly recommend users to enable the HDD mode to minimize the RAM consumption.

For an existing database (\*.mdb), HDD mode can also be enabled by making its HDD files via the command below. Then the \*.mdb.hdd would be generated and stored under the same directory as the database.

Example (make sure you are in "[example](#example-dataset)" path):
```
MetaDB-make-otu -d database.mdb
```

### Merge MSE databases

The command "MetaDB-merge" merges two existing databases (\*.mdb) into one.
**Usage:**

```
MetaDB-merge [-option] value

[Input and Output options]
	-1 The 1st database name [Required]
	-2 The 2nd database name [Required]
	-o Merged output database name, default is "database_merge.mdb"
	
[Other options]
	-h Help
```

Example: Here you can make another database named as "database_2.mdb"

```
MetaDB-merge -1 database.mdb -2 database_2.mbd -o database_merged
```

## Search the MSE database

### Search via Meta-Storms 2

The database is built by [MetaDB-make-otu](#build-a-mse-database) (\*.mdb). Meta-Storms 2 supports the index-based query, which features an extremely fast and constant search speed against very large microbiome databases.

The query sample(s) can be provided via either (*i*) [single sample](#single-sample-file-and-sample-list) (for a single sample in Parallel-META 3 format, by -i), or (*ii*) [single sample list](#single-sample-file-and-sample-list) (for multiple samples in Parallel- META 3 format, by -l with optional -p), or (*iii*) [OTU table](#otu-table) (for OTU table format, by -T).

We also recommend users to enable the HDD mode for large databases to minimize the RAM consumption (e.g., sample number > 10,000) (See [HDD mode](#hdd-mode)).

**Usage:**
```
	MetaDB-search-otu [Options] Value
	[Database options]
		-d Database file (*.mdb) [Required]
		-H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F
		-P (upper) Path for the HDD mode data files [Optional for '-H T']
	
	[Input options]
		-i Single input file name
	or
		-l Input filename list
		-p Input List file path prefix for '-l' [Optional for -l]
	or
		-T (upper) Input OTU table (*.Count)
	
	[Output options]
		-o Output file, default is "query.out"
	
	[Advanced options]
		-n Number of the matched sample(s), default is 10
		-m Minimum similarity of the matched sample(s), range (0.0 ~ 1.0], default is 0
		-e If enable the exhaustive search (low speed), T(rue) or F(alse), default is F
		-w Abundance weighted or unweighted, T(rue) or F(alse), default is T
	
	[Other options]
		-t CPU core number, default is auto
		-h Help
```
Example(make sure you are in "[example](#example-dataset)" path):
```
MetaDB-search-otu -d database.mdb -T taxa.OTU.Count -o query.out
```

Meta-Storms 2 also supports search by species or function with the commands "MetaDB-search-sp" and "MetaDB-search-func". 
### Search output

The search result contains a number of matches, each with its sample ID and its similarity score (always between 0 and 1) to the query. In the output, for each of the query samples, all of its matches are listed in tandem in a single line, e.g.

| **#**      | **Query** | **Match** | **Similarity** | **Match** | **Similarity** |
| ---------- | --------- | --------- | -------------- | --------- | -------------- |
| Query: | q_id_0    | ref_id_x  | 0.9823         | ref_id_y  | 0.9758         |
| Query: | q_id_1    | Ref_id_m  | 0.9541         | ref_id_n  | 0.9386         |

In the output above, the first query sample (q_id_0) matches against the reference sample (ref_id_x) with a similarity of 0.9823. In addition, q_id_0 also matches ref_id_v with a similarity of 0.9758. The number of matches is assigned by the parameter -n, and default is 10.

The similarity between query sample(s) and matched sample(s) is a phylogeny-based similarity that is computed using the Meta-Storms 2 scoring function. This algorithm takes the relative abundance of OTUs and their binary phylogeny between two samples as input, and output their quantitative similarities (always between 0 and 1). For high performance and parallel computing, this algorithm is optimized by non-recursive transformation, memory recycling and variable reallocation. Please also refer to "Meta-Storms: efficient search for similar microbial communities based on a novel indexing scheme and similarity score for metagenomic data, *Bioinformatics*, 2012" for details.

## Multiple classification based on search result

### Meta-data prediction

Important features of the query sample, such as the meta-data of habitat, status, etc., can potentially be predicted based on the meta-data of its matches. From a search output generated by [MetaDB-search](#search-the-mse-database), the meta-data of the query sample can be predicted by:

**Usage:**
```
	MetaDB-parse-meta [Options] Value
	[Input and Output options]
		-i Input file name (the output of MetaDB-search) [Required]
		-m Input meta-data file name (meta-data of the database in MetaDB-search) [Required]
		-l Meta-data column, default is 1 (exclude the ID column)
		-o Output file name, default is "query.out.meta"
	
	[Advanced options]
		-r Number of predictd meta-data, default is 1
		-b Base of the similarity in the input file, default is 0
		-n Max number of matches in the input file, default is 10
		-s Number of skipped matches in the input file, default is 0
	
	[Other options]
		-h Help
```
Usage for the -s:

When the query sample has already been included in the database, the search result must contain the query samples itself as the top hit since they have the 100% similarity, which causes the bias in meta-data prediction. Here we can use the parameter –s 1 to exclude this top hit in the meta-data prediction to avoid such bias.

Example(make sure you are in "[example](#example-dataset)" path):
```
MetaDB-parse-meta -i query.out -m meta.txt -o query.out.meta
```
### Multiple classification output

[MetaDB-parse-meta](#meta-data-prediction) generates the predicted meta-data with the assigned scores (always between 0 and 1). In the output, for each of the query samples, all of its predicted meta-data are listed in tandem in a single line, e.g.

| **#ID**    | **Meta-data** | **Score** | **Meta-data** | **Score** |
| ---------- | ------------- | --------- | ------------- | --------- |
| q_id_0 | Healthy       | 0.75      | Disease       | 0.25      |
| q_id_1 | Disease       | 0.72      | Healthy       | 0.28      |

In the output above, the first query sample (q_id_0) is predicted as "Healthy" with a score of 0.75, and "Disease" with a score of 0.25. The predicted meta-data are sorted by their scores.                              

The number of predicted meta-data is assigned by parameter -r, and default is 1 (i.e., only reporting the predicted meta-data with the highest score).

## Microbiome Novelty Score (MNS) based on search results

### Calculate the Microbiome Novelty Score (MNS)

With the search output generated by [MetaDB-search](#search-the-mse-database), the Microbiome Novelty Score (MNS) of each sample can be calculated by:

**Usage:**

	MetaDB-parse-mns [Options] Value
	
	[Input and Output options]
		-i Input file name (the output of MetaDB-search) [Required]
		-o Output file name, default is "query.out.mns"
	
	[Advanced options]
		-b Base of the similarity in the input file, default is 0
		-n Max number of matches in the input file, default is 10
		-s Number of skipped matches in the input file, default is 0
	
	[Other options]
		-h Help
Usage for the -s:

When the query sample has already been included in the database, the search result must contain the query samples itself as the top hit since they have the 100% similarity, which causes bias in calculating the MNS. Here we can use the parameter –s 1 to exclude this top hit in calculating the MNS.

Example (make sure you are in "[example](#example-dataset)" path):
```
MetaDB-parse-mns -i query.out -o query.out.mns
```
### Microbiome Novelty Score (MNS) output

[MetaDB-parse-mns](#microbiome-novelty-score-mns-based-on-search-results) generates the MNS (always between 0 and 1) of each query sample in a single line, e.g.

| **#ID**    | **MNS** |
| ---------- | ------- |
| q_id_0 | 0.06    |
| q_id_1 | 0.12    |

 

In the output above, the first query sample (q_id_0) reports a MNS of 0.06.

## Input and output file format

Meta-Storms 2 accepts the two alternative formats as input.

### Single sample file and sample list

A single sample is the OTU/species/KO information of a single microbiome sample profiled by [Parallel-META 3](https://github.com/qibebt-bioinfo/parallel-meta), [QIIME](http://qiime.org/) or [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) from the 16S/shotgun sequences (refer to [Pre-computing](#pre-computing) for details). It is a plain-text file. An example of the single sample is below:


| **#Database_OTU** | **Count** |
| ----------------- | --------- |
| OTU_1             | 10        |
| OTU_2             | 17        |
| OTU_3             | 38        |

A sample list is a plain-text file for listing multiple samples (by –l) as Meta-Storms 2 input, which consists of two columns: the sample IDs and the directories of samples’ single input files, e.g. 

| Sample_1 | /home/data/single_sample/Sample_2/classification.txt |
| -------- | ---------------------------------------------------- |
| Sample_2 | /home/data/single_sample/Sample_2/classification.txt |

The directory can be either absolute directory or relative directory. Meta-Storms 2 also provides an optional parameter –p to add a prefix for the all the directories in the sample list in case of a relative directory is preferred.

### OTU/species/KO table

An OTU/species/KO table is a plain-text file that contains the features (OTU/species/KO) and their richness for each of multiple samples. An example of the OTU table is bellow:

| **#Sample_ID** | **OTU_1** | **OTU_2** | **OTU_3** | **OTU_4** | **OTU_5** |
| -----------| --------- | --------- | --------- | --------- | --------- |
| Sample_1   | 10        | 17        | 38        | 2         | 2         |
| Sample_2   | 0         | 5         | 57        | 0         | 0         |
| Sample_3   | 2         | 35        | 7         | 0         | 0         |
| Sample_4   | 58        | 30        | 23        | 3         | 0         |
| Sample_5   | 95        | 5         | 5         | 4         | 0         |


The output file format can be found at [Search output](#search-output).

## Supplementary

The test code and datasets for reproducing the results of manuscript "Multiple-Disease Detection and Classification across Cohorts
via Microbiome Search" is available here ([Linux X86_64](http://bioinfo.single-cell.cn/Released_Software/meta-storms/test_package/test_package_linux.tar.gz) / [Mac OS X](http://bioinfo.single-cell.cn/Released_Software/meta-storms/test_package/test_package_mac.tar.gz), ~ 892 MB). See “Readme.txt" in the package for usage and details.

## Citation
X. Su*, G. Jing, D. McDonald, H. Wang, Z. Wang, A. Gonzalez, Z. Sun, S. Huang, J. Navas, R. Knight* and J. Xu*. Identifying and predicting novelty in microbiome studies, *mBio*, 2018.
## Contact

Any problem please contact MSE development team:

**JING Gongchao**&nbsp;&nbsp;&nbsp;&nbsp;Email: jinggc@qibebt.ac.cn
