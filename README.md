# Meta-Storms 2

![Version](https://img.shields.io/badge/Version-std%202.2.1-brightgreen.svg)
![Release date](https://img.shields.io/badge/Release%20date-Jan.%204%2C%202019-brightgreen.svg)



## Contents

- [Introduction](#introduction)
- [Download](#download)
- [Packages](#packages)
- [System Requirement and dependency](#system-requirement-and-dependency)
- [Installation guide](#installation-guide)
- [Notice before use](#notice-before-use)
- [Pre-computing](#pre-computing)
- [Example dataset](#example-dataset)
- [Build a MSE database](#build-a-mse-database)
- [Search the MSE database](#search-the-mse-database)
- [Multiple classification based on search result](#multiple-classification-based-on-search-result)
- [Microbiome Novelty Score (MNS) based on search results](#microbiome-novelty-score-mns-based-on-search-results)
- [Microbiome Attention Score (MAS) based on search results](#microbiome-attention-score-mas-based-on-search-results)
- [File Format](#file-formatsupplementary)
- [Contact](#contact)

## Introduction

Meta-Storms 2 is the standalone implementation of the Microbiome Search Engine (MSE; http://mse.single-cell.cn). MSE is a search engine designed to efficiently search a database of microbiome samples and identify similar samples based on phylogenetic or functional relatedness. Meta-Storms 2 consists of the following steps: (i) creating a database composed of reference microbiome samples, and (ii) searching for similar samples in the database with given query microbiome sample(s) via phylogenetic similarities. Meta-Storms 2 relies on an advanced indexing algorithm, providing a fast, and constant, search speed in very large databases. 

## Download

The latest release is available at:

<http://mse.single-cell.cn/>

## Packages

At present, Meta-Storms 2 provides two alternative packages for installation.

### Prebuilt binary package (recommended)

Meta-Storms 2 prebuilt binary package, with all the tools integrated, is available for Linux (64 bit) and Mac OS X.

### Source code package

Meta-Storms 2 source code package is also available for building and installation for other Unix/Linux/Mac OS X based operating systems.

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

a. Extract the package:
```
tar -xzvf meta-storms-2-std-bin.tar.gz
```

b. Install by installer:
```
cd meta-storms-2-std
source install.sh
```

The package should take less than 1 minute to install on a computer with the specifications recommended above.

### Tips for Automatic Installation

1. Please "cd meta-storms-2-std” directory, before running the automatic installer.

2. The automatic installer configures the environment variables to the default configuration specified in the file of "\~/.bashrc" or "\~/.bash_profile". If you prefer to configure the environment variables to other configuration file, please choose the option of manual installation below.

3. If the environment variables are not activated automatically, please enable them manually by running the command "source \~/.bashrc".

4. If the automatic installer fails, Meta-Storms can still be installed manually by the following options.

### Manual Installation

If the automatic installer fails, Meta-Storms 2 can still be installed manually.

a. Extract the package:
```
tar –xzvf meta-storms-2-std-src.tar.gz
```

b. Configure the environment variables (the default environment variable configuration file is "\~/.bashrc”):

```
export MetaStorms=Path to Meta-Storms 2
export PATH=”$PATH:$MetaStorms/bin/”
source ~/.bashrc
```

c. Compile the source code (this is required **only** when installing the source code package):
```
cd meta-storms-2-std
make
```

## Notice before use

1. For source code package based installation, please make sure proper versions of compilers have been installed: gcc 4.4 or higher for Linux / gcc-8 or higher for Mac OS X (refer to Software Requirements).

2. Meta-Storms 2 optionally accepts microbiome sample(s) that are pre-processed by [Parallel-META 3 ](http://bioinfo.single-cell.cn/parallel-meta.html)(version 3.2 or hihger) or [QIIME](http://qiime.org/) (version 1.9.1). Meta-Storms 2 software can also accept OTU tables (refer to[ **File format**](#file-formatsupplementary)). However, if starting from DNA sequences, [Parallel-META 3](http://bioinfo.single-cell.cn/parallel-meta.html)and [QIIME](http://qiime.org/) are the recommended software for converting amplicon sequences to OTU tables (refer to [**Pre-computing**](#pre-computing)).

3. Make sure that Meta-Storms 2 has the write permission in the output path.

4. We strongly recommend reading this manual carefully before using Meta-Storms 2.

## Pre-computing

To use Meta-Storms, all sequences of microbiome samples must be pre-computed and profiled against the Greengenes database (version 13-8) by [Parallel-META 3](http://bioinfo.single-cell.cn/parallel-meta.html)(version 3.2 or higher) or [QIIME](http://qiime.org/) (version 1.9.1). Then the profiling results will be used as input to Meta-Storms 2.

### Pre-computing by Parallel-META 3

For a give sequence file (FASTA or FASTQ format, eg. sample1.fa), to convert the sequences to OTUs by [Parallel-META 3](http://bioinfo.single-cell.cn/parallel-meta.html):

​	***PM-parallel-meta -f F -r sample1.fa -o sample1.out***

Then the output file *sample1.out/classification.txt* is qualified as the input for Meta-Storms 2 (refer to [**Single sample**](#single-sample-file-and-sample-list)). For multiple samples as input, samples should be listed in the sample list (refer to [**Sample list**](#single-sample-file-and-sample-list)).

### Pre-computing by QIIME

For a give sequence file (FASTA format, eg. sample1.fa), to convert the sequences to OTUs by [QIIME](http://qiime.org/):

​	***pick_otus.py -m uclust_ref --suppress_new_clusters -i sample1.fa -o sample1.out***

​	***MetaDB-parse-qiime-otu -i sample1.out/sample1_otus.txt -o sample1.out/classification.txt***

Then the output file *sample1.out/classification.txt* is qualified as the input for Meta-Storms 2 (refer to [**Single sample**](#single-sample-file-and-sample-list)**)**. For multiple samples as input, samples should be listed in the sample list (refer to [**Sample list**](#single-sample-file-and-sample-list)).

## Example dataset

Here we provide a demo dataset with 20 human oral microbiome samples in two different healthy statuses from *Huang, et al., 2014**. The pre-computing result (in the [**OTU table**](#otu-table)format and derived from Parallel-META 3 and the meta-data are in the "[**example**](#example-dataset)” folder in the installation package. We use this dataset to demonstrate all the following example commands.

Please change your work directory to the "**example**” folder by

​	***cd example***

\* Huang, S., et al., *Predictive modeling of gingivitis severity and susceptibility via oral microbiota*. ISME J, 2014. 8(9): p. 1768-80.

## Build a MSE database

### Build a MSE database by OTU

The command of **MetaDB-make-otu** builds a new MSE database for Meta-Storms 2 based search from the given samples. Samples are listed in either *(i)* [**single sample list**](#single-sample-file-and-sample-list)(for Parallel- META 3 format, by -i or -l with optional –p,), or *(ii)* [**OTU table**](#otu-table)(OTU table format, by -T). It outputs a database file (**.mdb*).

**Usage:**

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
Example (make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-make-otu -T taxa.OTU.Count -o database***

### Build a MSE database by function

The command of **MetaDB-make-func** builds a new MSE database for Meta-Storms 2 based search from the given samples. Samples are listed in either *(i)* [**single sample list**](#single-sample-file-and-sample-list)(for Parallel- META 3 format, by -i or -l with optional –p,), or *(ii)* [**KO table**](#ko-table)(KO table format, by -T). It outputs a database file (**.mdbf*).

**Usage:**

	MetaDB-make-func [Options] Value
	
	[Input options]
		-i or -l Input filename list
		-p List file path prefix for '-i' or '-l' [Optional for -i and -l]
	or
		-T (upper) Input KO table (*.KO.Count)
	or
		-d (*.mdbf) Make the HDD mode data files for a database
	
	[Output options]
		-o Output database name, default is "database.mdbf"
		-H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F
	
	[Other options]
		-h Help
Example (make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-make-func -T taxa.KO.Count -o database***

### Build a MSE database by species

The command of **MetaDB-make-sp** builds a new MSE database for Meta-Storms 2 based search from the given samples. Samples are listed in either *(i)* [**single sample list**](#single-sample-file-and-sample-list)(for Parallel- META 3 format, by -i or -l with optional –p,), or *(ii)* [**OTU table**](#otu-table)(OTU table format, by -T). It outputs a database file (**.mdbs*).

**Usage:**

	MetaDB-make-sp [Options] Value
	
	[Input options]
		-i or -l Input filename list
		-p List file path prefix for '-i' or '-l' [Optional for -i and -l]
	or
	  	-T (upper) Input OTU table (*.OTU.Count)
	or
	  	-d (*.mdbs) Make the HDD mode data files for a database
	
	[Output options]
	  	-o Output database name, default is "database.mdbs"
	  	-H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F
	
	[Other options]
	  	-h Help
Example (make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-make-sp -T taxa.OTU.Count -o database***

#### HDD mode

The HDD (Hard Drive Disk) mode uses the re-encoding technique to minimize the RAM usage for database search (although the mode is slower). When the HDD mode is enabled via –H t, **MetaDB-make-out/func/sp** will generate accessory data named as **.mdb.hdd* under the same directory of the output database (**.mdb*). For extremely large databases (e.g., sample number > 10,000), we strongly recommend users to enable the HDD mode to minimize the RAM consumption.*

For an existing database (**.mdb*), HDD mode can also be enabled by making its HDD files via the command below. Then the **.mdb.hdd* would be generated and stored under the same directory as the database.

Example (make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-make-otu -d database.mdb***

​	***MetaDB-make-func -d dat  abase.mdbf***

​	***MetaDB-make-sp -d database.mdbs***

### Merge MSE databases

The command of **MetaDB-merge** merges two existing databases (**.mdb*) into one.

**Usage:**

​		

```
MetaDB-merge [-option] value

[Input and Output options]
	-1 The 1st database name [Required]
	-2 The 2nd database name [Required]
	-o Merged output database name, default is "database_merge.mdb"
	
[Other options]
	-h Help
```

Example: Here you can make another database named as "*database_2.mdb*"

​	***MetaDB-merge -1 database.mdb -2 database_2.mbd -o database_merged***


## **Search the MSE database**

### Search via Meta-Storms 2 by OTU

Query sample(s) should also be pre-computed by [Parallel-META 3 ](http://bioinfo.single-cell.cn/parallel-meta.html) or [QIIME](http://qiime.org/) using the Greengenes database as reference (refer to [**Pre-computing**](#pre-computing)). The database is built by **MetaDB-make-otu** (**.mdb*). Meta-Storms 2 supports the index-based query, which features an extremely fast and constant search speed against very large microbiome databases.

The query sample(s) can be provided via either (*i*) [**single sample**](#single-sample-file-and-sample-list)(for a single sample in Parallel-META 3 format, by -i), or (*ii*) [**single sample list **](#single-sample-file-and-sample-list)(for multiple samples in Parallel- META 3 format, by -l with optional -p), or (*iii*) [**OTU table**](#otu-table)(for OTU table format, by -T).

We also recommend users to enable the HDD mode for large databases to minimize the RAM consumption (e.g., sample number > 10,000) (See [**HDD mode**](#hdd-mode)).

**Usage:**

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
Example(make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-search-otu -d database.mdb -T taxa.OTU.Count -o query.out***

### Search via Meta-Storms 2 by function

Query sample(s) should also be pre-computed by [Parallel-META 3](http://bioinfo.single-cell.cn/parallel-meta.html)or [QIIME](http://qiime.org/) using the Greengenes database as reference (refer to [**Pre-computing**](#pre-computing)). The database is built by **MetaDB- make-func** (**.mdbf*). Meta-Storms 2 supports the index-based query, which features an extremely fast and constant search speed against very large microbiome databases.

The query sample(s) can be provided via either (*i*) [**single sample**](#single-sample-file-and-sample-list)(for a single sample in [Parallel-META 3 ](http://bioinfo.single-cell.cn/parallel-meta.html) format, by -i), or (*ii*) [**single sample list**](#single-sample-file-and-sample-list)(for multiple samples in Parallel- META 3 format, by -l with optional -p), or (*iii*) [**KO table**](#ko-table)(for KO table format, by -T).

We also recommend users to enable the HDD mode for large databases to minimize the RAM consumption (e.g., sample number > 10,000) (See [**HDD mode**](#hdd-mode)).

**Usage:**

	MetaDB-search-func [Options] Value
	[Database options]
		-d Database file (*.mdbf) [Required]
		-H (upper) If enable the HDD mode (low RAM usage), T(rue) or F(alse), default is F
		-P (upper) Path for the HDD mode data files [Optional for '-H T']
	
	[Input options]
		-i Single input file name
	or
		-l Input filename list
		-p Input List file path prefix for '-l' [Optional for -l]
	or
		-T (upper) Input KO table (*.KO.Count)
	
	[Output options]
		-o Output file, default is "query.func.out"
	
	[Advanced options]
		-n Number of the matched sample(s), default is 10
		-m Minimum similarity of the matched sample(s), range (0.0 ~ 1.0], default is 0
		-e If enable the exhaustive search (high accuracy, low speed), T(rue) or F(alse), default is F
		-M (upper) Distance Metric, 0: Meta-Storms-Func; 1: Cosine; 2: Euclidean; 3: Jensen-Shannon, 4: Bray-Curtis; default is 0
	
	[Other options]
		-t CPU core number, default is auto
		-h Help	
Example(make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-search-func -d database.mdbf -T taxa.KO.Count -o query.func.out***

### Search via Meta-Storms 2 by species

Query sample(s) should also be pre-computed by [Parallel-META 3 ](http://bioinfo.single-cell.cn/parallel-meta.html) or [QIIME](http://qiime.org/)using the Greengenes database as reference (refer to [**Pre-computing**](#pre-computing)). The database is built by **MetaDB-make-sp** (**.mdbs*). Meta-Storms 2 supports the index-based query, which features an extremely fast and constant search speed against very large microbiome databases.

The query sample(s) can be provided via either (*i*) [**single sample**](#single-sample-file-and-sample-list)(for a single sample in [Parallel-META 3 ](http://bioinfo.single-cell.cn/parallel-meta.html) format, by -i), or (*ii*) [**single sample list**](#single-sample-file-and-sample-list)(for multiple samples in Parallel- META 3 format, by -l with optional -p), or (*iii*) [**species table**](#species-table)(for species table format, by -T).

We also recommend users to enable the HDD mode for large databases to minimize the RAM consumption (e.g., sample number > 10,000) (See [**HDD mode**](#hdd-mode)).

**Usage:**

	MetaDB-search-sp [Options] Value
	[Database options]
		-d Database file (*.mdbs) [Required]
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
Example(make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-search-sp -d database.mdbs -T taxa.sp.Count -o query.sp.out***

### Search output

[**MetaDB-search**](#search-the-mse-database)generates a number of matches, each with its sample ID and its similarity score (always between 0 and 1) to the query. In the output, for each of the query samples, all of its matches are listed in tandem in a single line, e.g.

| **#**      | **Query** | **Match** | **Similarity** | **Match** | **Similarity** |
| ---------- | --------- | --------- | -------------- | --------- | -------------- |
| **Query:** | q_id_0    | ref_id_x  | 0.9823         | ref_id_y  | 0.9758         |
| **Query:** | q_id_1    | Ref_id_m  | 0.9541         | ref_id_n  | 0.9386         |

In the output above, the first query sample (q_id_0) matches against the reference sample (ref_id_x) with a similarity of 0.9823. In addition, q_id_0 also matches ref_id_v with a similarity of 0.9758. The number of matches is assigned by the parameter -n, and default is 10.

The similarity between query sample(s) and matched sample(s) is a phylogeny-based similarity that is computed using the Meta-Storms scoring function. This algorithm takes the relative abundance of OTUs and their binary phylogeny between two samples as input, and output their quantitative similarities (always between 0 and 1). For high performance and parallel computing, this algorithm is optimized by non-recursive transformation, memory recycling and variable reallocation. Please also refer to "Meta-Storms: efficient search for similar microbial communities based on a novel indexing scheme and similarity score for metagenomic data, *Bioinformatics*, 2012” for details.

## Multiple classification based on search result

### Meta-data prediction

Important features of the query sample, such as the meta-data of habitat, status, etc., can potentially be predicted based on the meta-data of its matches. From a search output generated by [**MetaDB-search**](#search-the-mse-database), the meta-data of the query sample can be predicted by:

**Usage:**

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
Usage for the -s:

When the query sample has already been included in the database, the search result must contain the query samples itself as the top hit since they have the 100% similarity, which causes the bias in meta-data prediction. Here we can use the parameter –s 1 to exclude this top hit in the meta-data prediction to avoid such bias.

Example(make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-parse-meta -i query.out -m meta.txt -o query.out.meta***

### Multiple classification output

[**MetaDB-parse-meta**](#meta-data-prediction)generates the predicted meta-data with the assigned scores (always between 0 and 1). In the output, for each of the query samples, all of its predicted meta-data are listed in tandem in a single line, e.g.

| **#ID**    | **Meta-data** | **Score** | **Meta-data** | **Score** |
| ---------- | ------------- | --------- | ------------- | --------- |
| **q_id_0** | Healthy       | 0.75      | Disease       | 0.25      |
| **q_id_1** | Disease       | 0.72      | Healthy       | 0.28      |

In the output above, the first query sample (q_id_0) is predicted as "Healthy” with a score of 0.75, and "Disease” with a score of 0.25. The predicted meta-data are sorted by their scores.                              

The number of predicted meta-data is assigned by parameter -r, and default is 1 (i.e., only reporting the predicted meta-data with the highest score).

## Microbiome Novelty Score (MNS) based on search results

### Calculate the Microbiome Novelty Score (MNS)

With the search output generated by [**MetaDB-search**](#search-the-mse-database), the Microbiome Novelty Score (MNS) of each sample can be calculated by:

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

Example (make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-parse-mns -i query.out -o query.out.mns***

### Microbiome Novelty Score (MNS) output

[**MetaDB-parse-mns**](#microbiome-novelty-score-mns-based-on-search-results)generates the MNS (always between 0 and 1) of each query sample in a single line, e.g.

| **#ID**    | **MNS** |
| ---------- | ------- |
| **q_id_0** | 0.06    |
| **q_id_1** | 0.12    |

 

In the output above, the first query sample (q_id_0) reports a MNS of 0.06.

## Microbiome Attention Score (MAS) based on search results

### Calculate the Microbiome Attention Score (MAS)

With the search output generated by [**MetaDB-search**](#search-the-mse-database), the Microbiome Attention Score (MAS) of each sample can be calculated by:

**Usage:**

	MetaDB-parse-mas [Options] Value
	
	[Input and Output options]
		-i Input file name (the output of MetaDB-search) [Required]
		-o Output file name, default is "query.out.mas"
	
	[Advanced options]
		-b Base of the similarity in the input file, default is 0
		-n Max number of matches in the input file, default is 10
		-s Number of skipped matches in the input file, default is 0
		-m Input Meta-data file name to eliminate dup-bias [Optional]
		-l Meta-data column, default is 1 (exclude the ID column) [for -m]
	
	[Other options]
		-h Help

Usage for the -s:

When the query sample has already been included in the database, the search result must contain the query samples itself as the top hit since they have the 100% similarity, which causes bias in calculating the MAS. Here we can use the parameter –s 1 to exclude this top hit in calculating the MAS.

Example (make sure you are in "[**example**](#example-dataset)” path):

​	***MetaDB-parse-mas -i query.out -o query.out.mas***

### Microbiome Attention Score (MAS) output

[**MetaDB-parse-mas**](#microbiome-attention-score-mas-based-on-search-results)generates the MAS of each query sample in a single line, e.g.

| **#ID**    | **MAS**  |
| ---------- | -------- |
| **q_id_0** | 25.05690 |
| **q_id_1** | 17.95170 |

In the output above, the first query sample (q_id_0) reports a MAS of 25.05690.

## File format(supplementary)

Meta-Storms 2 accepts the alternative two formats as input.

### Single sample file and sample list

A single sample is the OTUs and taxonomy information of a single microbiome sample profiled by [Parallel-META 3](http://bioinfo.single-cell.cn/parallel-meta.html)or [QIIME](http://qiime.org/)from the amplicon sequences (refer to [**Pre-computing**](#pre-computing) for details). It is a plain-text file, normally named as "*classification.txt*”. An example of the single sample is below:

| **#Database_OTU** | **Count** |
| ----------------- | --------- |
| OTU_1             | 10        |
| OTU_2             | 17        |
| OTU_3             | 38        |

A sample list is a plain-text file for listing multiple samples (by –l) as Meta-Storms 2 input, which consists of two columns: the sample IDs and the directories of samples’ "*classification.txt*” files, e.g. 

| **Sample_1** | /home/data/single_sample/Sample_2/classification.txt |
| ------------ | ---------------------------------------------------- |
| **Sample_2** | /home/data/single_sample/Sample_2/classification.txt |

The directory can be either absolute directory or relative directory. Meta-Storms 2 also provides an optional parameter –p to add a prefix for the all the directories in the sample list in case of a relative directory is preferred.

### OTU table

An OTU table is a plain-text file that contains the OTUs and their sequence numbers for each of multiple samples. An example of the OTU table is bellow

| **#Sample_ID** | **OTU_1** | **OTU_2** | **OTU_3** | **OTU_4** | **OTU_5** |
| -------------- | --------- | --------- | --------- | --------- | --------- |
| **Sample_1**   | 10        | 17        | 38        | 2         | 2         |
| **Sample_2**   | 0         | 5         | 57        | 0         | 0         |
| **Sample_3**   | 2         | 35        | 7         | 0         | 0         |
| **Sample_4**   | 58        | 30        | 23        | 3         | 0         |
| **Sample_5**   | 95        | 5         | 5         | 4         | 0         |

### KO table

An KO table is a plain-text file that contains the KOs and their sequence numbers for each of multiple samples. An example of the KO table is bellow

| **#Sample_ID** | **KO_1** | **KO_2** | **KO_3** | **KO_4** | **KO_5** |
| -------------- | --------- | --------- | --------- | --------- | --------- |
| **Sample_1**   | 10        | 17        | 38        | 2         | 2         |
| **Sample_2**   | 0         | 5         | 57        | 0         | 0         |
| **Sample_3**   | 2         | 35        | 7         | 0         | 0         |
| **Sample_4**   | 58        | 30        | 23        | 3         | 0         |
| **Sample_5**   | 95        | 5         | 5         | 4         | 0         |

### Species table

An species table is a plain-text file that contains the Species and their sequence numbers for each of multiple samples. An example of the species table is bellow

| **#Sample_ID** | **sp_1** | **sp_2** | **sp_3** | **sp_4** | **sp_5** |
| -------------- | --------- | --------- | --------- | --------- | --------- |
| **Sample_1**   | 10        | 17        | 38        | 2         | 2         |
| **Sample_2**   | 0         | 5         | 57        | 0         | 0         |
| **Sample_3**   | 2         | 35        | 7         | 0         | 0         |
| **Sample_4**   | 58        | 30        | 23        | 3         | 0         |
| **Sample_5**   | 95        | 5         | 5         | 4         | 0         |
## Contact

Any problem please contact MSE development team:

​	**JING Gongchao**&nbsp;&nbsp;&nbsp;&nbsp;Email: jinggc@qibebt.ac.cn
