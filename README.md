# biomake #

Makefile library for bioinformatics programs, with a focus on
next-generation sequencing.

# Disclaimer #

Use the make rules at your own risk, and make sure you understand them
before running any commands. I take no responsibility if you'd happen
to run a **make clean** in an inappropriate location, removing
precious data in the process. You have been warned!

# Introduction  #

The makefiles contain general recipies for commonly used
bioinformatics programs. The use cases reflect the needs I've had and
do by no means have a comprehensive coverage. Nevertheless, many
commands are so commonly used that the recipes may be of general
interest.

## Installation ##

Clone the repository https://github.com/percyfal/biomake to an
appropriate location.

# Usage #

The intended usage is that the user first creates a Makefile for use
with a particular dataset/problem. Thereafter, **include** statements
are used to include Makefiles of interest:

	#-*- makefile -*-
	# Makefile example

	# Include bwa-specific Makefile
	include /path/to/Makefile.bwa

	# Each Makefile has configurable variables
	BWA_THREADS=4
	
Each Makefile has a set of configurable variables that can be
customized via the command line or in the "calling" Makefile. I have
tried to name the variables sensibly, providing namespaces following
the format PROGRAM_OPTIONS (e.g. PICARD_OPTIONS) or
PROGRAM_SUBPROGRAM_OPTIONS (e.g. GATK_SELECTSNPVARIANTS_OPTIONS).
Similarly, for increased granularity, specific options may be named by
appending the name of that option as in PROGRAM_OPTION_SPECIFIC (e.g.
ANGSD_OPTION_ANC).

There are also a set of higher level options that govern general
behaviour, such as THREADS (number of threads) and JAVA_MEM (java
memory). Where applicable, these options also have program specific
options that can be overridden (such as BWA_THREADS and
GATK_JAVA_MEM).

Unfortunately, the only way to know which configuration variables are
available is to sift through the relevant Makefiles.

# Pipelines #

NB: the pipeline makefiles are tailored for use with data generated at
SciLife, Stockholm. This means that data should be organized as
follows (example data available from repository
https://github.com/percyfal/ngs.test.data):

    |-- P001_101_index3
    |   |-- 120924_AC003CCCXX
    |   |   |-- P001_101_index3_TGACCA_L001_R1_001.fastq.gz
    |   |   |-- P001_101_index3_TGACCA_L001_R2_001.fastq.gz
    |   |   `-- SampleSheet.csv
    |   |-- 121015_BB002BBBXX
    |   |   |-- P001_101_index3_TGACCA_L001_R1_001.fastq.gz
    |   |   |-- P001_101_index3_TGACCA_L001_R2_001.fastq.gz
    |   |   `-- SampleSheet.csv
    |-- P001_102_index6
    |   `-- 120924_AC003CCCXX
    |       |-- P001_102_index6_ACAGTG_L002_R1_001.fastq.gz
    |       |-- P001_102_index6_ACAGTG_L002_R2_001.fastq.gz
    |       `-- SampleSheet.csv

## Haloplex ##

Copy the file **example/Makefile.pipeline.halo** to the root directory
where project data resides. Uncomment relevant sections and set the
variables. See the included Makefile (**Makefile.halo**) for further
options. 

The pipeline currently does the following:

1. makes flowcell targets *.sort.rg.bam (phony target flowcells)
    1. adapter trimming
	2. resyncing mates
	3. alignment with bwa mem
	4. sorting and adding of read group information

2. makes sample targets *.sort.merge.realign.recal.clip.bam (phony target samples)
	1. calls raw genotypes to identify realignment target regions
	2. indel realignment
	3. base recalibration
	4. read clipping

3. merges data and makes target all.filtered.eval_metrics (phony target all)
	1. merges samples to all.bam
	2. genotyping
	3. variant filtration with Haloplex-specific hard filters
	4. variant evaluation

TODO: add fastqc and picard metrics 

### Example commands ###

#### 1. Make entire pipeline ####

Run with '-n' flag to monitor commands:

	make -n all
	
Make all target with one of the following commands:

	make all
	make all.filtered.eval_metrics
	
#### 2. Make flowcell targets ####

	make flowcells
	
#### 3. Make sample targets ####

	make samples
	
#### 4. Clean up ####

The clean target removes everything.

	make clean
	
#### 5. Running a specific sample ####

	make all SAMPLES=P001_101_index3

#### 6. Running a specific flowcell ####

TODO.
