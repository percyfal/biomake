# biomake #


Makefile library for bioinformatics programs, with a focus on next-generation sequencing

# Usage #

# Pipelines #

NB: the pipeline makefiles are tailored for use with data generated at
SciLife, Stockholm. This means that data should be organized as
follows (example data available from repository
https://github.com/percyfal/ngs.test.data):

>    |-- P001_101_index3
>    |   |-- 120924_AC003CCCXX
>    |   |   |-- P001_101_index3_TGACCA_L001_R1_001.fastq.gz
>    |   |   |-- P001_101_index3_TGACCA_L001_R2_001.fastq.gz
>    |   |   `-- SampleSheet.csv
>    |   |-- 121015_BB002BBBXX
>    |   |   |-- P001_101_index3_TGACCA_L001_R1_001.fastq.gz
>    |   |   |-- P001_101_index3_TGACCA_L001_R2_001.fastq.gz
>    |   |   `-- SampleSheet.csv
>    |-- P001_102_index6
>    |   `-- 120924_AC003CCCXX
>    |       |-- P001_102_index6_ACAGTG_L002_R1_001.fastq.gz
>    |       |-- P001_102_index6_ACAGTG_L002_R2_001.fastq.gz
>    |       `-- SampleSheet.csv

## Haloplex ##

Copy the file example/Makefile.pipeline.halo to the root directory
where project data resides. Uncomment relevant sections and set the
variables. See the included Makefile (Makefile.halo) for further
options. Example commands follow.

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

### 1. Make entire pipeline ###

Run with '-n' flag to monitor commands:

	make -n all
	
Make all target with one of the following commands:

	make all
	make all.filtered.eval_metrics
	
### 2. Make flowcell targets ###	

	make flowcells
	
### 3. Make sample targets ###	

	make samples
	
### 4. Clean up ###	

The clean target removes everything.

	make clean
	
	
