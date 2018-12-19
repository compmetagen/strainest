CWL implementation of the StrainEst pipeline
============================================
We have implemented a basic version of the StrainEst pipeline using the 
Common Workflow Language standard (https://www.commonwl.org/)

Usage:

  .. code-block:: sh
  
        cwl-runner ~/path-to-workflow/strainest.cwl --reference_dir reference_dir \
        --reference_basename reference_basename \
        --bowtie2_read1 bowtie2_read1 \
        --bowtie2_read2 bowtie2_read2 \
        ---snv snv  \
        --strainest_est_output_dir_name strainest_est_output_dir_name

where:
"reference_dir" is the directory where the bowtie2-indexed reference database is located;

"reference_basename" is the basename of the bowtie2-indexed database;

"bowtie2_read1" is the forward read file in fastq or fastq.gz format

"bowtie2_read2" is the reverse read file in fastq or fastq.gz format

"snv" is the reference SNV file

"strainest_est_output_dir_name" is the name of the output directory.

The CWL implementation of the pipeline performs the following steps: 
i) alignment of metagenomic reads on the reference database suing bowtie2; 
ii) conversion of the sam file into bam, sorting and indexing;
iii) estimation of the relative abundace of strains using the "strainest est" 
subcommand.