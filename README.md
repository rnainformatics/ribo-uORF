# Source code for modified RibORF and related programs are available at https://github.com/rnainformatics/ribo-uORF.

> The source code is provided by Liuqi (biolq668@gmail.com). For questions and comments, please contact Liuqi or submit an issue on github.

This source code repository can be separated into three main sections：

> - 1: Modified_RibORF
> - 2: uORFscan: tools for uORF identification from user-loaded Ribo-seq datasets
> - 3: UTR5var: tools for investigating the effect of variation on 5’ UTR regions

### 1: Modified_RibORF

- `ORFannotate.pl`

- `ribORF.pl`

- `uorf_lib_filter.R`

- `collapse_reads_md.pl`

  Convert fastq to collapsed fasta:

  ```perl
  Usage: perl fq2collapedFa.pl -i <fq> -o <fa>
  <fq> : fastq file without adaptor (.fq .fastq).
  <fa> : output fasta file (.fa .fasta).
  ```

- `fq2collapedFa.pl`

### 2: uORFscan: tools for uORF identification from user-loaded Ribo-seq datasets

####  main program:

- `run.pl`   

  ```txt
  Usage: perl run.pl <species> <minlength> <maxlength>
  	<mismatch> <maxMultiMapping> 
  	<detectOffset> <deDup> <ORFpvalue> 
  	<jobid> <uploadfile>
  	
  <species>: Select the reference genome 
  <minlength>: Shortest RPF Length of interval (nt)
  <maxlength>: Longest RPF Length of interval (nt)
  <mismatch>: Allowed mismatch in RPF mapping
  <maxMultiMapping>: Max. of multiple-mapping
  <detectOffset>: Whether detect p-site offset automatically
  <deDup>: whether remove unique molecular identifiers (UMI) for the Ribo-seq data which used UMI to differentiate biological duplicates from PCR duplicates
  <ORFpvalue>: Score cutoff for active translated uORFs
  <jobid>: Random string (16 characters)
  <uploadfile>: fasta file (.fa .fasta)
  ```

#### related programs:

- `bamExpander.pl`
- `fastqparse.pl`
- `filterLength.pl`
- `filterUniqueFa.pl`
- `filter_stat.pl`
- `get_plot_highcharts.R`
- `loadMoreMake.pl`
- `parseFeatures.pl`
- `parsePrice.pl`
- `ribORF.parrel.pl`
- `ribo-meta_web_single.R`

### 3: UTR5var: tools for investigating the effect of variation on 5’ UTR regions

####  main program:
- `run_var.pl` 

  ```txt
  Usage: perl run_var.pl <species> <jobid> <email>
  
  <species>: Select the reference genome
  <jobid>: Random string (16 characters)
  <email>: Get an notification when the job is done 
  ```
####  related programs:
- `checkvcf.pl`

