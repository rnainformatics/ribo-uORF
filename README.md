# Source code for modified RibORF and related programs are available at https://github.com/rnainformatics/ribo-uORF.

> The source code is provided by Liuqi (biolq668@gmail.com). For questions and comments, please contact Liuqi or submit an issue on github.

This source code repository can be separated into three main sections：

> - 1: Modified_RibORF
> - 2: uORFscan
> - 3: UTR5var

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

### 2: uORFscan

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
  <jobid>: random string (16 characters)
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

### 3: UTR5var

####  main program:
- `run_var.pl` 

  ```txt
  Usage: perl run_var.pl <species> <jobid> <email>
  ```
####  related programs:
- `checkvcf.pl`

