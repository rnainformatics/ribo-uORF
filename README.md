# Source code for modified RibORF and related programs are available at https://github.com/rnainformatics/ribo-uORF.

> The source code is provided by Liuqi (biolq668@gmail.com). For questions and comments, please contact Liuqi or submit an issue on github.

### 1: Modified_RibORF

- `ORFannotate.pl`

- `ribORF.pl`

- `uorf_lib_filter.R`

- `collapse_reads_md.pl`

  Convert fastq to collapsed fasta:
      Usage: `perl fq2collapedFa.pl -i <fq> -o <fa>`
      `<fq>` : fastq file without adaptor (.fq .fastq).
      `<fa>` : output fasta file (.fa .fasta).

- `fq2collapedFa.pl`

### 2: uORFscan

####  main program:

- `run.pl`   

  Usage:
  
  ```perl
  perl run.pl <species> <minlength> <maxlength>
  <mismatch> <maxMultiMapping> 
  <detectOffset> <deDup> <ORFpvalue> 
  <email> <jobid> <uploadfile>
  ```

#### related program:

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

####  related program:
- `run_var.pl` 

  Usage:

  ```perl
  perl run_var.pl <species> <jobid> <email>
  ```
####  other program:
- `checkvcf.pl`

