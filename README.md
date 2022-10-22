# Source code for modified RibORF and related programs are available at https://github.com/rnainformatics/ribo-uORF.

> The source code is provided by Liuqi (biolq668@gmail.com). For questions and comments, please contact Liuqi or submit an issue on github.

http://rnainformatics.org.cn/RiboUORF/ 

### 1: `ORFannotate.pl`

### 2: `ribORF.pl`

### 3: `uorf_lib_filter.R`

### 4: `collapse_reads_md.pl`

Convert fastq to collapsed fasta:
    Usage: `perl fq2collapedFa.pl -i <fq> -o <fa>`
    `<fq>` : fastq file without adaptor (.fq .fastq).
    `<fa>` : output fasta file (.fa .fasta).

### 5: `fq2collapedFa.pl`

### 6：uORFtools

```
bamExpander.pl
fastqparse.pl
filterLength.pl
filterUniqueFa.pl
filter_stat.pl
get_plot_highcharts.R
loadMoreMake.pl
parseFeatures.pl
parsePrice.pl
ribORF.parrel.pl
ribo-meta_web_single.R
run.pl   # main program
```

### 7：UTR5var

```
checkvcf.pl
run_var.pl # main program
```

