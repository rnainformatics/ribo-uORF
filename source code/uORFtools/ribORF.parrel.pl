use Data::Dumper;
use Parallel::ForkManager;

if ($#ARGV < 3) {
	print "##ribORF predication"."\n";
	print "usage: perl ribORF.pl -f readCorrectedFile -c candidateORFFile -o outputDir [-l orfLengthCutoff] [-r orfReadCutoff] [-p predictPvalueCutoff]"."\n";
	print "-f readCorrectedFile: input read mapping file;"."\n";
	print "-s offsetfile: input read mapping offset file for offset correction;"."\n";
	print "-d candidateORFdir: directory of candidate ORFs in genePred format;"."\n";
	print "-o outputDir: output directory, with files reporting testing parameters and predicted translating probability;"."\n";
	print "-l orfLengthCutoff [optional]: cutoff of ORF length (nt), default: 6;"."\n";
	print "-r orfReadCutoff [optional]: cutoff of supported read number, default: ;"."\n";
	print "-p predictPvalueCutoff [optional]: cutoff used to select predicted translated ORF, default: 0.7."."\n";
	exit;
}


use Getopt::Std;

### get arguments ###
my %args; 
getopt("fdsolrp",\%args);
my $readfile=$args{f};
if (not -e $readfile) {
	print "No corrected read mapping file"."\n";
    exit;
}
my $offFile=$args{s};
if (not -e $offFile) {
	print "No corrected read mapping file"."\n";
    exit;
}
my $orfdir=$args{d};
if (not -e $orfdir) {
	print "No candidate ORF annotation file"."\n";
    exit;
}
my $outputDir=$args{o};
if (! $outputDir) {
	print "No output directory"."\n";
    exit;
}
my $orfLengthCutoff;
my $orfReadCutoff;
my $predictPvalueCutoff;
if (exists ($args{l})) {
	$orfLengthCutoff=$args{l};
} else {
	$orfLengthCutoff=6;
}
if (exists ($args{r})) {
	$orfReadCutoff=$args{r};
} else {
	$orfReadCutoff=11;
}
if (exists ($args{p})) {
	$predictPvalueCutoff=$args{p};
} else {
	$predictPvalueCutoff=0.5;
}

system  "mkdir -p $outputDir/temp" if ! -d "$outputDir/temp";


my %dist;
open (AN, "$offFile");
while (<AN>) {
	chomp;
	my @s=split /\t/, $_;
	$dist{$s[0]}=$s[1];
}
close AN;


my %read;
my %chrs;
open(IN, "samtools view -F4 $readfile |") || die "$!";
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\D+/, $s1[5];
	my @s3=split /\d+/, $s1[5];
	my $len=0;
	for (my $i=0; $i<$#s3; $i++) {
		if ($s3[$i+1] eq 'M') {
			$len+=$s2[$i];
		}
	}
	if (exists ($dist{$len})) {
		my $asite=$dist{$len};
		
		my $loc=$s1[3];
		my $ra=0;
		my $ind=1;
		if ($s1[1]==0) {
			for (my $i=0; $i<$#s3 && $ind==1; $i++) {
				if ($s3[$i+1] eq 'M') {
					if ($s2[$i] >= ($asite+1)) {
						$loc+=$asite;
						$ind=0;
					} else {
						$loc+=$s2[$i];
						$asite-=$s2[$i];
					}
				} elsif ($s3[$i+1] eq 'N') {
					$loc+=$s2[$i];
				} elsif ($s3[$i+1] eq 'D') {
					$loc+=$s2[$i];
				}
			}
		} else {
			for (my $i=0; $i<$#s3; $i++) {
				$loc+=$s2[$i];
			}
			$loc--;
			for (my $i=($#s3-1); $i>=0 && $ind==1; $i--) {
				if ($s3[$i+1] eq 'M') {
					if ($s2[$i] >= ($asite+1)) {
						$loc-=$asite;
						$ind=0;
					} else {
						$loc-=$s2[$i];
						$asite-=$s2[$i];
					}
				} elsif ($s3[$i+1] eq 'N') {
					$loc-=$s2[$i];
				} elsif ($s3[$i+1] eq 'D') {
					$loc-=$s2[$i];
				}
			}
		}
		#my $out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$loc."\t"."50"."\t"."1M"."\t".$s1[6]."\t".$s1[7]."\t".$s1[8]."\t".substr($s1[9], $asite-1, 1)."\t".substr($s1[10], $asite-1, 1);
		$out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$loc."\t"."50"."\t"."1M"."\t".$s1[6]."\t".$s1[7]."\t".$s1[8]."\t".substr($s1[9], $asite-1, 1)."\tA";
		$chrs{$s1[2]} = 1;
		for (my $i=11; $i<=$#s1; $i++) {
			$out.="\t".$s1[$i];
		}
		
        my @readinfo = split(/_x/,$s1[0]); 
		
		#my @s5=split /\D+/, $s1[5];
		#my @s6=split /\d+/, $s1[5];
		my @s5=split /\D+/, '1M';
		my @s6=split /\d+/, '1M';
		if ($s1[1] == 0) {
			my $k=$s1[2].":"."+".":".($loc-1);
			if(not exists $read{$k}) {
			   $read{$k} = $readinfo[1];
			} else {
			   $read{$k} += $readinfo[1];
			}
			
		#} elsif ($s1[1] == 16) {
		} else {
			for (my $i=0; $i<$#s6; $i++) {
				$loc+=$s5[$i];
			}
			my $k=$s1[2].":"."-".":".($loc-2);
			if(not exists $read{$k}) {
			   $read{$k} = $readinfo[1];
			} else {
			   $read{$k}+=$readinfo[1];
			}
		}
	}	
}
close IN;
=pod
if($org eq "cel_WBcel235") {
   	my @chrom = ("chrII","chrIII","chrIV","chrM","chrV","chrX");
} elsif($org eq "dme_BDGP6") {
   	my @chrom = ("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY");
} elsif($org eq "dre_GRCz11") {
   	my @chrom = ("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr23","chr24","chr25","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9");
} elsif($org eq "hg38") {
   	my @chrom = ("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY");
} elsif($org eq "mm10") {
   	my @chrom = ("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY");
} elsif($org eq "rn6") {
   	my @chrom = ("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY");
}

foreach my $chr @chrom {
	$chrs{$chr} = 1;
}
=cut

my $pm = new Parallel::ForkManager(4); 
foreach my $chr (sort sort_by_chr keys %chrs) {
    print $chr,"\n";
	my $orfFile = "$orfdir/$chr.candidateORF.genepred.txt";
	if(!-e $orfFile) {
	   print "[Warning] $orfFile is not existed!\n";
       next;
	}
	my $pid = $pm->start and next;
	analysis_mul($orfFile, $chr, \%read);
	$pm->finish;
}
$pm->wait_all_children;

my @all_out_chr  = glob "$outputDir/temp/*.input.parameters.txt";
merge_files(\@all_out_chr,"$outputDir/all.input.parameters.txt");


sub analysis_mul {
	my ($orfFile,$chr,$read) = @_;
	open (AN, "$orfFile");
	open (OUT, ">$outputDir/temp/$chr.input.parameters.txt");
	#print OUT "orfID"."\t"."chrom"."\t"."strand"."\t"."codon5"."\t"."codon3"."\t"."length"."\t"."readNum"."\t"."f1"."\t"."f2"."\t"."f3"."\t"."entropy"."\t"."MAXentropy"."\t"."PME"."\t"."codonNum"."\t"."f1max"."\n";
	my $dist=3;
	while (<AN>) {
		chomp;
		my @s1=split /\t/, $_;
		my @s2=split /,/, $s1[8];
		my @s3=split /,/, $s1[9];
		my $len1=0;
		my $tot=0;
		my @val;
		my @per;
		my @post;
		my @fra;
		for (my $i=0; $i<=2; $i++) {
			$per[$i]=0;
		}
		for (my $i=0; $i<=$#s2; $i++) {
			for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
				if ($j >=$s1[5] && $j < $s1[6]) {
					push @post, $j;
				}	
			}
		}
		if ($s1[2] eq '-') {
			@post=reverse(@post);
		}
		if ($#post >= 5) {
		#splice @post, 0, 3;
		#splice @post, -3;
		for (my $m=0; $m<=$#post; $m++) {
			$len1++;
			my $k=$s1[1].":".$s1[2].":".$post[$m]; 
			if (exists ($read{$k})) {
				$tot+=$read{$k};
				$val[int(($len1-1)/$dist)]+=$read{$k};
				$per[($len1-1)%3]+=$read{$k};
			}
		}
		for (my $m=0; $m<=$#post; $m+=3) {
			my $k1=$s1[1].":".$s1[2].":".$post[$m];
			my $k2=$s1[1].":".$s1[2].":".$post[$m+1];
			my $k3=$s1[1].":".$s1[2].":".$post[$m+2];
			my $n1=0;
			my $n2=0;
			my $n3=0;
			if (exists ($read{$k1})) {
				$n1=$read{$k1};
			}
			if (exists ($read{$k2})) {
				$n2=$read{$k2};
			}
			if (exists ($read{$k3})) {
				$n3=$read{$k3};
			}
			if (($n1+$n2+$n3) > 0) {
				$fra[0]++;
				if ($n1 > $n2 && $n1 > $n3) {
					$fra[1]++;
				}
			}
		}
		if ($tot >= $orfReadCutoff && $len1 >= $orfLengthCutoff) {
			my $ent=0;
			my $ten=0;
			my $a=int(($len1+2)/$dist);
			my $b=$tot;
			my $t1=int(($a+$b-1)/$b);
			my @val2;
			for ($i=0; $i<=$#val; $i++) {
				if ($val[$i] > 0) {
					$val2[int($i/$t1)]+=$val[$i];
				}
			}
			for ($i=0; $i<=$#val2; $i++) {
				if ($val2[$i] > 0) {
					my $p=$val2[$i]/($tot);
					$ent+=$p*log(1/$p);
				}
			}
			my $t2=int($a/$t1);
			my $d1=int($b/$t2);
			my $d2=$b%$t2;
			my @va;
			for (my $i=0; $i<$t2; $i++) {
				$va[$i]=$d1;
			}
			for (my $j=0; $j<$d2; $j++) {
				$va[$j]++;
			}
			for (my $i=0; $i<=$#va; $i++) {
				if ($va[$i] > 0) {
					$p=$va[$i]/($tot);
					$ten+=$p*log(1/$p);
				}
			}
			my $per;
			if ($ten == 0) {
				$per=1;
			} else {
				$per=$ent/$ten;
			}
			if ($tot > 0) {
				my $out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$s1[5]."\t".$s1[6]."\t".$len1."\t".$tot;
				$out.="\t".sprintf("%.3f", $per[0]/$tot);
				$out.="\t".sprintf("%.3f", $per[1]/$tot);
				$out.="\t".sprintf("%.3f", $per[2]/$tot);
				$out.="\t".sprintf("%.3f", $ent);
				$out.="\t".sprintf("%.3f", $ten);
				$out.="\t".sprintf("%.3f", $per);
				$out.="\t".$fra[0];
				$out.="\t".sprintf("%.3f", $fra[1]/$fra[0]);
				print OUT $out."\n";
			}
		}
		}
	}
	close AN;
	close OUT;
	print "$chr complelted!\n";
}


open (RI, ">$outputDir/ribORF.learning.R");
print RI "library(data.table)"."\n";
print RI "A <- read.table (\"$outputDir/all.input.parameters.txt\", sep=\"\\t\", header=T)"."\n";
print RI "B1 <- A[grepl('canonical', A[,1]),]"."\n";
print RI "B2 <- A[grepl('internal', A[,1]),]"."\n";
print RI "C1 <- data.frame(cbind(B1, gr=rep(1, nrow(B1))))"."\n";
print RI "C2 <- data.frame(cbind(B2, gr=rep(0, nrow(B2))))"."\n";
print RI "T1 <- rbind(C1[sample(nrow(C1), min(nrow(C1)/3, 1000)),], C2[sample(nrow(C2), min(nrow(C2)/3, 2000)),])"."\n";
print RI "T2 <- rbind(C1[sample(nrow(C1), min(nrow(C1)/3, 1000)),], C2[sample(nrow(C2), min(nrow(C2)/3, 2000)),])"."\n";
print RI "pred <- glm(gr~f1+PME+f1max, data=T1)"."\n";
print RI "E1 <- predict(pred, T2)"."\n";
print RI "E2 <- cbind(E1, T2\$gr)"."\n";
print RI "a <- seq(-0.5,1.5,0.005)"."\n";
print RI "fpr <- array()"."\n";
print RI "tpr <- array()"."\n";
print RI "out <- matrix(NA, nrow=length(a), ncol=7)"."\n";
print RI "for (i in 1:length(a)) {"."\n";
print RI "tp <- sum(E2[,1] > a[i] & E2[,2]>0)"."\n";
print RI "fp <- sum(E2[,1] > a[i] & E2[,2]==0)"."\n";
print RI "tn <- sum(E2[,1] <= a[i] & E2[,2]==0)"."\n";
print RI "fn <- sum(E2[,1] <= a[i] & E2[,2]>0)"."\n";
print RI "fpr[i] <- fp/(fp+tn)"."\n";
print RI "tpr[i] <- tp/(tp+fn)"."\n";
print RI "out[i,] <- c(a[i], tp, fp, tn, fn, fpr[i], tpr[i])"."\n";
print RI "}"."\n";
print RI "pdf (file=\"$outputDir/plot.ROC.curve.pdf\")"."\n";
print RI "plot(fpr, tpr, col=0, main=\"ROC Curve\", xlab=\"False positive rate\", ylab=\"True positive rate\")"."\n";
print RI "lines(fpr, tpr,col=2, lwd=3)"."\n";
print RI "dev.off()"."\n";
print RI "colnames(out) <- c(\"cutoff\", \"True.pos\", \"False.pos\", \"True.neg\", \"False.neg\", \"False.pos.rate\", \"True.pos.rate\")"."\n";
print RI "write.table (out, \"$outputDir/stat.cutoff.txt\", sep=\"\\t\", quote=F, row.names=F)"."\n";
print RI "pred.pvalue <- sprintf(\"%.4f\", predict(pred, A))"."\n";
print RI "A <- as.data.table(A)"."\n";
print RI "A[, orfscore := mapply(function(t1,t2,t3,tot) log2(chisq.test(c(t1 * tot, t2 * tot, t3 * tot),p=c(1/3,1/3,1/3))\$statistic), f1, f2, f3, readNum)]"."\n";
print RI "A[, orfscore := round(orfscore,4)]"."\n"; 
print RI "A[f1 < f2 | f1 < f3, orfscore := orfscore * -1]"."\n"; 
print RI "strength <- fread(paste0('$orfdir','/uORF.merged.annotated.txt'),sep='\t',head=F)\n"; 
print RI "strength <- strength[,.(orfID=V1,strength=V7)]\n"; 
print RI "A <- merge(A,strength,by='orfID',all.x=T)\n"; 
print RI "out2 <- data.frame(cbind(A, pred.pvalue))"."\n";
#print RI "out2 <- out2[grepl('uORF', out2[,1]),]"."\n";
print RI "write.table (out2, \"$outputDir/pred.pvalue.parameters.txt\", sep=\"\\t\", quote=F, row.names=F)"."\n";
close RI;

my $com1="Rscript $outputDir/ribORF.learning.R";
system ($com1);


open (IN, "$outputDir/pred.pvalue.parameters.txt");
my %cluster;
<IN>;
while (<IN>) {
	chomp;
	my @s1=split /\t/, $_;
	if ($s1[$#s1] > $predictPvalueCutoff) {
		my @s2=split /\|/, $s1[0];
		my @s3=split /:/, $s2[2];
		my $k1=$s2[0].":".$s3[2];
		$cluster{$k1}.=$s3[1].":".$s2[4].":".$s1[6]."~";
	}
}
close IN;

my %sel;
foreach my $k (sort keys %cluster) {
	my @s1=split /~/, $cluster{$k};
	my $a1;
	my $a2;
	for (my $i=0; $i<=$#s1; $i++) {
		my @s3=split /:/, $s1[$i];
		if ($s3[1] eq 'ATG') {
			$a1.=$s1[$i]."|";
		} else {
			$a2.=$s1[$i]."|";
		} 
	}
	my $a;
	if ($a1 =~ /\|/) {
		$a=$a1;
	} else {
		$a=$a2;
	} 
	my @s2=split /\|/, $a;
	my $b=0;
	my $d;
	my $ind=0;
	for (my $i=0; $i<=$#s2 && $ind==0; $i++) {
		my @s3=split /:/, $s2[$i];
		if ($i < $#s2) {
			my @s4=split /:/, $s2[$i+1];
			$b=$s4[2];
		} else {
			$b=0;
		}
		if ($s3[2] > $b) {
			$ind=1;
			$d=$s3[0].":".$s3[1];
		}
	}
	$sel{$k.":".$d}=1;	
}

open (AN, "$outputDir/pred.pvalue.parameters.txt");
open (OUT, ">$outputDir/repre.valid.pred.pvalue.parameters.txt");
my $line=<AN>;
print OUT $line;
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\|/, $s1[0];
	my @s3=split /:/, $s2[2];
	my $k=$s2[0].":".$s3[2].":".$s3[1].":".$s2[4];
	if (exists ($sel{$k})) {
		print OUT $_."\n";
	}
}
close AN;
close OUT;

open (AN, "$orfFile");
open (OUT, ">$outputDir/repre.valid.ORF.genepred.txt");
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\|/, $s1[0];
	my @s3=split /:/, $s2[2];
	my $k=$s2[0].":".$s3[2].":".$s3[1].":".$s2[4];
	if (exists ($sel{$k})) {
		print OUT $_."\n";
	}
}
close AN;
close OUT;

#system("rm -fr $outputDir/all.input.parameters.txt");

sub merge_files {
	my ($file_ref,$output) = @_;
	my %file = ();
	foreach my $file(@{$file_ref}) {
		my ($chr) = $file =~/(\w+?)\.input.parameters.txt/;
		$file{$chr}= $file;
	}
	my @files = map{$file{$_}} sort {$a cmp $b} keys %file;
	#print Dumper(\@files);
	system("echo \"orfID\tchrom\tstrand\tcodon5\tcodon3\tlength\treadNum\tf1\tf2\tf3\tentropy\tMAXentropy\tPME\tcodonNum\tf1max\n\" > $output");
	system("cat @files >> $output");
	
	#foreach my $file(@files){
	#	unlink "$file";
	#}
}

sub sort_by_chr {  
    my $typeA=0;
	my $typeB=0;
	my $aa;
	my $bb;
    if($a =~ /chr(\d+)/i){
	   $aa = $1;
	   $typeA = 1;
	}
	if($b =~ /chr(\d+)/i){
	   $bb = $1;
	   $typeB = 1;
	}
	if($typeA == 1 and $typeB == 1) {
	   return $aa<=>$bb;
	} else {
	   return $a cmp $b;
	}
}
