#!/usr/bin/perl -w

use Getopt::Std;
use Parallel::ForkManager;
use Data::Dumper;
use JSON qw( to_json );

if ($#ARGV < 3) {
	print "##ribORF predication"."\n";
	print "usage: perl ribORF.pl -f readCorrectedFile -s offsetFile -c candidateORFFile -o outputDir [-l orfLengthCutoff] [-r orfReadCutoff] [-p predictPvalueCutoff]"."\n";
	print "-f readCorrectedFile: input read mapping file;"."\n";
	print "-s offsetfile: input read mapping offset file for offset correction;"."\n";
	print "-d candidateORFFile: candidate ORFs in genePred format;"."\n";
	print "-o outputDir: output directory, with files reporting testing parameters and predicted translating probability;"."\n";
	print "-t numberOfThreads [optional]: number of threads to calling ORF, default: 3."."\n";
	print "-l orfLengthCutoff [optional]: cutoff of ORF length (nt), default: 6;"."\n";
	print "-r orfReadCutoff [optional]: cutoff of supported read number, default: ;"."\n";
	print "-p predictPvalueCutoff [optional]: cutoff used to select predicted translated ORF, default: 0.5."."\n";
	exit;
}


my $orftype = "uORF";
### get arguments ###
my %args; 
getopt("fdsolrp",\%args);
my $readfile=$args{f};
if (! $readfile) {
	print "No corrected read mapping file"."\n";
    exit;
}
my $offFile=$args{s};
if (! $offFile) {
	print "No corrected read mapping file"."\n";
    exit;
}
my $orfdir=$args{d};
if (! $orfdir) {
	print "No candidate ORF annotation file"."\n";
    exit;
}
my $outputDir=$args{o};
if (! $outputDir) {
	print "No output directory"."\n";
    exit;
}

system  "mkdir -p $outputDir/temp" if ! -d "$outputDir/temp";

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
if (exists ($args{t})) {
	$threads=$args{t};
} else {
	$threads=3;
}


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
#open (IN, "$readfile");
open(IN, "samtools view -F 4 $readfile |");
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


my $pm = new Parallel::ForkManager($threads); 
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
print RI "pdf (file=\"$outputDir/plot.ROC.curve.$orftype.pdf\")"."\n";
print RI "plot(fpr, tpr, col=0, main=\"ROC Curve\", xlab=\"False positive rate\", ylab=\"True positive rate\")"."\n";
print RI "lines(fpr, tpr,col=2, lwd=3)"."\n";
print RI "dev.off()"."\n";
print RI "colnames(out) <- c(\"cutoff\", \"True.pos\", \"False.pos\", \"True.neg\", \"False.neg\", \"False.pos.rate\", \"True.pos.rate\")"."\n";
print RI "write.table (out, \"$outputDir/stat.cutoff.$orftype.txt\", sep=\"\\t\", quote=F, row.names=F)"."\n";
print RI "pred.pvalue <- sprintf(\"%.4f\", predict(pred, A))"."\n";
print RI "A <- as.data.table(A)"."\n";
print RI "A[, orfscore := mapply(function(t1,t2,t3,tot) log2(chisq.test(c(t1 * tot, t2 * tot, t3 * tot),p=c(1/3,1/3,1/3))\$statistic), f1, f2, f3, readNum)]"."\n";
print RI "A[f1 < f2 | f1 < f3, orfscore := orfscore * -1]"."\n";  
print RI "A <- cbind(A, pred.pvalue)"."\n";
print RI "strength <- fread(paste0('$orfdir','/uORF.merged.annotated.txt'),sep='\\t',head=F)\n"; 
print RI "strength <- strength[,.(orfID=V1,strength=V7)]\n"; 
print RI "A <- merge(A,strength,by='orfID',all.x=T)\n"; 
print RI "A[is.na(strength),strength:='']\n";
print RI "A[, c('strength','pred.pvalue') := .(pred.pvalue,strength)]\n"; 
print RI "setnames(A,c('strength','pred.pvalue'),c('pred.pvalue','strength'))\n";  
#print RI "A <- A[grep(\"uORF|extension|polycistronic\",orfID)]"."\n";
print RI "fwrite(A, file=\"$outputDir/pred.pvalue.parameters.$orftype.txt\", sep=\"\\t\")"."\n";
close RI;

my $com1="/public/home/liyiliang/anaconda3/envs/R4.1.3/bin/Rscript $outputDir/ribORF.learning.R";
system ($com1);


open (IN, "$outputDir/pred.pvalue.parameters.$orftype.txt");
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

#print Dumper(\%cluster);
my %sel;
foreach my $k (sort keys %cluster) {
	my @s1=split /~/, $cluster{$k};
	my $a1 = "";
	my $a2 = "";
	for (my $i=0; $i<=$#s1; $i++) {
		my @s3=split /:/, $s1[$i];
		if ($s3[1] eq 'ATG') {
			$a1.=$s1[$i]."|";
		} else {
			$a2.=$s1[$i]."|";
		} 
	}
	#print Dumper(\@s1);
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


my @forJson;
open (AN, "$outputDir/pred.pvalue.parameters.$orftype.txt");
open (JSON, ">$outputDir/pred.pvalue.parameters.filter.$orftype.txt") or die "cannot write";
open (NUM, ">$outputDir/pred.pvalue.parameters.lineNum.$orftype.txt") or die "cannot write";
open (OUT, ">$outputDir/repre.valid.pred.pvalue.parameters.$orftype.txt");
my $line=<AN>;
print OUT $line;
my $n = 0;
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /\|/, $s1[0];
	my @s3=split /:/, $s2[2];
	my $k=$s2[0].":".$s3[2].":".$s3[1].":".$s2[4];
	if($s1[$#s1] > $predictPvalueCutoff) {
		#my %record = ('ORF_ID' => $s1[0], 'Chr'=> $s1[1], 'Strand'=> $s1[2], 'CodonStartPos'=> $s1[3], 'CodonEndPos'=> $s1[4], 'OrfLength'=> $s1[5], 'ReadsNum'=> $s1[6], 'P-value'=> $s1[16]);
		#push @forJson, \%record;
		print JSON $_,"\n"; 
		$n++;
	}	
		
	if (exists ($sel{$k})) {
		print OUT $_."\n";
		#$n++;
	}
}
close AN;
close OUT;

#my $numberOfPage = int($n/20)+1;
#my %oder = ("last_page" => $numberOfPage, "data" => \@forJson);
#my $json = to_json(\%oder);
#print JSON "$json\n";
print NUM $n;
close NUM;
close JSON;

#system("rm -fr $outputDir/temp");
#system("rm -fr $outputDir/all.input.parameters.txt");

=pod
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
=cut


sub analysis_mul {
    my ($orfFile,$chr,$read) = @_;
	open (AN, "$orfFile");
	open (OUT, ">$outputDir/temp/$chr.input.parameters.txt");
	#print OUT "orfID"."\t"."chrom"."\t"."strand"."\t"."codon5"."\t"."codon3"."\t"."length"."\t"."readNum"."\t"."f1"."\t"."f2"."\t"."f3"."\t"."entropy"."\t"."MAXentropy"."\t"."PME"."\t"."codonNum"."\t"."f1max"."\n";
	my $dist=3;
	my $n=0;
	my $control = 0;
	while (<AN>) {
	chomp;
	$n++;
	#if($n % 1000 == 0){
	#   print $n,"\n";
	#}
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[8];
	my @s3=split /,/, $s1[9];
	my $len1=0;
	my $tot=0;
	my %val;
	my @per;
	my @post;
	my @fra = (0,0);
	for (my $i=0; $i<=2; $i++) {
		$per[$i]=0;
	}
	next if $s1[6] eq "";
	my $start = 0;
	my $end = $#s2;
	for (my $i=0; $i<=$#s2-1; $i++) {
	    if( $s1[5] >= $s2[$i] and $s1[5] < $s2[$i+1] ) {
		   $start = $i;
		   if($s1[5]> $s2[$i] and $s1[5] < $s3[$i]) {
		       $s2[$i] = $s1[5];
		   }
		}
		if( $s1[6] >= $s3[$i] and $s1[5] < $s3[$i+1] ) {
		   $end = $i+1;
		}
	}
	for (my $i=$start; $i<=$end; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			if ($j >=$s1[5] && $j < $s1[6]) {
				push @post, $j;
			}
			if($j >= $s1[6]){
			    last;
			}
		}
	}

	if ($s1[2] eq '-') {
		@post=reverse(@post);
	}
	if ($#post >= 5) {
	#splice @post, 0, 3;
	#splice @post, -3;
	#for (my $m=0; $m<=$#post; $m++) {
	#	$len1++;
	#	my $k=$s1[1].":".$s1[2].":".$post[$m]; 
	#	if (exists ($read{$k})) {
	#		$tot+=$read{$k};
	#		$val[int(($len1-1)/$dist)]+=$read{$k};
	#		$per[($len1-1)%3]+=$read{$k};
	#	}
	#}
	
		for (my $m=0; $m<=$#post-2; $m+=3) {
			my $k1=$s1[1].":".$s1[2].":".$post[$m];
			my $k2=$s1[1].":".$s1[2].":".$post[$m+1];
			my $k3=$s1[1].":".$s1[2].":".$post[$m+2];
			
			if (exists ($read{$k1})) {
				$tot+=$read{$k1};
				#if(not exists $val{int($m/$dist)}) {
				#   $val{int($m/$dist)}
				#}
				$val{int($m/$dist)}+=$read{$k1};
				$per[0]+=$read{$k1};
			}
			if (exists ($read{$k2})) {
				$tot+=$read{$k2};
				$val{int(($m+1)/$dist)}+=$read{$k2};
				$per[1]+=$read{$k2};
			}
			if (exists ($read{$k3})) {
				$tot+=$read{$k3};
				$val{int(($m+2)/$dist)}+=$read{$k3};
				$per[2]+=$read{$k3};
			}
			
			if(scalar keys %val == 0){
			  next;
			}
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
		
		if($fra[0] == 0 and $fra[1] = 0){
		   next;
		}

		my $len1 = scalar @post;
		if ($tot >= $orfReadCutoff && $len1 >= $orfLengthCutoff) {
			my $ent=0;
			my $ten=0;
			my $a=int(($len1+2)/$dist);
			my $b=$tot;
			my $t1=int(($a+$b-1)/$b);
			my %val2; 
			#print $a,"\t",$b,"\n";
			#for (my $i=0; $i<=$#val; $i++) {
			#	if (defined $val[$i] and $val[$i] > 0) {
			#		$val2[int($i/$t1)]+=$val[$i];
			#	}
			#}
			
			foreach my $i (sort keys %val){
					$val2{int($i/$t1)}+=$val{$i};
			}

			#for (my $i=0; $i<=$#val2; $i++) {
			#	if (defined $val2[$i] and $val2[$i] > 0) {
			#		my $p=$val2[$i]/($tot);
			#		$ent+=$p*log(1/$p);
			#	}
			#}
			foreach my $i (sort keys %val2){
					my $p=$val2{$i}/($tot);
					$ent+=$p*log(1/$p);
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
					my $p=$va[$i]/($tot);
					$ten+=$p*log(1/$p);
				}
			}
			my $per;
			if ($ten == 0) {
				$per=1;
			} else {
				$per=$ent/$ten;
			}

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
	close AN;
	close OUT;
}


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
