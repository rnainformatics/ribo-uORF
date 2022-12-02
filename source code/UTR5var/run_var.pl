#!/usr/bin/perl -w

my ($ref, $jobid)=@ARGV;

my $samplenames = "NA";

my $dbbasepath =   "./db/annotation/";
my $utr5bed = "$dbbasepath/$ref.utr5.bed.nonchr.bed"; 

#$minlen ||= 25;

my $assembly;
my $species;
if($ref eq "hg19") {
	$assembly = "GRCh37";
	$species = "homo_sapiens";
} elsif($ref eq "hg38"){
	$assembly = "GRCh38";
	$species = "homo_sapiens";
} elsif($ref eq "mm10") {
	$assembly = "GRCm38";
	$species = "mus_musculus";
}

$resultdir = "./data/results/$jobid";

if($samplenames eq "NA") {
   $samplenames = $jobid;
}

&get_log_status("Data importing","Importing",10);
my $rec = 0;
open VCF, "./data/results/$jobid/$jobid.vcf" || die "$!";
while(<VCF>) {
   chomp;
   next if (/^#/);
   $rec++;
}
close VCF;
if($rec == 0){
   &get_log_status("No variation locating on the 5UTR region","Error",20);
   exit;
}
system("sed 's/^chr//i' ./data/results/$jobid/$jobid.vcf > ./data/results/$jobid/$jobid.nochr.vcf");

&get_log_status("Filtering the variation on 5UTR","Filtering",20);
system("vcftools --vcf ./data/results/$jobid/$jobid.nochr.vcf --bed $utr5bed --out ./data/results/$jobid/utr5 --recode --keep-INFO-all");

#system("intersectBed -a ./data/results/$jobid/$jobid.nochr.vcf -b $utr5bed -header > ./data/results/$jobid/utr5.recode.vcf");

##filter the rRNA defived reads
&get_log_status("Annotate the variations using VEP","Annotating",50);
system("vep  --cache --dir_cache /public/liuqi/ensembl-vep/vep_cache --format vcf --fork 6 -i ./data/results/$jobid/utr5.recode.vcf --species $species --assembly $assembly --plugin UTRannotator,/public/liuqi/ensembl-vep/plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt --dir_plugins /public/liuqi/ensembl-vep/plugins/UTRannotator --tab --host ensembldb.ensembl.org --port 5306 --force_overwrite -o ./data/results/$jobid/vep.out");

#if($ref eq "hg38" || $ref eq "mm10"){
#	system("Rscript ./program/utr.annotation.R $species $jobid");
#}
system("Rscript ./program/parseVep.R $jobid $ref");
my $nrecord = 0;
open SEQ, "$resultdir/vep.effected.txt" || die "$!";
<SEQ>;
while(<SEQ>) {
   chomp; 
   $nrecord++;
}
close SEQ;
if($nrecord < 2){
   &get_log_status("No effect detected in your data","Error",80);
   exit;
}

system("cat ./data/results/$jobid/vep.out_summary.html | sed 's/\\/mnt\\/data\\/software\\/ensembl-vep\\/cache\\///g' | sed 's/.\\/data\\/results\\///g' > ./data/results/$jobid/vep.out_summary_rev.html");       
&get_log_status("Job completed. Thanks for using Ribo-uORF!","Job complete",100);


sub get_log_status()
{
    my $info = $_[0];
	my $status = $_[1];
	my $percent = $_[2];
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon +=1;
	my $logdate=sprintf("%04d-%02d-%02d %02d:%02d:%02d\n",$year,$mon,$mday,$hour,$min,$sec);
	$logdate = "$info\t$logdate";
	print($logdate);
	system("echo '{\"type\":\"varAnnote\",\"species\":\"$ref\",\"sample\":\"$samplenames\",\"current\":\"$status\",\"percent\":\"$percent\",\"message\":\"$info\"}' > $resultdir/status.json");
	system("echo '$logdate' >> $resultdir/log.txt");
}

