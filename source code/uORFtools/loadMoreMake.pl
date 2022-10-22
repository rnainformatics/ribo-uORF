#!/usr/bin/perl -w
#make 'load more' and 'pdf links' html elements
use strict;

my $jobid = $ARGV[0];

if (! -d "./data/results/$jobid/html") {
  system  "mkdir -p ./data/results/$jobid/html" ;
}



my @qctypes = ("length", "frame", "period");
foreach my $qctype (@qctypes) {
	my $picfile = "./data/results/$jobid/qc_$qctype.png";
	my $pdffile = "./data/results/$jobid/qc_$qctype.pdf";
	my $pdfHtmlfile = "./data/results/$jobid/html/qc_$qctype.html";
	pdfMaker($pdffile,$pdfHtmlfile);
}


my @loci = ("E", "P", "A", "+1", "+2", "+3");
my %loci_raname = ('E'=>'E', 'P'=>'P', 'A'=>'A', '+1'=>'+1 (Relative to A site)', '+2'=>'+2 (Relative to A site)', '+3'=>'+3 (Relative to A site)');
open OUT, ">./data/results/$jobid/html/loadmore_EAP.html" || die "$!";
print OUT '<div class="cbp-loadMore-block1">\n';
for(my $i=0; $i<@loci; $i++) {
my $locus = $loci[$i];
my $picfile = "./data/results/$jobid/usage_graph_position_$locus.png";
my $pdffile = "./data/results/$jobid/usage_graph_position_$locus.pdf";
my $pdfHtmlfile = "./data/results/$jobid/html/usage_graph_position_$locus.html";
pdfMaker($pdffile,$pdfHtmlfile);
next if($i<3);
my $html_ele = <<DATA;
<div class="cbp-item $locus logos">
					<div class="cbp-caption">
						<div class="cbp-caption-defaultWrap">
							<img src="$picfile" alt=""> </div>
						<div class="cbp-caption-activeWrap">
							<div class="cbp-l-caption-alignCenter">
								<div class="cbp-l-caption-body">
									<a href="$pdfHtmlfile" class="cbp-singlePage cbp-l-caption-buttonLeft btn red uppercase btn red uppercase embed-link" rel="nofollow">View PDF</a>
									<a href="$picfile" class="cbp-lightbox cbp-l-caption-buttonRight btn red uppercase btn red uppercase" data-title="World Clock Widget<br>by Paul Flavius Nechita">view larger</a>
								</div>
							</div>
						</div>
					</div>
					<div class="cbp-l-grid-projects-title text-center">$loci_raname{$locus} site</div>
				</div>
DATA
print OUT $html_ele."\n";
}
print OUT '</div>\n';
close OUT;



my @aas = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");
open OUT, ">./data/results/$jobid/html/loadmore_AA.html" || die "$!";
print OUT '<div class="cbp-loadMore-block1">\n';
for(my $i=0; $i<@aas; $i++) {
my $aa = $aas[$i];
my $picfile = "./data/results/$jobid/aa_graph.$aa.png";
my $pdffile = "./data/results/$jobid/aa_graph.$aa.pdf";
my $pdfHtmlfile =  "./data/results/$jobid/html/aa_graph.$aa.html";
pdfMaker($pdffile,$pdfHtmlfile);
next if($i<8);
my $html_ele = <<DATA;
<div class="cbp-item $aa logos">
					<div class="cbp-caption">
						<div class="cbp-caption-defaultWrap">
							<img src="$picfile" alt=""> </div>
						<div class="cbp-caption-activeWrap">
							<div class="cbp-l-caption-alignCenter">
								<div class="cbp-l-caption-body">
									<a href="$pdfHtmlfile" class="cbp-singlePage cbp-l-caption-buttonLeft btn red uppercase btn red uppercase" rel="nofollow">View PDF</a>
									<a href="$picfile" class="cbp-lightbox cbp-l-caption-buttonRight btn red uppercase btn red uppercase" data-title="World Clock Widget<br>by Paul Flavius Nechita">view larger</a>
								</div>
							</div>
						</div>
					</div>
					<div class="cbp-l-grid-projects-title uppercase text-center">$aa</div>
				</div>
DATA
if(-e $picfile) {
	print OUT $html_ele."\n";
}
}
print OUT '</div>\n';
close OUT;


my @codons = ("AAA", "AAC", "AAG", "AAU", "ACA", "ACC", "ACG", "ACU", "AGA", "AGC", "AGG", "AGU", "AUA", "AUC", "AUG", "AUU", "CAA", "CAC", "CAG", "CAU", "CCA", "CCC", "CCG", "CCU", "CGA", "CGC", "CGG", "CGU", "CUA", "CUC", "CUG", "CUU", "GAA", "GAC", "GAG", "GAU", "GCA", "GCC", "GCG", "GCU", "GGA", "GGC", "GGG", "GGU", "GUA", "GUC", "GUG", "GUU", "UAC", "UAU", "UCA", "UCC", "UCG", "UCU", "UGC", "UGG", "UGU", "UUA", "UUC", "UUG", "UUU");
open OUT, ">./data/results/$jobid/html/loadmore_CODON1.html" || die "$!";
print OUT '<div class="cbp-loadMore-block1">\n';
for(my $i=0; $i<@codons; $i++) {
my $codon = $codons[$i];
my $picfile = "./data/results/$jobid/codon_graph.$codon.png";
my $pdffile = "./data/results/$jobid/codon_graph.$codon.pdf";
my $pdfHtmlfile = "./data/results/$jobid/html/codon_graph.$codon.html";
pdfMaker($pdffile,$pdfHtmlfile);
next if($i<12);
my $html_ele = <<DATA;
<div class="cbp-item $codon logos">
					<div class="cbp-caption">
						<div class="cbp-caption-defaultWrap">
							<img src="$picfile" alt=""> </div>
						<div class="cbp-caption-activeWrap">
							<div class="cbp-l-caption-alignCenter">
								<div class="cbp-l-caption-body">
									<a href="$pdfHtmlfile" class="cbp-singlePage cbp-l-caption-buttonLeft btn red uppercase btn red uppercase" rel="nofollow">View PDF</a>
									<a href="$picfile" class="cbp-lightbox cbp-l-caption-buttonRight btn red uppercase btn red uppercase" data-title="World Clock Widget<br>by Paul Flavius Nechita">view larger</a>
								</div>
							</div>
						</div>
					</div>
					<div class="cbp-l-grid-projects-title uppercase text-center">$codon</div>
				</div>
DATA
if(-e $picfile) {
	print OUT $html_ele."\n";
}
}
print OUT '</div>\n';
close OUT;


@codons = ("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT");
open OUT, ">./data/results/$jobid/html/loadmore_CODON2.html" || die "$!";
print OUT '<div class="cbp-loadMore-block1">\n';
for(my $i=0; $i<@codons; $i++) {
my $codon = $codons[$i];
my $picfile = "./data/results/$jobid/metacodon_graph.$codon.png";
my $pdffile = "./data/results/$jobid/metacodon_graph.$codon.pdf";
my $pdfHtmlfile = "./data/results/$jobid/html/metacodon_graph.$codon.html";
pdfMaker($pdffile,$pdfHtmlfile);
next if($i<12);
my $html_ele = <<DATA;
<div class="cbp-item $codon logos">
					<div class="cbp-caption">
						<div class="cbp-caption-defaultWrap">
							<img src="$picfile" alt=""> </div>
						<div class="cbp-caption-activeWrap">
							<div class="cbp-l-caption-alignCenter">
								<div class="cbp-l-caption-body">
									<a href="$pdfHtmlfile" class="cbp-singlePage cbp-l-caption-buttonLeft btn red uppercase btn red uppercase" rel="nofollow">View PDF</a>
									<a href="$picfile" class="cbp-lightbox cbp-l-caption-buttonRight btn red uppercase btn red uppercase" data-title="World Clock Widget<br>by Paul Flavius Nechita">view larger</a>
								</div>
							</div>
						</div>
					</div>
					<div class="cbp-l-grid-projects-title uppercase text-center">$codon</div>
				</div>
DATA
if(-e $picfile) {
	print OUT $html_ele."\n";
}
}
print OUT '</div>\n';
close OUT;
 
 
sub pdfMaker {
	my $infilename = shift;
	my $outfilename = shift;
	open HTML, ">$outfilename" || die "$!";
	my $pdfhtml = <<DATA;
<style>
.pdfobject-container {
		width: 100%;
		max-width: 1000px;
		height: 1000px;
		margin: 2em 0;
}
.pdfobject { border: solid 1px #666; }
</style>
<script src="./assets/global/plugins/pdfobject.min.js"></script>
<div id="pdf" class="pdfobject-container"></div>
<script>
var options = {
	pdfOpenParams: {
		pagemode: "thumbs",
		navpanes: 0,
		toolbar: 0,
		statusbar: 0,
		view: "FitV"
	}
};
var myPDF = PDFObject.embed("$infilename", "#pdf", options);
</script>
DATA
    print HTML $pdfhtml;
	close HTML;
}


#'<div class="cbp-loadMore-block1"><div class="cbp-item web-design logos"><div class="cbp-caption"><div class="cbp-caption-defaultWrap"><img src="./data/results/$jobid/pic"> </div><div class="cbp-caption-activeWrap"><div class="cbp-l-caption-alignCenter"><div class="cbp-l-caption-body"><a href="../assets/global/plugins/cubeportfolio/ajax/project1.html" class="cbp-singlePage cbp-l-caption-buttonLeft btn red uppercase" rel="nofollow">view PDF</a><a href="./data/results/$jobid/pic" class="cbp-lightbox cbp-l-caption-buttonRight btn red uppercase" data-title="Shopping Gallery<br>by Cosmin Capitanu">view larger</a></div></div></div></div><div class="cbp-l-grid-projects-title uppercase text-center">Shopping Gallery</div></div></div>'
#'<div class="cbp-loadMore-block2"><div class="cbp-item web-design logos"><div class="cbp-caption"><div class="cbp-caption-defaultWrap"><img src="./data/results/$jobid/pic" alt=""> </div><div class="cbp-caption-activeWrap"><div class="cbp-l-caption-alignCenter"><div class="cbp-l-caption-body"><a href="../assets/global/plugins/cubeportfolio/ajax/project1.html" class="cbp-singlePage cbp-l-caption-buttonLeft btn red uppercase" rel="nofollow">view PDF</a><a href="./data/results/$jobid/pic" class="cbp-lightbox cbp-l-caption-buttonRight btn red uppercase" data-title="Mountaineer<br>by Cosmin Capitanu">view larger</a></div></div></div></div><div class="cbp-l-grid-projects-title uppercase text-center">Mountaineer</div></div></div>'