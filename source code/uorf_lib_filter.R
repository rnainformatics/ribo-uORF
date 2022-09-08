suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(data.table))


args <- commandArgs(TRUE)
org <- args[1]



txdb <- loadDb(paste0("/public/home/liyiliang/ribo/ribotoolkit/db/annotation/",org,".gencode.sqlite"))
cdsRange <- cdsBy(txdb,use.names=T)
cds5site <- as.data.table(cdsRange)

peakfile = paste0("all.uORF.utr5.bed")
peak_site = fread(peakfile,head=F,sep="\t")
setnames(peak_site,  c("chr","ustart","uend","ID","strand","numCDS"))

peak_range = GRanges(seqnames=peak_site[,chr], ranges=IRanges(peak_site[,ustart], peak_site[,uend]),strand=peak_site[,strand])
cn <- colnames(peak_site)
if (length(cn) > 3) {
	for (j in c(4)) {
		mcols(peak_range)[[cn[j]]] <- peak_site[, eval(as.name(cn[j]))]
	}
}
overlaps = findOverlaps(peak_range, cdsRange, minoverlap=3)
varclashm = as.data.table(overlaps)
peak_site[, queryHits := .I]
peak_site <- peak_site[!queryHits %in% varclashm$queryHits]
peak_site[,c("queryHits") := NULL]
	

if(org == "hg38") {
	obj_chr = c("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")
} else if(org == "mm10") {
	obj_chr = c("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")
} else if(org == "rn6") {
	obj_chr = c("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")
} else if(org == "dre_GRCz11") {
	obj_chr = c("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr1","chr20","chr21","chr22","chr23","chr24","chr25","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM")
} else if(org == "dme_BDGP6") {
	obj_chr = c("chr2L","chr2R","chr3L","chr3R","chr4","chrM","chrX","chrY");
} else if(org == "cel_WBcel235") {
	obj_chr = c("chrI","chrII","chrIII","chrIV","chrM","chrV","chrX");
}

libfiles <- list.files("./",pattern="*.candidateORF.genepred.txt$")
for(peakfile in libfiles) {
	#peakfile = paste0(i,".candidateORF.genepred.txt")
	uorf_site = fread(peakfile,head=F,sep="\t")
	setnames(uorf_site, c("ID","chr","strand","txstart","txend","cdsStart","cdsEnd","numCDS","blockStart","blockEnd"))
	uorf_site <- uorf_site[ID %in% peak_site$ID]
	fwrite(uorf_site,file=paste0(peakfile,".uORF"),col.names=F,sep="\t",quote=F)
}

