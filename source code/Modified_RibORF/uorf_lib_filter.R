suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(data.table))


args <- commandArgs(TRUE)
org <- args[1]


txdb <- loadDb(paste0(org,".gencode.sqlite"))
cdsRange <- cdsBy(txdb,use.names=T)
cds5site <- as.data.table(cdsRange)

tr2gene <- fread(paste0(org,".txlens.txt"),head=T,sep="\t")
tr2gene <- tr2gene[,.(tx_name,gene_id)]

peakfile = paste0(org,"/all.uORF.utr5.bed")
peak_site = fread(peakfile,head=F,sep="\t")
setnames(peak_site,  c("chr","ustart","uend","ID","strand","numCDS"))

peak_site2 <- peak_site[,  tstrsplit(ID, "|", fixed=TRUE)]
peak_site3 <- peak_site2[,  tstrsplit(V3, ":", fixed=TRUE)]
peak_site5 <- peak_site2[,  tstrsplit(V1, ":", fixed=TRUE)]

peak_site_gene <- cbind(peak_site, peak_site5[,.(tx_name=V1)])
peak_site_gene <- cbind(peak_site_gene, peak_site3[,.(trStart=V2,trEnd=V3)])
peak_site_gene <- merge(peak_site_gene,tr2gene,by="tx_name",all.x=T)

peak_site_gene[,uORFid := paste0(gene_id,"_",chr,"_",ustart,"_",uend)]
setorder(peak_site_gene,uORFid)
peak_site_gene <- peak_site_gene[,.(chr,ustart,uend,ID,strand,numCDS,tx_name, trStart, trEnd, gene_id,uORFid)]
fwrite(peak_site_gene,file=paste0(org,"/all.uORF.utr5.gene.bed"),sep="\t")
peak_site_gene_genomic_bed <- peak_site_gene[,.(chr, ustart, uend, uORFid, numCDS, strand)]
peak_site_gene_genomic_bed <- unique(peak_site_gene_genomic_bed,by="uORFid")
fwrite(peak_site_gene_genomic_bed,file=paste0(org,"/all.uORF.utr5.gene.genomic.bed6"),sep="\t",col.names=F)
peak_site_gene_transcript_bed <- peak_site_gene[,.(tx_name, trStart, trEnd, uORFid, "*", "*")]
peak_site_gene_transcript_bed <- unique(peak_site_gene_transcript_bed,by="uORFid")
fwrite(peak_site_gene_transcript_bed,file=paste0(org,"/all.uORF.utr5.gene.transcript.bed6"),sep="\t",col.names=F)

peak_site <- cbind(peak_site, peak_site3[,.(V2,V3)])
peak_site[,V2:=as.numeric(V2)]
peak_site[,V3:=as.numeric(V3)]
peak_site[,orfLen:=V3-V2]
peak_site_short <- peak_site[orfLen < 18]
peak_site_short[,type:="too_short"]
peak_site <- peak_site[orfLen >= 18]

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
peak_site <- peak_site[queryHits %in% varclashm$queryHits]
peak_site[,c("queryHits") := NULL]
peak_site[,type:="overlap_CDS"]

peak_site <- rbind(peak_site_short, peak_site)
fwrite(peak_site,file=paste0(org,"/all.uORF.filterd.bed"),col.names=F,sep="\t",quote=F)



