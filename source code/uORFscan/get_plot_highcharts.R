suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))
 
 
#main
args <- commandArgs(TRUE)
jobid <- args[1]

#jobid  <- "j8VIMjDT0eJRp8CM"

options(bitmapType='cairo')

color <- c("#7cb5ec", "#3a3a3c", "#90ed7d", "#f7a35c", "#8085e9", "#f15c80", "#e4d354", "#2b908f", "#f45b5b", "#91e8e1")

resultdir = paste0("./data/results/",jobid)

if(file.exists(paste0(resultdir,"/cleanStat.txt"))) {
	aa <- fread(paste0(resultdir,"/cleanStat.txt"),sep="\t")
	aa$Type <- factor(aa$Type, levels=c("rRNA","tRNA","snRNA","Clean"))
	aa[,Unique := Unique/sum(Unique)]
	aa[,Total := Total/sum(Total)]
	plotlist <- list()
	lab <- paste0(aa$Type," (",round(aa$Unique*100,2),"%",")") 
	p1 <- ggdonutchart(aa, "Unique", label = rep("",length(aa$Type)), fill = "Type", lab.pos="in", color = "white", palette = c('#FF6347','#F4A460','#FFA500','#00FF7F'),legend.title="")
	plotlist[['unique']] <- p1
	lab <- paste0(aa$Type," (",round(aa$Total*100,2),"%",")") 
	p2 <- ggdonutchart(aa, "Total", label = rep("",length(aa$Type)), fill = "Type", lab.pos="in", color = "white", palette = c('#FF6347','#F4A460','#FFA500','#00FF7F'),legend.title="")
	plotlist[['total']] <- p2
	glist <- lapply(plotlist, ggplotGrob)
	ggsave(paste0(resultdir, "/", file="cleanStat.pdf"), width = 8, height =5, marrangeGrob(grobs = glist, nrow=1, ncol=2, top=NULL))
	ggsave(paste0(resultdir, "/", file="cleanStat.pdf"), width = 8, height =5, marrangeGrob(grobs = glist, nrow=1, ncol=2, top=NULL))
}

if(file.exists(paste0(resultdir,"/geneBiotypeStat.txt"))) {
	aa <- fread(paste0(resultdir,"/geneBiotypeStat.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Biotype))
	ymax <- max(aa$Frequency)+(max(aa$Frequency)/5)
	ggbarplot(aa, "Biotype", "Frequency", fill = "Biotype", color = "Biotype",palette = palettes, legend = "none", orientation = "horiz",label = TRUE, label.pos = "out", lab.vjust=0.5, lab.hjust=-0.1, xlab = "", ylab = "") + ylim(0,ymax) 
	ggsave(paste0(resultdir,"/geneBiotypeStat.pdf"), width = 8, height = 8)
	ggsave(paste0(resultdir,"/geneBiotypeStat.png"), width = 8, height = 8)
}

if(file.exists(paste0(resultdir,"/geneFeatureStat.txt"))) {
	aa <- fread(paste0(resultdir,"/geneFeatureStat.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Feature))
	ymax <- max(aa$Frequency)+(max(aa$Frequency)/10)
	ggbarplot(aa, "Feature", "Frequency", fill = "Feature", color = "Feature",palette = palettes, legend = "none", label = TRUE, label.pos = "out", xlab = "", ylab = "RPF frequency") + ylim(0,ymax)
	ggsave(paste0(resultdir,"/geneFeatureStat.pdf"), width = 8, height = 6)
	ggsave(paste0(resultdir,"/geneFeatureStat.png"), width = 8, height = 6)
}

if(file.exists(paste0(resultdir,"/rawLengthStat.txt"))) {
	aa <- fread(paste0(resultdir,"/rawLengthStat.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Length))
	aa$Length <- as.factor(aa$Length)
	aa$Percentage <- aa$Percentage*100
	ymax <- max(aa$Percentage)+max(aa$Percentage)/10
	label_text <- paste0(round(aa$Percentage,1),"%")
	ggbarplot(aa, "Length", "Percentage", fill = "Length", color = "Length",palette = palettes, legend = "none", label = label_text, label.pos = "out", xlab = "", ylab = "Percentage") + ylim(0,ymax)
	ggsave(paste0(resultdir,"/rawLengthStat.pdf"), width = 8, height = 6)
	ggsave(paste0(resultdir,"/rawLengthStat.png"), width = 8, height = 6)
}

if(file.exists(paste0(resultdir,"/RPFlengthStat.txt"))) {
	aa <- fread(paste0(resultdir,"/RPFlengthStat.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Percentage))
	ymax <- max(aa$Percentage)+(max(aa$Percentage)/10)
	aa$Length <- as.factor(aa$Length)
	ymax <- max(aa$Percentage)+max(aa$Percentage)/10
	label_text <- paste0(round(aa$Percentage,2),"%")
	ggbarplot(aa, "Length", "Percentage", fill = "Length", color = "Length",palette = palettes, legend = "none", label = label_text, label.pos = "out", xlab = "", ylab = "Percentage") + ylim(0,ymax)
	ggsave(paste0(resultdir,"/RPFlengthStat.pdf"), width = 8, height = 6)
	ggsave(paste0(resultdir,"/RPFlengthStat.png"), width = 8, height = 6)
}

if(file.exists(paste0(resultdir,"/ORFstat.txt"))) {
	aa <- fread(paste0(resultdir,"/ORFstat.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Count))
	ymax <- max(aa$Count)+(max(aa$Count)/10)
	ggbarplot(aa, "Type", "Count", fill = "Type", color = "Type",palette = palettes, legend = "none", label = TRUE, label.pos = "out", xlab = "", ylab = "ORF frequency") + ylim(0,ymax) + rotate_x_text()
	ggsave(paste0(resultdir,"/ORFstat.pdf"), width = 8, height = 6)
	ggsave(paste0(resultdir,"/ORFstat.png"), width = 8, height = 6)
}

if(file.exists(paste0(resultdir,"/ORFtypeStat.txt"))) {
	aa <- fread(paste0(resultdir,"/ORFtypeStat.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Count))
	ymax <- max(aa$Count)+(max(aa$Count)/10)
	ggbarplot(aa, "Type", "Count", fill = "Type", color = "Type",palette = palettes, legend = "none", label = TRUE, label.pos = "out", xlab = "", ylab = "ORF frequency") + ylim(0,ymax) + rotate_x_text()
	ggsave(paste0(resultdir,"/ORFtypeStat.pdf"), width = 8, height = 6)
	ggsave(paste0(resultdir,"/ORFtypeStat.png"), width = 8, height = 6)
}

if(file.exists(paste0(resultdir,"/frameStatByLength.txt"))) {
	aa <- fread(paste0(resultdir,"/frameStatByLength.txt"),sep="\t")
	palettes <- rep(color, length.out = length(aa$Frequency))
	ymax <- max(aa$Frequency)+(max(aa$Frequency)/10)
	aa$Frame <- as.factor(aa$Frame)
	ggbarplot(aa, "Length", "Frequency", fill = "Frame", color = "Frame",palette = palettes, legend = "bottom", xlab = "", ylab = "RPF frequency", position = position_dodge(0.9)) + ylim(0,ymax) + scale_x_continuous(name="", breaks=seq(min(aa$Length),max(aa$Length),1), labels = seq(min(aa$Length),max(aa$Length),1))
	ggsave(paste0(resultdir,"/frameStatByLength.pdf"), width = 12, height = 6)
	ggsave(paste0(resultdir,"/frameStatByLength.png"), width = 8, height = 6)
}

resultTxt <- dir(resultdir, full.names = TRUE, pattern = ".txt$")
resultTxt <- dir(resultdir, full.names = TRUE, pattern = ".csv$")
resultPng <- dir(resultdir, full.names = TRUE, pattern = ".png$")
resultPdf <- dir(resultdir, full.names = TRUE, pattern = ".pdf$")
resultdir2zip <- c(resultTxt, resultPdf, resultPng)
excludeFile <- c("activeORF_collapsed.txt","featurecount.unique.txt","featurecount.exp.txt","ribocode.config.txt","activeORF_ORFs_category.pdf")
excludeFile <- paste0(resultdir,"/",excludeFile)
resultdir2zip <- setdiff(resultdir2zip,excludeFile)
zip(zipfile = paste0(resultdir,"/",jobid,".results.zip"), files = resultdir2zip, flags="-q")