suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicFeatures))
#suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
#library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library("zoo"))
suppressPackageStartupMessages(library("signal"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library(Rsamtools))

set_offset_variable <- function(reads, offset_list) {
	seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
	reads[, read_number := .I]
	reads[, strand := as.character(strand)]
	reads[, transcript_coordinate := as.integer(0)]
	
	# For each read length, apply the corresponding offset to those reads
	for (i in names(offset_list)) {
		len <- as.integer(i)
		offset <- as.integer(offset_list[[i]])
		reads[strand == '+' & match_len == len, transcript_coordinate:=as.integer(start+offset)]
		reads[strand == '-' & match_len == len, transcript_coordinate:=as.integer(end-offset)]
	}
	reads <- reads[transcript_coordinate != 0]
	return(reads)
}

write.wig <- function(coverage_table, prefix){
	options(scipen=999) # Avoid scientific notation
	# Set visualization options
	track_options <- " type=bedGraph color=0,0,255 autoScale=on maxHeightPixels=250:200:50 windowingFunction=none"
	# Write positive strand
	file_fw <- paste0(prefix, '_fw.wig')
	track_line_fw <- paste0("track name=", file_fw, track_options)
	cat(track_line_fw, '\n', file = file_fw)
	write.table(coverage_table[strand == '+', list(chrom, transcript_coordinate-1, transcript_coordinate, coverage)],
	append = T, col.names = F, file = file_fw, row.names = F, sep = '\t', quote= F)
	# Write negative strand
	file_rc <- paste0(prefix, '_rc.wig')
	track_line_rc <- paste0("track name=", file_rc, track_options)
	cat(track_line_rc, '\n', file = paste0(prefix, '_rc.wig'))
	write.table(coverage_table[strand == '-', list(chrom, transcript_coordinate-1, transcript_coordinate, coverage)], append = T, col.names = F, file = file_rc, row.names = F, sep = '\t', quote = F) 
	options(scipen=0) # Reset to default
}

plotlength <- function (reads, cl = 100)
{
    seqlength=reads$match_len
    dist <- table(factor(seqlength, levels = c(min(seqlength):max(seqlength))))
    dist <- data.frame(Length = names(dist), Percentage = as.vector((dist/sum(dist)) * 100))
    xmin <- quantile(seqlength, (1 - cl/100)/2)
    xmax <- quantile(seqlength, 1 - (1 - cl/100)/2)
    p <- ggplot(dist, aes(as.numeric(as.character(Length)), Percentage)) + 
       geom_bar(stat = "identity", fill = "gray80") + labs(title = "", 
     x = "Read length", y = "Count (%)") + theme_bw(base_size = 18) +
        theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits = c(xmin - 
        0.5, xmax + 0.5), breaks = seq(xmin + ((xmin)%%2), xmax, by = floor((xmax - xmin)/(xmax - xmin))))
	fwrite(dist, file=paste0(resultdir, "/", file="RPFlengthStat.txt"), sep="\t")
	fwrite(dist, file=paste0(resultdir, "/", file="RPFlengthStat.csv"))
    ret_list <- list()
    #ret_list[["plot"]] <- p
    return(p)
}


frame_psite_length = function (df, samplename="sample", region = "all", cl = 100, length_range = NULL) 
{
    if (!region %in% c("all", "cds", "5utr", "3utr")) {
        warning("region is invalid. Set to default \"all\"\n")
        region = "all"
    }
	if (region == "all") {
		#df <- subset(df, start_pos != 0 & stop_pos != 0)
		df = df[start_pos != 0 & stop_pos != 0,]
		df[,frame:=site_from_start%%3]
		minl <- quantile(df$match_len, (1 - cl/100)/2)
		maxl <- quantile(df$match_len, 1 - (1 - cl/100)/2)
		if (length(length_range) != 0) {
			if (!inherits(length_range, "numeric") & !inherits(length_range, 
			  "integer")) {
			  warning("length_range is invalid. Confidence interval is used\n")
			}
			else {
			  minl <- min(length_range)
			  maxl <- max(length_range)
			}
		}
		if (exists("final_frame_df")) {
           rm(final_frame_df)
	    }
		lenmax <- length(minl:maxl)
		#frame_df <- as.data.table(table(factor(df$match_len, 
        #        levels = minl:maxl), df$frame, factor(df$psite_region, 
        #        levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", 
        #          "CDS", "3' UTR"))))
		frame_df = df[,sum(freq),by=.(match_len,frame, psite_region)] 
		frame_df$match_len = factor(frame_df$match_len,levels = minl:maxl)	
		frame_df$psite_region = factor(frame_df$psite_region,levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", "CDS", "3' UTR"))		
        setnames(frame_df, c("length", "frame", "region", "count"))
		sum_count = frame_df[,list(sum(count)),by=c("length","region")]
        #sum_count <- as.vector(by(frame_df$count, frame_df[, c("length", "region")], function(x) sum(x)))
		sum_count = sum_count$V1
        frame_df[,perc := (count/c(rep(sum_count[1:lenmax], 
                3), rep(sum_count[(lenmax + 1):(2 * lenmax)], 
                3), rep(sum_count[(2 * lenmax + 1):(3 * lenmax)], 
                3))) * 100]
	} else {		
		df = df[site_region == region,]
		df[,frame:=site_from_start%%3]
		minl <- quantile(df$match_len, (1 - cl/100)/2)
		maxl <- quantile(df$match_len, 1 - (1 - cl/100)/2)
		if (length(length_range) != 0) {
			if (!inherits(length_range, "numeric") & !inherits(length_range, 
			  "integer")) {
			  warning("length_range is invalid. Confidence interval is used\n")
			}
			else {
			  minl <- min(length_range)
			  maxl <- max(length_range)
			}
		}
		lenmax <- length(minl:maxl)
		#frame_df <- as.data.table(table(factor(df$match_len, 
		#	levels = minl:maxl), df$frame))
		frame_df = df[,sum(freq),by=.(match_len,frame)] 	
		frame_df$match_len = factor(match_len,levels = minl:maxl)	
		setnames(frame_df, c("length", "frame", "count"))
		sum_count = frame_df[,list(sum(count)),by=c("length")]
		#sum_count <- as.vector(by(frame_df$count, frame_df$length, function(x) sum(x)))
		sum_count = sum_count$V1
		frame_df$perc <- (frame_df$count/rep(sum_count, 3)) * 100
	}
	
	frame_df[,samp := samplename]
	if (exists("final_frame_df")) {
		final_frame_df <- rbind(final_frame_df, frame_df)
	} else {
		final_frame_df <- frame_df
	}
	rm(frame_df)

    final_frame_df[is.na(perc), perc := 0]
	final_frame_df$region =factor(final_frame_df$region,levels = c("5' UTR", "CDS", "3' UTR"))
    plot <- ggplot(final_frame_df, aes(frame, as.numeric(as.character(length)))) + 
        geom_tile(aes(fill = perc)) + scale_fill_gradient("Normalized read number ", 
        low = "white", high = "#104ec1", breaks = c(min(final_frame_df$perc), 
            (max(final_frame_df$perc) - min(final_frame_df$perc))/2, 
            max(final_frame_df$perc)), labels = c(round(min(final_frame_df$perc)), 
            round((max(final_frame_df$perc) - min(final_frame_df$perc))/2), 
            round(max(final_frame_df$perc)))) + labs(x = "Frame", 
        y = "Read length") + theme_bw(base_size = 20) + theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()) + theme(legend.position = "top") + 
        scale_y_continuous(limits = c(minl - 0.5, maxl + 0.5), 
            breaks = seq(minl + ((minl)%%2), maxl, by = max(2, 
                floor((maxl - minl)/7))))
    if (region == "all") {
        plot <- plot + facet_grid(samp ~ region)
    } else {
        plot <- plot + facet_wrap(~samp, ncol = 3)
    }
    ret_list <- list()
    #ret_list[["df"]] <- final_frame_df
    #ret_list[["plot"]] <- plot
    return(plot)
}


psite <- function(data, flanking = 6, start = TRUE, extremity = "auto", cl = 99) {
    names <- names(data)
    offset <- NULL
    n <- "sample"
    dt <- data
    lev <- sort(unique(dt$qwidth))
    
    if(start == T | start == TRUE){
      base <- 0
      dt[, site_dist_end5 := start - start_pos]
      dt[, site_dist_end3 := end - start_pos]
    } else {
      base <- -5
      dt[, site_dist_end5 := start - stop_pos - base]
      dt[, site_dist_end3 := end - stop_pos - base]
    }
	#if(species == "sce_R64" | species == "ecoli_k12" | species == "bsu_168" | species == "pfu_dsm_3638" | species == "hsa_NRC1") {
	#	site_sub <- dt[site_dist_end5 <= 2*flanking]
	#} else {
	#	site_sub <- dt[site_dist_end5 <= -flanking & site_dist_end3 >= flanking - 1]
	#}
	site_sub <- dt
    minlen <- min(site_sub$qwidth)
    maxlen <- max(site_sub$qwidth)
    t <- table(factor(site_sub$qwidth, levels = lev))
    
    # offset
    offset_temp <- data.table(qwidth = as.numeric(as.character(names(t))), percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    offset_temp[, around_site := "T"][percentage == 0, around_site := "F"]
    offset_temp5 <- site_sub[, list(offset_from_5 = as.numeric(names(which.max(table(site_dist_end5))))), by = qwidth]
    offset_temp3 <- site_sub[, list(offset_from_3 = as.numeric(names(which.max(table(site_dist_end3))))), by = qwidth]
    merge_allx <- function(x, y) merge(x, y, all.x = TRUE, by = "qwidth")
    offset_temp  <-  Reduce(merge_allx, list(offset_temp, offset_temp5, offset_temp3))
    
    # adjusted offset
    adj_off <- function(dt_site, dist_site, add, bestoff){
      temp_v <- dt_site[[dist_site]]
      t <- table(factor(temp_v, levels = seq(min(temp_v) - 2, max(temp_v) + add)))
      t[1:2] <- t[3] + 1
      locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
      adjoff <- locmax[which.min(abs(locmax - bestoff))]
      ifelse(length(adjoff) != 0, adjoff, bestoff)
    } 
    
    best_from5_tab <- offset_temp[, list(perc = sum(percentage)), offset_from_5
                                  ][perc == max(perc)]
    best_from3_tab <- offset_temp[, list(perc = sum(percentage)), offset_from_3
                                  ][perc == max(perc)]
    
    if(extremity == "auto" &
       ((best_from3_tab[, perc] > best_from5_tab[, perc] &
         as.numeric(best_from3_tab[, offset_from_3]) <= minlen - 2) |
        (best_from3_tab[, perc] <= best_from5_tab[, perc] &
         as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1)) |
       extremity == "3end"){
      best_offset <- as.numeric(best_from3_tab[, offset_from_3])
      line_plot <- "from3"
      cat(sprintf("best offset: %i nts from the 3' end\n", best_offset))
      adj_tab <- site_sub[, list(corrected_offset_from_3 = adj_off(.SD, "site_dist_end3", 0, best_offset)), by = qwidth]
      offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "qwidth")
      offset_temp[is.na(corrected_offset_from_3), corrected_offset_from_3 := best_offset
                  ][, corrected_offset_from_5 := -corrected_offset_from_3 + qwidth - 1]
    } else {
      if(extremity == "auto" &
         ((best_from3_tab[, perc] <= best_from5_tab[, perc] &
           as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1) |
          (best_from3_tab[, perc] > best_from5_tab[, perc] &
           as.numeric(best_from3_tab[, offset_from_3]) > minlen - 2)) |
         extremity == "5end"){
        best_offset <- as.numeric(best_from5_tab[, offset_from_5])
        line_plot <- "from5"
        cat(sprintf("best offset: %i nts from the 5' end\n", -best_offset))
        adj_tab <- site_sub[, list(corrected_offset_from_5 = adj_off(.SD, "site_dist_end5", 1, best_offset)), by = qwidth]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "qwidth")
        offset_temp[is.na(corrected_offset_from_5), corrected_offset_from_5 := best_offset
                    ][, corrected_offset_from_5 := abs(best_offset)
                      ][, corrected_offset_from_3 := abs(corrected_offset_from_5 - qwidth + 1)]
      }
    }
    
    t <- table(factor(dt$qwidth, levels = lev))
    offset_temp[!is.na(offset_from_5), offset_from_5 := abs(offset_from_5)
                ][, total_percentage := as.numeric(format(round((as.vector(t)/sum(as.vector(t))) * 100, 3), nsmall=4))
                  ][, percentage := as.numeric(format(round(percentage, 3), nsmall=4))
                    ][, sample := n]
    
    setcolorder(offset_temp, c("qwidth", "total_percentage", "percentage", "around_site", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    if(start == TRUE | start == T){
      setnames(offset_temp, c("qwidth", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    } else {
      setnames(offset_temp, c("qwidth", "total_percentage", "stop_percentage", "around_stop", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    }
    
    # plot
      options(warn=-1)
      minlen <- ceiling(quantile(site_sub$qwidth, (1 - cl/100)/2))
      maxlen <- ceiling(quantile(site_sub$qwidth, 1 - (1 - cl/100)/2))
	  plotlist <- list()
	  resultlist <- list()
      for (len in minlen:maxlen) {
        site_temp <- dt[site_dist_end5 %in% seq(-len + 1, 0) & qwidth == len]
        site_tab5 <- data.table(table(factor(site_temp$site_dist_end5, levels = (-len + 1) : (len))))
        site_temp <- dt[site_dist_end3 %in% seq(0, len - 2) & qwidth == len]
        site_tab3 <- data.table(table(factor(site_temp$site_dist_end3, levels = (-len) : (len - 2))))
        setnames(site_tab5, c("distance", "reads"))
        setnames(site_tab3, c("distance", "reads"))
        site_tab5[, distance := as.numeric(as.character(site_tab5$distance))
                  ][, extremity := "5' end"]
        site_tab3[, distance := as.numeric(as.character(site_tab3$distance))
                  ][, extremity := "3' end"]
        final_tab <- rbind(site_tab5[distance <= 0], site_tab3[distance >= 0])
        final_tab[, extremity := factor(extremity, levels = c("5' end", "3' end"))]
        final_tab[,length:=len]
        p <- ggplot(final_tab, aes(distance, reads, color = extremity)) +
          geom_line() +
          geom_vline(xintercept = seq(floor(min(final_tab$distance)/3) * 3, floor(max(final_tab$distance)/3) * 3, 3), linetype = 2, color = "gray90") +
          geom_vline(xintercept = 0, color = "gray50") +
          geom_vline(xintercept = - offset_temp[qwidth == len, offset_from_5], color = "#D55E00", linetype = 2, size = 1.1) +
          geom_vline(xintercept = offset_temp[qwidth == len, offset_from_3], color = "#56B4E9", linetype = 2, size = 1.1) +
          geom_vline(xintercept = - offset_temp[qwidth == len, corrected_offset_from_5], color = "#D55E00", size = 1.1) +
          geom_vline(xintercept = offset_temp[qwidth == len, corrected_offset_from_3], color = "#56B4E9", size = 1.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - len, xmax = -flanking , fill = "#D55E00", alpha = 0.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - 1 , xmax = len - flanking - 1, fill = "#56B4E9", alpha = 0.1) +
          labs(x = "Distance from start (nt)", y = "RPF frequency", title = paste("length=", len, " nt", sep = ""), color= "Extremity") +
          theme_bw(base_size = 20) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        
        if(line_plot == "from3"){
          p <- p + geom_vline(xintercept = best_offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset - len + 1, color = "black", linetype = 3, size = 1.1)
        } else {
          p <- p + geom_vline(xintercept = best_offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset + len - 1, color = "black", linetype = 3, size = 1.1)
        }
        
        p <- p + scale_x_continuous(limits = c(min(final_tab$distance), max(final_tab$distance)),
                             breaks = seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5), 
                             labels = as.character(seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5) + base))
        plotlist[[as.character(len)]] <- p
		resultlist[[as.character(len)]] <- final_tab 
        #ggsave(paste(subplot_dir, "/", len, ".", plot_format, sep = ""), plot = p, width = 15, height = 5, units = "in")
      }
      options(warn=0)

    dt[, c("site_dist_end5", "site_dist_end3") := NULL]
    offset <- rbind(offset, offset_temp)
	results <- rbindlist(resultlist)
	fwrite(results, file=paste0(resultdir, "/", file="psite.txt"), sep="\t")
	fwrite(results, file=paste0(resultdir, "/", file="psite.csv"))
	glist <- lapply(plotlist, ggplotGrob)
	ggsave(paste0(resultdir, "/", file="psite_v2.pdf"), width = 8, height = 2.6*length(plotlist), marrangeGrob(grobs = glist, nrow=length(plotlist), ncol=1, top=NULL))
	#return(plotlist)
}
       
psite_plot <- function(data, flanking = 30, cl = 99) {
    names <- names(data)
    offset <- NULL
    n <- "sample"
    dt <- data
    lev <- sort(unique(dt$qwidth))
   
    if(species == "ecoli_k12" | species == "bsu_168" | species == "pfu_dsm_3638" | species == "hsa_NRC1") {
		dt[, site_dist_end5 := end - start_pos]
		dt[, site_dist_end3 := end - stop_pos]
	} else {
		dt[, site_dist_end5 := start - start_pos]
		dt[, site_dist_end3 := start - stop_pos]
	}

	site_sub5 <- dt[site_dist_end5 <= 50 & site_dist_end5 >= -25]
	offset_temp <- site_sub5[, list(offset_from_5 = as.numeric(names(which.max(table(site_dist_end5))))), by = qwidth]
    # plot
    options(warn=-1)
    #minlen <- ceiling(quantile(site_sub5$qwidth, (1 - cl/100)/2))
    #maxlen <- ceiling(quantile(site_sub5$qwidth, 1 - (1 - cl/100)/2))
	minlen<- min(site_sub5$qwidth)
	maxlen<- max(site_sub5$qwidth)
	plotlist <- list()
	resultlist <- list()
    for (len in minlen:maxlen) {
        site_temp <- dt[qwidth == len]
        site_tab5 <- data.table(table(factor(site_temp$site_dist_end5, levels = -25 : 50)))
        setnames(site_tab5, c("distance", "reads"))		
        site_tab5[, distance := as.numeric(as.character(site_tab5$distance))]
        final_tab <- site_tab5
        final_tab[,length:=len]
        p <- ggplot(final_tab, aes(distance, reads)) +
          geom_line() +
          geom_vline(xintercept = seq(floor(min(final_tab$distance)/3) * 3, floor(max(final_tab$distance)/3) * 3, 3), linetype = 2, color = "gray90") +
          geom_vline(xintercept = 0, color = "gray50") +
          #geom_vline(xintercept = offset_temp[qwidth == len, offset_from_5], color = "#D55E00", linetype = 2, size = 1.1) +
		  #geom_vline(xintercept = offset_list[as.character(len)], color = "black", linetype = 3, size = 1.1) + 
          labs(x = "Distance from start (nt)", y = paste0("RPF frequency (",len," nt)"), title = "", color= "Extremity") +
          theme_bw(base_size = 16) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        
        p <- p + scale_x_continuous(limits = c(min(final_tab$distance), max(final_tab$distance)),
                             breaks = seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5), 
                             labels = as.character(seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5)))
        plotlist[[as.character(len)]] <- p
		resultlist[[as.character(len)]] <- final_tab 
    }
    options(warn=0)

	results <- rbindlist(resultlist)
	fwrite(results, file=paste0(resultdir, "/", file="psite.txt"), sep="\t")
	fwrite(results, file=paste0(resultdir, "/", file="psite.csv"))
	glist <- lapply(plotlist, ggplotGrob)
	ggsave(paste0(resultdir, "/", file="psite.pdf"), width = 8, height = 3*length(plotlist), marrangeGrob(grobs = glist, nrow=length(plotlist), ncol=1, top=NULL))
	
	# plot for stop codon
	site_sub3 <- dt[site_dist_end3 <= 25 & site_dist_end3 >= -50]
	offset_temp <- site_sub3[, list(offset_from_5 = as.numeric(names(which.max(table(site_dist_end3))))), by = qwidth]
    
    options(warn=-1)
	minlen<- min(site_sub3$qwidth)
	maxlen<- max(site_sub3$qwidth)
	plotlist <- list()
	resultlist <- list()
	for (len in minlen:maxlen) {
        site_temp <- dt[qwidth == len]
        site_tab5 <- data.table(table(factor(site_temp$site_dist_end3, levels = -50 : 25)))
        setnames(site_tab5, c("distance", "reads"))		
        site_tab5[, distance := as.numeric(as.character(site_tab5$distance))]
        final_tab <- site_tab5
        final_tab[,length:=len]

        p <- ggplot(final_tab, aes(distance, reads)) +
          geom_line() +
          geom_vline(xintercept = seq(floor(min(final_tab$distance)/3) * 3, floor(max(final_tab$distance)/3) * 3, 3), linetype = 2, color = "gray90") +
          geom_vline(xintercept = 0, color = "gray50") +
          #geom_vline(xintercept = offset_temp[qwidth == len, offset_from_5], color = "#D55E00", linetype = 2, size = 1.1) +
		  #geom_vline(xintercept = offset_list[as.character(len)], color = "black", linetype = 3, size = 1.1) + 
          labs(x = "Distance from start (nt)", y = paste0("RPF frequency (",len," nt)"), title = "", color= "Extremity") +
          theme_bw(base_size = 16) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        
        p <- p + scale_x_continuous(limits = c(min(final_tab$distance), max(final_tab$distance)),
                             breaks = seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5), 
                             labels = as.character(seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5)))
        plotlist[[as.character(len)]] <- p
		resultlist[[as.character(len)]] <- final_tab 
    }
    options(warn=0)
	results <- rbindlist(resultlist)
	fwrite(results, file=paste0(resultdir, "/", file="psite_stopcodon.txt"), sep="\t")
	fwrite(results, file=paste0(resultdir, "/", file="psite_stopcodon.csv"))
	glist <- lapply(plotlist, ggplotGrob)
	ggsave(paste0(resultdir, "/", file="psite_stopcodon.pdf"), width = 8, height = 3*length(plotlist), marrangeGrob(grobs = glist, nrow=length(plotlist), ncol=1, top=NULL))
	#return(plotlist)
}	   

frame_psite = function (df, samplename="sample", region = "all", length_range = "all") 
{ 
    if (!region %in% c("all", "cds", "5utr", "3utr")) {
        warning("region is invalid. Set to default \"all\"\n")
        region = "all"
    }
    if (!identical(length_range, "all") & !inherits(length_range, 
        "numeric") & !inherits(length_range, "integer")) {
        warning("length_range is invalid. Set to default \"all\"\n")
        length_range = "all"
    }
	if (exists("final_frame_df")) {
      rm(final_frame_df)
	}
	if (region == "all") {
		if (length_range[1] != "all") {
			df <- df[match_len %in% length_range,]
		}
		df <- df[start_pos != 0 & stop_pos != 0,]
		df[,frame := site_from_start%%3]
		#frame_df <- as.data.table(table(df$frame, factor(df$psite_region, levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", "CDS", "3' UTR"))))		
		frame_df = df[,sum(freq),by=.(frame, psite_region)] 
		frame_df$psite_region = factor(frame_df$psite_region,levels = c("5utr", "cds", "3utr"), labels = c("5' UTR", "CDS", "3' UTR"))			
		setnames(frame_df, c("frame", "region", "count"))
		for(ele in c("5' UTR", "CDS", "3' UTR")) {
			if(!ele %in% frame_df$region){
				frame_tem = data.table(frame=c(0,1,2),region=rep(ele,3),count=rep(0,3))
				frame_df <- rbindlist(list(frame_df,frame_tem))
			}
			
		}
		#frame_df[, perc := count/rep(by(count, region, function(x) sum(x)), each = 3) * 100]
		#frame_df[, perc := count/rep(by(count, region, function(x) sum(x)), each = 3)]
		frame_df[,perc := count/sum(count),by="region"]
		frame_df[is.nan(perc),sum:=0]
	} else {
		df <- data[site_region == region,]
		df[,frame := site_from_start%%3]
		#frame_df <- as.data.table(table(df$frame))		
		frame_df = df[,sum(freq),by=.(frame)] 		
		setnames(frame_df, c("frame", "count"))
		#frame_dfperc <- (frame_df$count/sum(frame_df$count)) * 100
		frame_dfperc <- (frame_df$count/sum(frame_df$count))
	}
	
	frame_by_length_dt <- df[psite_region=="cds",.(freq = sum(freq)),by=.(width,frame)]
	setorderv(frame_by_length_dt, c("width","frame"))
	setnames(frame_by_length_dt,c("Length","Frame","Frequency"))
	
	frame_df$samp <- samplename
	if (exists("final_frame_df")) {
		final_frame_df <- rbind(final_frame_df, frame_df)
	} else {
		final_frame_df <- frame_df
	}
	rm(frame_df)

	final_frame_df$region =factor(final_frame_df$region,levels = c("5' UTR", "CDS", "3' UTR"))
    if (region == "all") {
        plottitle_region <- NULL
    } else {
        if (region == "5utr") {
            plottitle_region <- "Region: 5' UTR"
        }
        if (region == "cds") {
            plottitle_region <- "Region: CDS"
        }
        if (region == "3utr") {
            plottitle_region <- "Region: 3' UTR"
        }
    }
    if (length_range[1] == "all") {
        plottitle_range <- NULL
    } else {
        if (min(length_range) == max(length_range)) {
            plottitle_range <- paste("Read length: ", min(length_range), 
                " nts", sep = "")
        } else {
            if (identical(length_range, min(length_range):max(length_range)) | 
                identical(length_range, seq(min(length_range), 
                  max(length_range), 1))) {
                plottitle_range <- paste("Read lengths: ", min(length_range), 
                  "-", max(length_range), " nts", sep = "")
            } else {
                plottitle_range <- paste("Read lengths: ", paste(length_range, 
                  collapse = ","), " nts", sep = "")
            }
        }
    }
	final_frame_df$frame <- as.factor(final_frame_df$frame)
    plot <- ggplot(final_frame_df, aes(x = frame, y = perc, fill=frame)) + 
        geom_bar(stat = "identity") + theme_bw(base_size = 16) + 
        labs(x = "", y = "Normalized read number", title = paste(plottitle_region, 
		 plottitle_range, sep = "; ")) + scale_fill_manual(values=c("#7cb5ec","#3a3a3c","#90ed7d"))
    if (region == "all") {
        plot <- plot + facet_grid(. ~ region)
    } else {
        plot <- plot + facet_wrap(~samp, ncol = 3)
    }
	fwrite(final_frame_df, file=paste0(resultdir, "/", file="frameStat.txt"), sep="\t")
	fwrite(frame_by_length_dt, file=paste0(resultdir, "/", file="frameStatByLength.txt"), sep="\t")
	fwrite(final_frame_df, file=paste0(resultdir, "/", file="frameStat.csv"))
	fwrite(frame_by_length_dt, file=paste0(resultdir, "/", file="frameStatByLength.csv"))
    #ret_list <- list()
    #ret_list[["df"]] <- final_frame_df
    #ret_list[["plot"]] <- plot
    return(plot)
}


metaprofile_psite = function (df, samplename="Sample", scale_factors = NULL, length_range = "all", transcripts = NULL, utr5l = 25, cdsl = 120, utr3l = 25, plot_title = "auto")
{
    #rownames(annotation) <- as.character(annotation$transcript)
    #l.transcripts <- rownames(annotation)[which(annotation$l_utr5 >= utr5l & annotation$l_cds >= 2 * (cdsl + 1) & annotation$l_utr3 >= utr3l)]		
	#df = df[utr5_len>=utr5l & cds_len >= 2 * (cdsl + 1) & utr3_len >= utr3l,]
	df = df[cds_len >= 2 * (cdsl + 1),]
    if (length(transcripts) != 0) {
        df["transcriptID" %in% transcripts,]
    }
	ntr <- length(df$transcriptID)
	min_reads = 20
	df = df[gene_cov>min_reads]
	df[,normal_value := coverage/gene_cov]
	
    if (!identical(length_range, "all") & !inherits(length_range, 
        "numeric") & !inherits(length_range, "integer")) {
        warning("length_range is invalid. Set to default \"all\"\n")
        length_range = "all"
    }
	if (identical(length_range, "all")) {
		start.sub <- df[site_from_start %in% seq(-utr5l, cdsl), ]
		stop.sub <- df[site_from_stop %in% seq(-cdsl, utr3l), ]
	} else {
		start.sub <- df[site_from_start %in% seq(-utr5l, cdsl) & match_len %in% length_range, ]
		stop.sub <- df[site_from_stop %in% seq(-cdsl, utr3l) & match_len %in% length_range, ]
	}	
	start.tab = start.sub[,list(sitesum = sum(normal_value),sitecount =.N),by=site_from_start]
	start.tab[,normal_value2:=sitesum]
	start.tab = start.tab[order(site_from_start)]
	#region = seq(-utr5l, cdsl)
	#missin = region[which(! region %in% start.tab$site_from_start)]
	#data.table(site_from_start=missin,V1=numeric(length(missin)))	
	#start.tab <- as.data.table(table(factor(start.sub$site_from_start, levels = -utr5l:cdsl)))
	start.tab[,c("sitesum","sitecount") := NULL]
	
	region <- data.table(site_from_start=-utr5l:cdsl)
	start.tab <- merge(region,start.tab,all.x=T,by="site_from_start")
	start.tab[is.na(start.tab)] <- 0
	start.tab[,reg := "start"]
	
	#stop.tab <- as.data.table(table(factor(stop.sub$site_from_stop, levels = -cdsl:utr3l)))
	stop.tab = stop.sub[,list(sitesum = sum(normal_value),sitecount =.N),by=site_from_stop]
	stop.tab[,normal_value2:=sitesum]
	stop.tab = stop.tab[order(site_from_stop)]
	stop.tab[,c("sitesum","sitecount") := NULL]
	
	region <- data.table(site_from_stop=-cdsl:utr3l)
	stop.tab <- merge(region,stop.tab,all.x=T,by="site_from_stop")
	stop.tab[is.na(stop.tab)] <- 0
	stop.tab[,reg := "stop"]
	
	setnames(start.tab,c("distance", "reads", "reg"))
	setnames(stop.tab,c("distance", "reads", "reg"))	
	samp.tab <- rbind(start.tab, stop.tab)
	
	if (! is.null(scale_factors)) {
		samp.tab[,reads := reads * scale_factors]
	}
	if (exists("final.tab.psitemetaprofile")) {
		final.tab.psitemetaprofile$reads <- final.tab.psitemetaprofile$reads + samp.tab$reads
	} else {
		final.tab.psitemetaprofile <- samp.tab
	}

    final.tab.psitemetaprofile$reg <- factor(final.tab.psitemetaprofile$reg, 
        levels = c("start", "stop"), labels = c("Distance from start (nt)", "Distance from stop (nt)"))
    linestart <- data.frame(reg = rep(c("Distance from start (nt)", 
        "Distance from stop (nt)"), times = c(length(c(rev(seq(-3, -utr5l, -3)), seq(3, cdsl, 3))), length(c(rev(seq(-2, -cdsl, -3)), seq(1, utr3l, 3))))), line = c(rev(seq(-3, -utr5l, -3)), seq(3, cdsl, 3), rev(seq(-2, -cdsl, -3)), seq(1, utr3l, 3)))
    linered <- data.frame(reg = c("Distance from start (nt)", "Distance from stop (nt)"), line = c(0, 1))
		
    if (length(plot_title) != 0 && plot_title == "auto") {
        if (identical(length_range, "all")) {
            plot_title <- paste(paste(samplename, collapse = "+"), " (", ntr, " tr)", sep = "")
        } else {
            if (min(length_range) == max(length_range)) {
                plot_title <- paste(paste(samplename, collapse = "+"), " (", ntr, " tr) - Read length: ", min(length_range), " nts", sep = "")
            } else {
                if (identical(length_range, min(length_range):max(length_range)) | identical(length_range, seq(min(length_range),  max(length_range), 1))) {
                  plot_title <- paste(paste(samplename, collapse = "+"), " (", ntr, " tr) - Read lengths: ", min(length_range), "-", max(length_range), " nts", sep = "")
                } else {
                  plot_title <- paste(paste(samplename, collapse = "+"), " (", ntr, " tr) - Read lengths: ", paste(length_range, collapse = ","), " nts", sep = "")
                }
            }
        }
    }
    plot <- ggplot(final.tab.psitemetaprofile, aes(as.numeric(as.character(distance)), 
        reads)) + geom_line(size = 1.05, color = "#7cb5ec") + 
        geom_vline(data = linered, aes(xintercept = line), linetype = 1, 
            color = "red") + labs(x = "", y = "Normalized read number", title = "") + 
        theme_bw(base_size = 16) + theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) + facet_grid(. ~ 
        reg, scales = "free", switch = "x") + theme(strip.background = element_blank(), 
        strip.placement = "outside") + theme(plot.title = element_text(hjust = 0.5)) + 
        geom_vline(data = linestart, aes(xintercept = line), 
            linetype = 3, color = "gray60")
	fwrite(final.tab.psitemetaprofile, file=paste0(resultdir, "/", file="psite_metaprofile.txt"), sep="\t")
	fwrite(final.tab.psitemetaprofile, file=paste0(resultdir, "/", file="psite_metaprofile.csv"))
    output <- list()
    #output[["plot"]] <- plot
    ##output[["df"]] <- final.tab.psitemetaprofile
    return(plot)
}


#' Main function to transform alignment data to ribosome occupancy/coverage data
#' @return Nothing, instead it writes out occupancy data in a .tsv file and in 2 .wig IGV tracks (+ and - strands)
get_coverage <- function(bamfile, offset_list, txlens, uorfile) {
	if(ribotool == "price") {
		uorfs <- read_uORF(uorfile)
	} else {
		uorfs <- read_uORF_ribORF(uorfile)
	}
	resultdir <- dirname(bamfile)	
	txlens2 <- data.table::copy(txlens)	
	param <- ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE, reverseComplement = FALSE, what = c("qname"))
	reads <- as.data.table(readGAlignments(bamfile, index=paste(bamfile, '.bai', sep = ''), param=param))
	reads[, match_len := cigarWidthAlongQuerySpace(cigar, after.soft.clipping =T)]	
	reads[,freq := gsub("^.*_x","",qname, perl=T)]
	reads[,freq := as.numeric(freq)]
	p1 = plotlength(reads,cl=100)
	
	txlens2[,c('tx_id','gene_id','nexon') := NULL]
	setnames(txlens2,c("tx_name"),c("transcriptID"))
	setnames(reads,c("seqnames"),c("transcriptID"))
	reads = merge(reads, txlens2, by='transcriptID')
	reads[,start_pos:=utr5_len+1]
	reads[,stop_pos:=utr5_len+cds_len]
	psite_plot(reads)
	#ggsave(filename=paste0(samplename,".length.plot.pdf"), plot=p[["plot"]])	
	#reads <- reads[match_len > 25 & match_len < 35]	
    reads <- reads[strand == '+']	
	# Get A-site location for all reads 
	reads <- set_offset_variable(reads, offset_list)
	#if(species == "rn6") {
	#	reads[,seqnames := gsub("\\.\\d+$","",seqnames, perl=T)]
	#}	
	
	reads[,site_from_start := ifelse(start_pos == 0, 0, transcript_coordinate - start_pos)]
	reads[,site_from_stop := ifelse(stop_pos == 0, 0, transcript_coordinate - stop_pos)]
	
	reads[stop_pos == 0, psite_region := "NA"]
	reads[site_from_start >= 0 & site_from_stop <= 0, psite_region := "cds"]
	reads[site_from_start < 0 & site_from_stop < 0,psite_region := "5utr"]
	reads[site_from_start > 0 & site_from_stop > 0,psite_region := "3utr"]
    #reads[,psite_region := ifelse(stop_pos == 0, NA, ifelse(site_from_start >= 0 & site_from_stop <= 0, "cds", ifelse(site_from_start < 0 & site_from_stop < 0, "5utr", "3utr")))]
			
    p2 = frame_psite(reads)
	ggsave(paste0(resultdir,'/frameStat_qc.pdf'), p2, width = 8, height = 5)
	ggsave(paste0(resultdir,'/frameStat_qc.png'), p2, width = 8, height = 5)
	
	coverage_table <- reads[, sum(freq), by=list(transcriptID,strand,transcript_coordinate)]
	gene_coverage <- reads[, list(gene_cov=sum(freq)), by=list(transcriptID)]
	# Release memory	
	rm(reads)
	#gc()
	
	setnames(coverage_table, c('V1'), c('coverage'))	
	#cds_utr_len <- setDT(read.table('./cds_fiveutr.loc.txt', header=T,sep="\t"))
	#setnames(cds_utr_len, c('gene'), c('transcriptID'))
	coverage_table = merge(coverage_table, txlens2, by='transcriptID')
	coverage_table[, CDS_coordinate := transcript_coordinate - utr5_len]	
	coverage_cds <- data.table::copy(coverage_table)	
	#coverage_table <- coverage_table[CDS_coordinate > -60]
	#coverage_table <- coverage_table[CDS_coordinate < cdslen+60]
	
	#get RPF count for mRNA based TE
	#coverage_table_for_count <- coverage_table[CDS_coordinate > 45 & CDS_coordinate < cds_len-15]	
	#coverage_count = coverage_table_for_count[, sum(coverage), by='transcriptID']
	#setnames(coverage_count,c("tx_name","count"))
	#coverage_count <- merge(coverage_count, txlensMax[,.(tx_name,gene_id)], by="tx_name", all.x=T)
	#coverage_count <- coverage_count[!is.na(gene_id),.(gene_id,count)]
	#fwrite(coverage_count, file = paste0(resultdir, '/mrna.count.txt'), sep = '\t')
	
	coverage_table[,start_pos := utr5_len+1]
	coverage_table[,stop_pos := utr5_len+cds_len]
	coverage_table[,site_from_start := ifelse(stop_pos == 0, 0, transcript_coordinate - start_pos)]
	coverage_table[,site_from_stop := ifelse(stop_pos == 0, 0, transcript_coordinate - stop_pos)]	
    
	#coverage_table[,psite_region := ifelse(stop_pos == 0, NA, ifelse(site_from_start >= 0 & site_from_stop <= 0, "cds", ifelse(site_from_start < 0 & site_from_stop < 0, "5utr", "3utr")))]
	coverage_table[stop_pos == 0, psite_region := "NA"]
	coverage_table[site_from_start >= 0 & site_from_stop <= 0, psite_region := "cds"]
	coverage_table[site_from_start < 0 & site_from_stop < 0, psite_region := '5utr']
	coverage_table[site_from_start > 0 & site_from_stop > 0, psite_region := '3utr']
	
	coverage_table <- coverage_table[gene_coverage,gene_cov:=i.gene_cov,on='transcriptID'] 
	
	p3 <- metaprofile_psite(coverage_table)
	ggsave(paste0(resultdir,'/metaprofile_qc.pdf'), p3, width = 8, height = 5)
	ggsave(paste0(resultdir,'/metaprofile_qc.png'), p3, width = 8, height = 5)
	
	coverage_cds_count <- coverage_table[,sum(coverage), by=c('transcriptID','psite_region')]
	coverage_cds_count <- coverage_cds_count[psite_region=="cds"]
	coverage_cds_count <- coverage_cds_count[,.(transcriptID,cds_count=V1)]
	coverage_uORF <- merge(coverage_table,uorfs,by="transcriptID",all.y=T, allow.cartesian=TRUE)
	coverage_uORF[transcript_coordinate >=tr_start & transcript_coordinate < tr_end, psite_region := 'uORF']
	coverage_uORF <- coverage_uORF[psite_region == 'uORF']
	coverage_uORF[, uORF_cov := sum(coverage), by=c('transcriptID')]
	coverage_uORF[,uORF_len:=tr_end-tr_start+1]
	coverage_uORF <- unique(merge(coverage_uORF[,.(transcriptID,ID,tr_start,tr_end,cds_len,uORF_len,uORF_cov)],coverage_cds_count,by="transcriptID",all.x=T))
	coverage_uORF[,uORF_ratio := (uORF_cov/uORF_len)/(cds_count/cds_len)]

	if(ribotool == "price") {
		a <- fread(uorfile,head=F)
		setnames(a,c("chr","start","end","ID","type","strand","startcodon","depth","p_value"))
		a <- merge(a,coverage_uORF[,.(ID,transcriptID,tr_start,tr_end,uORF_cov,cds_count,uORF_ratio)],by="ID")
	} else {
		a <- fread(uorfile,head=T)
		a <- a[,.(chr=chrom, start=codon5, end=codon3, ID=orfID, readNum, PME, codonNum, orfscore, strength, p_value=pred.pvalue)]
		uorfs[,num := seq_len(.N),by=type]
		a <- merge(a,uorfs,by="ID")
		a <- merge(a,coverage_uORF[,.(ID,uORF_cov,cds_count,uORF_ratio)],by="ID")
		a[type=="outframe.overlap.uORF",type:="uoORF"]
		a <- a[type!="internal"]
		a[,ID:=paste0(type,"_",num)]
		a <- a[!is.na(cds_count)]
	}
	a <- a[!is.infinite(orfscore)]
	fwrite(a, file = paste0(resultdir, "/uORFs.tsv"),sep="\t")
	fwrite(a, file = paste0(resultdir, "/uORFs.csv"),sep=",")
	jsonlite::write_json(a, paste0(resultdir,"/uORFs.json"))

	# Write coverage data to a tab separated file
	coverage_cds[, totalcov := sum(coverage), by='transcriptID']
	coverage_cds[CDS_coordinate > 0 & CDS_coordinate <= cds_len, cdscov := sum(coverage), by='transcriptID']
	
	coverage_cds[, length_codon:=ceiling(cds_len/3)]
	#coverage_cds[CDS_coordinate > 0 & CDS_coordinate <= cds_len, codon := ceiling((CDS_coordinate)/3)]
	coverage_cds[, codon := ceiling((CDS_coordinate)/3)]
	coverage_cds[codon <= 0, codon := NA]
	
	cov_codon <- coverage_cds[, sum(coverage), by=c('transcriptID', 'codon')]
	coverage_cds[cov_codon, coverage_codon := i.V1, on=c('transcriptID', 'codon')]
	coverage_cds[is.na(coverage_cds$coverage_codon), coverage_codon:=0]
	
	fwrite(coverage_cds, file = paste0(resultdir, "/coverage_cds.tsv"))
	rm(coverage_cds)
	#gc()
	#res = list(p1,p2,p3)	
	#return(res)
}

read_uORF <- function(uorfile){
	transcriptRange <- exonsBy(txdb,by="tx",use.names=T)
	a <- fread(uorfile,head=F)
	#if(nrow(a) <= 1){
	#   get_log_status("No uORF detected in your data","Error",70);
	#   exit;
	#}
	a[,chr:=V1]
	a[chr=="chrMT",chr:="chrM"]
	a[,start:=V2]
	a[,end:=V3]
	setnames(a,"V6","strand")
	setnames(a,"V7","startCodon")
	#grl <- makeGRangesListFromDataFrame(a, split.field = "V4",names.field = "V4")

	b <- a[,.(chr,start,end,V4,V5,strand,startCodon)]
	b[,end:=start]
	c <- a[,.(chr,start,end,V4,V5,strand,startCodon)]
	c[,start:=end]

	gr1 <- makeGRangesFromDataFrame(b)
	gr2 <- makeGRangesFromDataFrame(c)

	overlaps <- mapToTranscripts(gr1, transcriptRange)
	varclashm <- as.data.table(overlaps)
	varclashm <- varclashm[,.(seqnames,start,strand,xHits)]
	setnames(varclashm, c("transcriptID","tr_start","strand","xHits"))
	b[, xHits := .I]
	b <- b[, .(ID = V4, Type=V5, startCodon, xHits)]
	b[,id2:=gsub("_start","",ID)]
	b[,tx:=gsub("_.*","",ID)]
	b[,tx:=gsub(".*\\|","",tx)]
	varclashm <- merge(varclashm, b, by="xHits", all.x=T)
	varclashm[,c("xHits") := NULL]
	varclashm <- varclashm[transcriptID==tx]

	overlaps <- mapToTranscripts(gr2, transcriptRange)
	varclashm2 <- as.data.table(overlaps)
	varclashm2 <- varclashm2[,.(seqnames,end,strand,xHits)]
	setnames(varclashm2, c("transcriptID","tr_end","strand","xHits"))
	c[, xHits := .I]
	c <- c[, .(ID = V4, Type=V5, startCodon, xHits)]
	c[,id2:=gsub("_end","",ID)]
	c[,tx:=gsub("_.*","",ID)]
	c[,tx:=gsub(".*\\|","",tx)]
	varclashm2 <- merge(varclashm2, c, by="xHits", all.x=T)
	varclashm2[,c("xHits") := NULL]
	varclashm2 <- varclashm2[transcriptID==tx]

	res <- merge(varclashm,varclashm2,by="id2",all.x=T)
	res <- res[,.(transcriptID=tx.x, tr_start, tr_end, strand=strand.x, ID=id2, type=Type.x, startCodon=startCodon.x)]
	res[]
	res[strand=="-", c("tr_start", "tr_end") := .(tr_end, tr_start)]
	return(res)
}
 

read_uORF_ribORF <- function(uorfile){
	#chr1    37357374        37357831        ENSMUST00000027287.10_uORF_1    uORF    +       ATG     31.1    0.00080838
	a <- fread(uorfile,head=T)
	a1 <- a[, tstrsplit(orfID, "|", fixed=TRUE)]
	a2 <- a1[, tstrsplit(V1, ":", fixed=TRUE)]
	setnames(a2,c("transcriptID","chr","strand"))
	a3 <- a1[,  tstrsplit(V3, ":", fixed=TRUE)]
	setnames(a3,c("tr_len","tr_start","tr_end"))
	a <- cbind(a[,.(orfID=orfID)],a1,a2,a3)
	b <- a[,.(transcriptID,tr_start,tr_end,strand,ID=orfID,type=V4,startCodon=V5)]
	b[,tr_start:=as.numeric(tr_start)]
	b[,tr_end:=as.numeric(tr_end)]
	b[,tr_end:=tr_end-1]
	return(b)
}

get_ncrna_coverage <- function(bamfile, offset_list, txlens) {	
	txlens2 <- data.table::copy(txlens)	
	param <- ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE, reverseComplement = FALSE, what = c("qname"))
	reads <- as.data.table(readGAlignments(bamfile, index=paste(bamfile, '.bai', sep = ''), param=param))
	reads[, match_len := cigarWidthAlongQuerySpace(cigar, after.soft.clipping =T)]	
	reads[,freq := gsub("^.*_x","",qname, perl=T)]
	reads[,freq := as.numeric(freq)]
	#reads <- reads[match_len > 25 & match_len < 35]	
    reads <- reads[strand == '+']	
	reads <- set_offset_variable(reads, offset_list)
	txlens2[,c('tx_id','gene_id','nexon') := NULL]
	setnames(txlens2,c("tx_name"),c("transcriptID"))
	setnames(reads,c("seqnames"),c("transcriptID"))
	reads = merge(reads, txlens2, by='transcriptID')
	reads[,start_pos:=utr5_len+1]
	reads[,stop_pos:=utr5_len+cds_len]
	reads[,site_from_start := ifelse(start_pos == 0, 0, transcript_coordinate - start_pos)]
	reads[,site_from_stop := ifelse(stop_pos == 0, 0, transcript_coordinate - stop_pos)]
	reads[stop_pos == 0, psite_region := "NA"]
	reads[site_from_start >= 0 & site_from_stop <= 0, psite_region := "cds"]
	reads[site_from_start < 0 & site_from_stop < 0,psite_region := "5utr"]
	reads[site_from_start > 0 & site_from_stop > 0,psite_region := "3utr"]
		
	coverage_table <- reads[, sum(freq), by=list(transcriptID,strand,transcript_coordinate)]
	rm(reads)
	#gc()
	
	setnames(coverage_table, c('V1'), c('coverage'))	
	coverage_table = merge(coverage_table, txlens2, by='transcriptID')
	coverage_table[, CDS_coordinate := transcript_coordinate - utr5_len]	
	coverage_ncRNA <- data.table::copy(coverage_table)	
	
	# Write coverage data to a tab separated file
	coverage_ncRNA[, totalcov := sum(coverage), by='transcriptID']
	coverage_ncRNA[CDS_coordinate > 0 & CDS_coordinate <= cds_len, cdscov := sum(coverage), by='transcriptID']	
	coverage_ncRNA[, length_codon:=ceiling(cds_len/3)]
	#coverage_ncRNA[CDS_coordinate > 0 & CDS_coordinate <= cds_len, codon := ceiling((CDS_coordinate)/3)]
	coverage_ncRNA[, codon := ceiling((CDS_coordinate)/3)]
	coverage_ncRNA[codon <= 0, codon := NA]
	
	cov_codon <- coverage_ncRNA[, sum(coverage), by=c('transcriptID', 'codon')]
	coverage_ncRNA[cov_codon, coverage_codon := i.V1, on=c('transcriptID', 'codon')]
	coverage_ncRNA[is.na(coverage_ncRNA$coverage_codon), coverage_codon:=0]
	#coverage_cds <- fread(paste0(resultdir, "/coverage_cds.tsv"))
	#coverage_all <- rbind(coverage_cds,coverage_ncRNA)
	#fwrite(coverage_all, file = paste0(resultdir, "/coverage_all.tsv"))
	fwrite(coverage_ncRNA, file = paste0(resultdir, "/coverage_ncRNA.tsv"))
	#gc()
}


###### DEFINE FUNCTIONS ######
as_occupancy_profiles <- function(coverage_cds) {
	
	tmp <- coverage_cds[!is.na(codon), .(list(codon), list(coverage_codon),list(integer(ceiling(length_codon[[1]])))), by=transcriptID]
	setnames(tmp, c('transcriptID', 'codon', 'coverage', 'occupancy'))

	profiles <- tmp[, {.(list(mapply(function(x,y,z) {tmp<-unlist(x); tmp[y]=z; return(tmp)}, x=occupancy, y=codon, z=coverage)))}, by=transcriptID]

	setnames(profiles, c('transcriptID', 'profile'))
	profiles[, profile:=sapply(profile, as.integer)]
	
	annot <- unique(coverage_cds[, -c('transcript_coordinate', 'CDS_coordinate', 'coverage', 'codon', 'coverage_codon')])
	output <- merge(profiles, annot, by='transcriptID')
	return(output)
}

wellexp_ingolia <- function(profile) {
	windows <- split(profile, ceiling(seq_along(profile)/5)) 
	sum_windows <- sapply(windows, sum)
	median_windows <- median(sum_windows)
	return(median_windows)
}

metagene_profile <- function(dt, prefix, window, min_reads, what='start') {

	dt <- dt[length_codon > window]
	
	dt[, sum := sapply(profile,sum)]
	
	dt <- dt[sum >= min_reads]
	dt[, profile_norm := mapply(function(x,y){x/y}, x=profile, y=sum)]
	
	if (what=='stop') dt[, profile_norm := lapply(profile_norm, rev)] 
	
	dt[, profile_nostop := lapply(profile_norm, function(x){head(x, -5)})]

	g <- dt[, mean(sum)]/nrow(dt)
	
	max_len <- max(sapply(dt$profile_nostop,length))
	metaprofile <- numeric(max_len)
	
	for (i in seq(1:max_len)) {
		v <- sapply(dt$profile_nostop, '[', i) 
		metaprofile[[i]] <- sum(v, na.rm = T)*g 
	}
	
	metaprofile <- c(0,metaprofile) # Add a 'false' 0 position with 0 counts to improve visualization
	meta_dt <- data.table(metaprofile) # Convert vector to data.table for ggplot
	meta_dt[,position:=(seq(meta_dt[,.N])-1)] # Add position data to the data.table for plotting
	meta_dt <- meta_dt[position < window] # Keep only positions to be plotted
	if (what=='stop') meta_dt[, position:=-position] # Reverse X axis for 'stop' metagene plot
	
	x_lab <- ifelse(what=='stop', 'Distance from stop (codons)', 'Distance from start (codons)') # X axis label
	
	# Make the graph using ggplot2
	ggplot(meta_dt, aes(position)) + xlab(x_lab) + ylab('Average ribosome occupancy') + geom_line(aes(y=metaprofile, colour='orange')) + theme(panel.background = element_blank())
	
	# Save metagene plot as TIFF image file
	filename <- ifelse(what=='stop', paste0(prefix, '_metagene_stop.pdf'),
	paste0(prefix, '_metagene_start.pdf'))
	ggsave(filename = filename, width = 7, height = 4)
	#return(rp)
	
}

codon_usage_count <- function(coverage_cds, d, normalization) {

	# Keep only in-frame reads
	coverage_cds <- coverage_cds[(CDS_coordinate-1)%%3 == 0]
	
	#coverage_cds <- coverage_cds[codon > 15 & codon < length_codon-5]
	coverage_cds <- coverage_cds[codon > 15]
	
	##add seq
	coverage_cds[transcript_seqs, cdsseq:=i.seq, on=.(transcriptID)]
	siteNames <- c("E","P","A","+1","+2","+3")
	names(siteNames) <- -1:4
	mutialinfor <- function(i) {
	   coverage_cds[strand=='+', V1 := mapply(function(seqs,codon) toupper(substr(seqs, codon, codon+2)), cdsseq, utr5_len+CDS_coordinate+(i*3))]
	   aatem <- coverage_cds[,list(V1,coverage)]
	   setnames(aatem, "V1",paste0('position_', siteNames[as.character(i)]))
	}
	require(parallel)
    aa <- mclapply(-1:4, mutialinfor, mc.cores = 6)
	coverage_cds[,cdsseq:=NULL]
	
	coverage_cds[, freq := 1]
	# Calculate sum of raw reads by codon sequence
	summary <- list()
	for (i in 1:length(aa)) {
		tem <- aa[[i]]
		pos <- names(tem)[which(names(tem) %like% 'position')]
		if(pos == "position_A") {
			codon_A_sites = tem[,position_A]
		}
		summary[[i]] <- setnames(tem[, sum(coverage), by=eval(as.character(pos))], c('codon', pos))
	}

	codon_usage <- Reduce(function(x,y) merge(x,y,all=T), summary) 
	codon_usage[is.na(codon_usage)] = 0
	codon_usage <- codon_usage[grepl('[ATGC]{3}', codon_usage$codon)]
	
	# Normalize by +1 to +3
	codon_usage[, baseline:=rowMeans(codon_usage[, paste0('position_', '+', c(1:3))])]
	norm_codon_usage<-sapply(codon_usage[, names(codon_usage) %like% 'position', with=F], function(x) x/codon_usage$baseline)
	codon_usage <- cbind(codon_usage[,1], norm_codon_usage)
	

	# Calculate read count for each codon as percentage
	codon_usage <- as.data.table(append(codon_usage, list(aminoacid=as.character(Biostrings::translate(DNAStringSet(codon_usage$codon)))), after = 1))
	
	# Keep CDS with at least 50 codons after removing beggining and end codons
	#coverage_cds <- coverage_cds[(length/3)>(2*d)+50][codon > d & codon < (ceiling(length/3)-d)]
	
	coverage_cds$position_A <- codon_A_sites
	if (normalization == 'gene_avg_density') {
		# Calculate total footprints for each CDS
		stats <- coverage_cds[, sum(coverage), by='transcriptID']
		times <- coverage_cds[, sum(freq), by='transcriptID']
		coverage_cds[stats, total_rpf := i.V1, on='transcriptID']
		coverage_cds[times, times := i.V1, on='transcriptID']
		coverage_cds[, coverage := coverage/(total_rpf/times)]
	}
	
	codon_occuracy_matrix = coverage_cds[,.(list(coverage)),by='position_A']
	setnames(codon_occuracy_matrix,c("position_A","V1"),c("codon","occupacy_metric"))
	codon_usage[codon_occuracy_matrix, occupacy_metric := i.occupacy_metric, on="codon"]
	
	# Output ordered results table
	return(codon_usage)
}


as_occupancy_profiles_bin <- function(cov_list) {
	peak_table <- cov_list
	peak_table <- peak_table[cds_len >= 100 & CDS_coordinate > 0 & CDS_coordinate <= cds_len]		
	#peak_table <- peak_table[(CDS_coordinate-1)%%3 == 0]	
	#peak_table = peak_table[name %in% objgenes]		
	#coverage_by_gene <- peak_table[, list(sumcount=.N), by=transcriptID]
	#peak_table= peak_table[coverage_by_gene,on="transcriptID"]
	#peak_table[,normalized := ceiling(coverage*10/(gene_coverage/sumcount)),]
	#coverage_table = merge(coverage_table, coverage_by_gene, by='transcriptID')	
	peak_table <- peak_table[length_codon>=100]	
	peak_table <- unique(peak_table[,.(transcriptID,strand,cdscov,length_codon, codon, coverage_codon)])
	peak_table[, rel_pos := codon/length_codon]
	#peak_table[, nor_coverage := (coverage/totalcov)]	
	#min_reads = 64
	#peak_table = peak_table[totalcov>min_reads]
	peak_table[, nor_coverage := coverage_codon/cdscov]
	peak_table <- peak_table[rel_pos >= 0 & rel_pos <= 1]
	peak_table[, bin_pos := ceiling(rel_pos*100)*0.01]
	
	peak_table <- peak_table[,sum(nor_coverage),by=c("transcriptID","bin_pos")]
	
	results <- peak_table[,list(mean(V1)),by=bin_pos]
	#results = peak_table
	#results = results[bin_pos>0.01 & bin_pos <0.98]
	setnames(results,"V1","Normalized_coverage")	
	#results[,eval(n) := coverage_codon/sum(coverage_codon)]	
	filename <- names(peak_table)
	return(results)
}

as_occupancy_profiles_for_start <- function(cov_list) {
	peak_table = cov_list	
	#peak_table = peak_table[cds_len!=0,]
	
	#peak_table <- peak_table[cds_len>=100]
	#peak_table <- peak_table[(CDS_coordinate-1)%%3 == 0]
	if(species == "sce_R64" | species == "ecoli_k12" | species == "bsu_168" | species == "pfu_dsm_3638" | species == "hsa_NRC1") {
		peak_table <- peak_table[length_codon>=100]	
	} else {
		peak_table <- peak_table[utr5_len > 60 & length_codon>=100]	
	}
	
	peak_table <- unique(peak_table[,.(transcriptID,strand,cdscov,length_codon, codon, coverage_codon)])
	
	cov_for_norm = peak_table[codon >0 & codon <=100, sum(coverage_codon), by=c('transcriptID')] 	
	peak_table[cov_for_norm, coverage_500 := i.V1, on=c('transcriptID')]
	peak_table[,coverage_500 := coverage_500/500]
	#peak_table = peak_table[name %in% objgenes]	
	###peak_table <- peak_table[length_codon >= 300]	
	#coverage_by_gene <- peak_table[, list(sumcount=.N), by=transcriptID]
	#peak_table= peak_table[coverage_by_gene,on="transcriptID"]
	#peak_table[,normalized := ceiling(coverage*10/(gene_coverage/sumcount)),]
	#coverage_table = merge(coverage_table, coverage_by_gene, by='transcriptID')		
	peak_table[, nor_coverage := coverage_codon/cdscov]
	peak_table[,codon_coordinate := codon]
	peak_table <- peak_table[codon > 0 & codon <= 100]
	peak_table = peak_table[, mean(nor_coverage), by=c('codon_coordinate')]
	#peak_table[is.na(V1),V1:=0]
	setnames(peak_table,"V1","Normalized_coverage")
	##peak_table = peak_table[transcriptID %in% txlens$tx_name,]
	return(peak_table)
}


as_occupancy_profiles_for_end <- function(cov_list) {
	peak_table = cov_list
	#peak_table = peak_table[cds_len!=0,]
	#peak_table <- peak_table[cds_len>=300]
	#peak_table <- peak_table[(CDS_coordinate-1)%%3 == 0]	
	if(species == "sce_R64" | species == "ecoli_k12" | species == "bsu_168" | species == "pfu_dsm_3638" | species == "hsa_NRC1") {
		peak_table <- peak_table[length_codon>=100]
	} else {
		peak_table <- peak_table[length_codon>=100 & utr3_len > 100]
	}
	peak_table <- unique(peak_table[,.(transcriptID,strand,cdscov,length_codon, codon, coverage_codon)])
		
	#peak_table = peak_table[name %in% objgenes]	
	###peak_table <- peak_table[length_codon >= 300]	
	#coverage_by_gene <- peak_table[, list(sumcount=.N), by=transcriptID]
	#peak_table= peak_table[coverage_by_gene,on="transcriptID"]
	#peak_table[,normalized := ceiling(coverage*10/(gene_coverage/sumcount)),]
	#coverage_table = merge(coverage_table, coverage_by_gene, by='transcriptID')		
	peak_table[, nor_coverage := coverage_codon/cdscov]
	
	peak_table <- peak_table[codon >= length_codon-100 & codon <= length_codon+50]
	peak_table[,codon_coordinate := codon-length_codon]
	peak_table = peak_table[, sum(nor_coverage), by=c('codon_coordinate')]
	peak_table <- peak_table[!is.na(V1)]
	setnames(peak_table,"V1","Normalized_coverage")
	##peak_table = peak_table[transcriptID %in% txlens$tx_name,]
	return(peak_table)
}


#for metacodon analysis
codon_occupancy_profiles <- function(cov_list, transcript_seqs){

    coverage_cds <- cov_list
	coverage_cds <- coverage_cds[CDS_coordinate <= cds_len & cds_len > 61, ]
	tmp <- coverage_cds[!is.na(codon), .(list(CDS_coordinate), list(coverage),list(integer(ceiling(cds_len[[1]])))), by=transcriptID]
	setnames(tmp, c('transcriptID', 'cdscor', 'coverage', 'occupancy'))
	profiles <- tmp[, {.(list(mapply(function(x,y,z) {tmp<-unlist(x); tmp[y]=z; return(tmp)}, x=occupancy, y=cdscor, z=coverage)))}, by=transcriptID]
	setnames(profiles, c('transcriptID', 'profile'))
	profiles[, profile:=sapply(profile, as.integer)]
	txlens2 <- data.table::copy(txlens)
	setnames(txlens2, "tx_name", "transcriptID")	
	profiles = merge(x=profiles,y=txlens2[,c("transcriptID","cds_len","utr5_len","tx_len")], by="transcriptID", all.x=T)
	#profiles = unique(profiles,by = c('transcriptID'))
	#profiles [,"codon_seqs" := mapply(function(seq,ft) {substring(x, seq(1,nchar(seq),2), seq(2,nchar(seq),2))}, cdsseq, utr5_len)]
	#profiles [,"codon_seqs" := mapply(function(seq,ft,cl) {laply(seq(ft,ft+cl,3), function(i) substr(seq, i, i+2))}, cdsseq, utr5_len, cds_len)]
	profiles[transcript_seqs, cdsseq:=i.seq, on=.(transcriptID)]	
	##filter
	#profiles = profiles[nchar(cdsseq) != tx_len, ]
	profiles[,"cds" := mapply(function(seq,i,j) {substr(seq,i+1,i+j)}, cdsseq, utr5_len, cds_len)]
	profiles = profiles[startsWith(profiles$cds,"atg"),] 	
	#profiles[,numx := sapply(cds, function(x){length(which(alphabet(translate(DNAString(x))) == "*"))})]	
	profiles[,"codon_seqs" := mapply(function(seq,i,j) {strsplit(gsub("([[:alnum:]]{3})", "\\1 ", substr(seq,i+1,i+j)), " ")[[1]]}, cdsseq, utr5_len, cds_len)]

	profiles[,c('cdsseq','cds') := NULL]		
	data_summary1 <- function(profile, codon_seqs, i, transcriptID) {
		cs = unlist(codon_seqs)
		x = unlist(profile)
		if(length(x)<61){
		   cat(paste0(i," ",length(x)))
		   return()
		}
		if(as.numeric(i) %% 1000 == 0) {
		   cat(paste0(i,"\n"))
		}
		m <- rollapply(x, 61, by = 3, c)
		ccc = split(m, row(m))
		#res = data.table(codons_seq=head(tail(cs, -30), -30),windows=ccc, transcriptID=transcriptID)
		csl = length(cs)-10
		if(length(11:csl) != length(ccc)) {
		    #cat(i)
            #cat("\n")
 			return()
		}
		res = data.table(codons_seq=cs[11:csl],windows=ccc, transcriptID=transcriptID)
		return(res)
	}	
	data_summary2 <- function(i, profile, codon_seqs, transcriptID) {
		cs = unlist(codon_seqs[i])
		x = unlist(profile[i])
		if(length(x)<61){
		   cat(paste0(i," ",length(x)))
		   return()
		}
		if(i %% 1000 == 0) {
		   cat(paste0(i,"\n"))
		}
		m <- rollapply(x, 61, by = 3, c)
		ccc = split(m, row(m))
		csl = length(cs)-10
		if(length(11:csl) != length(ccc)) {
 			return()
		}
		res = data.table(codons_seq=cs[11:csl],windows=ccc, transcriptID=transcriptID[i])
		return(res)
	}
	
    cct <- mcmapply(data_summary1, profiles$profile, profiles$codon_seqs, rownames(profiles), profiles$transcriptID, SIMPLIFY = FALSE, mc.cores = 6, mc.cleanup = TRUE)
	#cct <- mclapply(seq_along(profiles$profile), data_summary2, profiles$profile, profiles$codon_seqs, profiles$transcriptID, mc.cores = 6, mc.preschedule = FALSE)	
	rl = rbindlist(cct)
	rl[,sum:=sapply(windows,sum)]
	#rl = rl[sum>50,]
	rl = rl[sum>20,]
	#rl = merge(x=rl,y=txlens[,c("transcriptID","cds_len")], by="transcriptID", all.x=T)	

	rl[,normalizedValue := mapply(function(x){x = unlist(x);list(x*length(x)/sum(x))},x=windows)]		
	customfunc <- function(dt){
	  #library(plyr)
      #df <- ldply (l, data.frame)
	  q = do.call(rbind.data.frame, dt$normalizedValue)
	  m = apply(q,2,mean)
	  return(list(as.vector(m)))
	}
	attg = rl[, list(window_profile = customfunc(.SD)), by=codons_seq, .SDcols=c("normalizedValue")]
	setnames(attg,"window_profile","normalized_value")	
	return(attg)
}

get_log_status <- function(messages, status, percent){
	logdate <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
	messages <- gsub("\r?\n|\r", "", messages)
	messages <- gsub("\"", "'", messages)
	#if(status != "Error") {
		messages <- paste0(messages,"\t",logdate,"\n")
	#}
	writeLines(paste0('{"type":','"singlecase",','"species":','"',species,'",','"current":','"',status,'",','"sample":','"',samplename,'",','"percent":','"',percent,'"}'), paste0(resultdir,"/status.json"))
    write(messages, file=paste0(resultdir,"/log.txt"), append=TRUE, sep="\n")
}


#main
args <- commandArgs(TRUE)
species <- args[1]
jobid <- args[2]
samplename <- args[3]
 
#species <- "mm10"
#jobid <- "KHZAzY2YBVxS8P81"
#samplename <- "KHZAzY2YBVxS8P81"


theme_set(theme_cowplot())
options(bitmapType='cairo')
samplenames <- strsplit(unlist(strsplit(samplename, "[#]")),"[-]")

resultdir = paste0("./data/results/",jobid)
offsets=paste0(resultdir, "/offsets.conf.txt")
if(!file.exists(offsets)) {
  offsets <- "db/offsets.conf.txt";
}

bamfile=paste0(resultdir, "/mRNA.sort.bam")
txdb <- loadDb(paste0("./db/annotation/",species,".gencode.sqlite"))

if(file.exists(paste0(resultdir, "/repre.valid.pred.pvalue.parameters.txt"))) {
	uorfile <- paste0(resultdir, "/repre.valid.pred.pvalue.parameters.txt")
	ribotool <- 'riborf'
} else if(file.exists(paste0(resultdir, "/uorf.df.txt"))) {
	uorfile <- paste0(resultdir, "/uorf.df.txt")
	ribotool <- 'price'
}

###### GET RIBOSOME OCCUPANCY DATA FROM ALIGNMENT FILES ######
# Read offsets list from file
offsets_table <- read.table(offsets, header = T)
offset_list <-offsets_table$p_offset
names(offset_list) <- offsets_table$length

#txlens <- setDT(transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE))
load(paste0("./db/annotation/",species,".txlens.rda"))

# Process all datasets using get_coverage()
get_log_status("Get coverage from BAM file", "Ribosome profiling", 70)
get_coverage(bamfile = bamfile, offset_list, txlens, uorfile)

#glist <- lapply(ptt, ggplotGrob)
#ggsave(paste0(resultdir,"/qc_all.pdf"), width = 16, height = 16, marrangeGrob(grobs = glist, nrow=3, ncol=2, top=NULL, byrow=TRUE))
#ggsave(paste0(resultdir,"/qc_all.png"), width = 16, height = 16, marrangeGrob(grobs = glist, nrow=3, ncol=2, top=NULL, byrow=TRUE))


dataset <- paste0(resultdir, "/coverage_cds.tsv")
cov_list <- fread(dataset, header = T)

txlens <- txlens[cds_len>0]
txlens[, maxlen := max(cds_len), by = gene_id]
txlensMax <- txlens[cds_len==maxlen]
txlensMax <- txlensMax[!duplicated(gene_id)]

transcript_seqs <- read.fasta(paste0("./db/mRNA/",species,".txdb.fa"), seqtype = 'DNA', as.string = T)
transcript_seqs <- data.table(transcriptID=names(transcript_seqs), seq=as.character(transcript_seqs))


#get_log_status("Count RPF from BAM file", "Ribosome profiling", 75)
#get_count_rpf(bamfile = bamfile, offset_list, txlens)

# median_windows_filter <- lapply(occupancy, function(dt) dt[median_windows>2, transcriptID])
get_log_status("Codon related statistics", "Ribosome profiling", 80)
cov_list <- cov_list[transcriptID %in% txlensMax$tx_name]
cov_list <- cov_list[totalcov >20]


#occupancy profiles
occupancy <- as_occupancy_profiles(cov_list) 
occupancy[, mean_rle:= sapply(profile, function(x) {v <- rle(head(tail(x, -15), -5)); v$lengths[v$values == 0] <- 1; profile <- inverse.rle(v);mean(profile)})]
save(occupancy, file = paste0(resultdir,"/occupancy_profiles.RData"))



pdfFiles <- list.files(resultdir,pattern = "*.pdf$", full.names = T, recursive=T)
if(!dir.exists(paste0(resultdir,"/html"))) {
    dir.create(paste0(resultdir,"/html"))
}
for(i in 1:length(pdfFiles)) {
    pdfFile <- pdfFiles[i]
	pdfFileName <- basename(pdfFile)
    outHtmlFile <- paste0(resultdir,"/html/",sub(".pdf",".html",pdfFileName))
	html <- paste0("<style>\n.pdfobject-container {\n		width: 100%;\n		max-width: 1000px;\n		height: 1000px;\n		margin: 2em 0;\n}\n.pdfobject { border: solid 1px #666; }\n</style>\n<script src='./assets/global/plugins/pdfobject.min.js'></script>\n<div id='pdf' class='pdfobject-container'></div>\n<script>\nvar options = {\n	pdfOpenParams: {\n		pagemode: 'thumbs',\n		navpanes: 0,\n		toolbar: 0,\n		statusbar: 0,\n		view: 'FitV'\n	}\n};\nvar myPDF = PDFObject.embed('",pdfFile,"', '#pdf', options);\n</script>\n")
	writeLines(html,outHtmlFile)
}

#fwrite(codon_occupancy,file=paste0(resultdir,"/metacodon_freq.txt"),sep="\t")
