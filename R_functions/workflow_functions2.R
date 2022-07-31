family_freq <- function(TEgff){
  TE <- read.table(TEgff)
  TEfam <- TE$V10
  #remove motif
  no_motif <- str_replace(TEfam, "Motif:", "")
  #exclude everything with '('
  no_ssr <- no_motif[!str_detect(no_motif, "\\(")]
  #exclude "GA-rich" etc
  ltr_rmd <- no_ssr[!str_detect(no_ssr, "-rich")]
  fam_data <- data.frame(table(ltr_rmd))
  fam_data
}

##Reference Length##
ref_length <- function(ref){
  ref <- read.dna(ref, format = "fasta")
  chroms <- head(ref, 7) #excludes mt
  chrom_length <- lengths(chroms)
  sum(chrom_length)
}


##AT-richness##
##use AT rich gff from occultercut
# gives total bp of ATrich, number of blocks, and prop of genome
AT_rich <- function(ATgff, ref){
  gff <- read.table(ATgff, stringsAsFactors = FALSE, sep = "\t")
  total <- sum(gff$V5 - gff$V4 + 1)
  n <- nrow(gff)
  res <- list(bp = total, n = n)
  #adding in ref for proportions
  ref_len <- ref_length(ref)
  non_AT <- ref_len - res$bp
  prop_AT <-(res$bp/ref_len)*100 
  lineage <- str_extract(gff$V1[1], "[[:alpha:]]+")
  df <- data.frame(lineage, ref_len, res$n, res$bp, non_AT, prop_AT)
  names(df)[1:6] <- c("lineage", "ref_length", "AT_n", "AT_bp", "non_AT", "prop_AT")
  df
}
#use python script to run RIP index of AT rich regions.




###########################
####Class composition######

#updated 14/12/2021
class_freq_plot_all <- function(rm_df_aggregated_g){ #generate mega rm df using rbind. of read_rm_new files.
  c_rm_df <- complete(rm_df_aggregated_g, lineage, group, fill = list(mean = 0)) #will complete for u
  ggplot(c_rm_df, aes(x= reorder(group, desc(group)), fill = lineage)) +
    geom_bar(width = 0.8, color = "black", position = position_dodge()) +
    theme_gray(base_size = 18) +
    scale_fill_manual("", values = lineage_cols) +
    ylab("Number of elements") + xlab("TE Class") + coord_flip() 
} 

class_len_plot_all <- function(rm_df_aggregated_g){ #generate mega rm df using rbind. of read_rm_new files.
  c_rm_df <- complete(rm_df_aggregated_g, lineage, group, fill = list(mean = 0)) #will complete for u
  c_rm_df <- c_rm_df %>% group_by(lineage, group) %>% summarise(len = sum(tlen))
  ggplot(c_rm_df, aes(x= reorder(group, desc(group)), y = len, fill = lineage)) +
    geom_col(color = "black", position = position_dodge()) +
    #theme(axis.text.x = element_text(size 18, angle = 90),
    #      axis.title.y = element_blank())
    scale_fill_manual("", values = lineage_cols) +
    ylab("Length (Mb)") + xlab("TE Class") + coord_flip() +
    scale_y_continuous("Length (Mb)", label=function(x) x/1e6) + theme_gray(base_size = 18)
}

psub_plot <- function(mega_rm_df_g){
  c_rm_df <- complete(mega_rm_df_g, lineage, group, fill = list(mean = 0))
  ggplot(c_rm_df, aes(p_sub, fill = group)) +
    #geom_histogram(color = "black") +
    geom_histogram(colour = "black", size = 0.02) +
    #geom_histogram() +
    facet_grid(lineage ~ .) + 
    scale_fill_manual("", values = TE_cols) +
    xlab("Divergence from TE consensi (%)") + ylab("Number of elements") + 
    theme_gray(base_size=18) + xlim(0, 35) + 
    theme(legend.position = "bottom")
  
}

####Mean RIP of AT rich regions

mean_RIP <- function(idx.tsv){
  idx <- read.table(idx.tsv, header = TRUE)
}


fam_lens <- function(fam_lens.tsv){
  lens <- read.table(fam_lens.tsv, header = TRUE, sep = "\t",
                     comment.char = "")
  sep <- lens %>% separate(Family, c("family", "TE_class"),
                           "#")
  fam_lens <- subset(sep, !(TE_class %in% c("Simple_repeat", "Low_complexity", 
                                            "rRNA", "Satellite.dimer", "buffer")))
}

# read in RIP
RIP <- function(RIP_indices.tsv, fam_lens.tsv){
  fam_lens <- fam_lens(fam_lens.tsv)
  idx <- read.table(RIP_indices.tsv, header = TRUE)
  df <- left_join(idx, fam_lens)
  df <- subset(df, !(TE_class %in% c("Simple_repeat", "Low_complexity", 
                                     "rRNA", "Satellite.dimer", "buffer")))
}


# 1kb RIP
RIP_window <- function(RIP_indices.tsv, window.bed){
  RIP_idx <- read.table(RIP_indices.tsv, header = TRUE)
  window_bed <- read.table(window.bed)
  df <- cbind.data.frame(window_bed, RIP_idx)
  df <- df %>% rename(chrom = V1, start = V2, end = V3)
  df
}

AT_coords <- function(ATgff){
  AT <- read.table(ATgff)
  df <- data.frame("chrom" = AT$V1, "start" = AT$V4, "end" = AT$V5)
  df
}


AT_RIP_ideogram <- function(chrom_lens, AT, RIPwd){
  #chrom_lens <- rename(chrom_lens, len = Length)
  ggplot() + geom_rect(data = chrom_lens, aes(xmin=0, xmax=len, ymin=0, ymax=2), colour="black", fill="white") + 
    geom_rect(data = chrom_lens, aes(xmin=0, xmax=len-0, ymin=0, ymax=2), fill="#99d8c9") + 
    geom_rect(data = AT, aes(xmin = start, xmax = end, ymin = 0,
                             ymax = 2), fill = "#d2bfdd") +    
    # scale_fill_manual(values = TE_cols) + 
    facet_wrap( ~chrom, ncol = 1) + ylim(-2,2) +
    scale_x_continuous("Position (Mb)", label=function(x) x/1e6) + 
    scale_y_continuous(limits=c(-2,5)) +
    theme_minimal() +
    theme(panel.spacing = unit(1, "pt"),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          #axis.ticks.y  = element_blank(),
          strip.text = element_blank()
    ) +
    geom_line(RIPwd, mapping = aes(end, RIP_idx_2), size = 0.1)
}



read_rm_new <- function(RMout){
  out <- read_rm(RMout, include_secondary = FALSE)
  no_ssrs <- subset(out, !(tclass %in% c("Simple_repeat", "Low_complexity", 
                                         "rRNA", "Satellite", "buffer", "Satellite.dimer")))
  out_file <- mutate(no_ssrs, lineage = str_extract(no_ssrs$qname, "[[:alpha:]]+"))
  out_file$ID <- paste(out_file$lineage, out_file$ID, sep = "_")#to give unique ID number per run
  out_file
} # tested. works on all 3


read_rm_no_cleaning <- function(RMout){
  out <- read_rm(RMout, include_secondary = TRUE)
  #no_ssrs <- subset(out, !(tclass %in% c("Simple_repeat", "Low_complexity", 
  #                                      "rRNA", "Satellite", "buffer", "Satellite.dimer")))
  out_file <- mutate(out, lineage = str_extract(out$qname, "[[:alpha:]]+"))
  out_file$ID <- paste(out_file$lineage, out_file$ID, sep = "_")#to give unique ID number per run
  out_file
}


read_effector_P <- function(EPtsv){
  ep <- read.table(EPtsv, sep = "\t")
  ep$V1 <- str_replace(ep$V1, "[[:alnum:]]+_\\d+", "")
  ep$V1 <- str_replace(ep$V1, "\\-[[:alnum:]]+\\s[[:alnum:]]+_", "")
  ep <- rename(ep, gene_id = V1, prediction = V2, probability = V3)
  ep
} #tested on epo. works. 

read_SM <- function(SMtsv){
  ep <- read.table(SMtsv, sep = "\t")
  ep$V1 <- str_replace(ep$V1, "[[:alnum:]]+_", "")
  #ep$V1 <- str_replace(ep$V1, "\\-[[:alnum:]]+\\s[[:alnum:]]+_", "")
  ep <- rename(ep, gene_id = V1)
  ep
}


TEs <- function(RMout){
  TEs <- RMout[,c("qname", "qstart", "qend", "tclass")]
  TEs <- rename(TEs, start = qstart, end = qend, chrom = qname)
  TEs <- group_all_tclass(TEs)
  TEs
}


read_gene_gff <- function(gene_gff){
  gff <- read.table(gene_gff, sep = "\t")
  gff$V9 <- str_replace(gff$V9, "ID=[[:alnum:]]+_", "")
  gff$V9 <- str_replace(gff$V9, ";", "")
  gff <- rename(gff, gene_id = V9, start = V4, end = V5, seqid = V1)
  gff
  
}

effector_location <- function(gene_gff, EPtsv){
  genes <- read_gene_gff(gene_gff)
  ep <- read_effector_P(EPtsv)
  df <- inner_join(genes, ep, by = "gene_id")
  effector_location <- df[,c("seqid", "gene_id", "start", "end", "prediction", "probability")]
  effector_location <- rename(effector_location, chrom = seqid)
  effector_location
}


### File PATHS ######################################

closest_genes_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_closest", sep = "/")
  y <- paste(spp, "closest_d_GENES.gff3", sep = "_")
  paste(x, y, sep = "/")
}
#usage: closest_genes_path(spp = "Ety")

RMout_path <- function(spp){ 
  x <- paste(spp, "repeat_outs", sep = "/")
  y <- paste(spp, "TE.out", sep = "_")
  paste(x, y, sep = "/")
}

effector_path <- function(spp){
  x <- paste("Functional_Annotations", spp, "Proteins", "effectorP_results", sep = "/")
  y <- paste(spp, "effectorP.tsv", sep= "_")
  paste(x,y,sep = "/")
}

SM_path <- function(spp){
  x <- paste("Functional_Annotations", spp, sep = "/")
  y <- paste(spp, "secondary_metabolites.tsv", sep = "_")
  paste(x, y, sep = "/")
  #"Functional_Annotations/Ecl/Ecl_secondary_metabolites.tsv"
}

AT_path <- function(spp){
  x <- paste(spp, "AT_rich", sep = "/")
  y <- paste(spp, "AT_rich.gff", sep = "_")
  paste(x, y, sep = "/")
}

fam_lens_path <- function(spp){
  x <- paste(spp, "references", sep = "/")
  y <- paste(spp, "fam_lens.tsv", sep = "_")
  paste(x, y, sep = "/")
}

chrom_lens_path <- function(spp){
  x <- paste(spp, "references", sep = "/")
  y <- paste(spp, "chrom_lens.tsv", sep = "_")
  paste(x, y, sep = "/")
}

RIP_indices_path <- function(spp){
  x <- paste(spp, "RIP_files", sep = "/")
  y <- paste(spp, "con_RIP_indices.tsv", sep = "_")
  paste(x, y, sep = "/")
}

window_indices_path <- function(spp){
  x <- paste(spp, "references", sep = "/")
  y <- paste(spp, "1kb_indices.tsv", sep = "_")
  paste(x, "1kb_windows", y, sep = "/")
}

window_bed_path <- function(spp){
  x <- paste(spp, "references", sep = "/")
  y <- paste(spp, "1kb.bed", sep = "_")
  paste(x, "1kb_windows", y, sep = "/")
}

gene_gff_path <- function(spp){
  x <- paste(spp, "gene_gffs", sep = "/")
  y <- paste(spp, "genes.gff3", sep = "_")
  paste(x, y, sep = "/")
}


ref_path_Ecl <- function(x){
  "Ecl/references/reference_genome/Ecl1605_22_Epichloe_clarkii_1605_22_45692596_v2.fna"
}

ref_path_Ety <- function(x){
  "Ety/references/reference_genome/Ety1756_Epichloe_typhina_1756_33930528_v4.fna"
}

ref_path_Epo <- function(x){
  "Epo/references/reference_genome/Etp76_Epichloe_typhina_var_poae_NFe76_38327242_v1.fna"
}


ref_path <- function(spp){
  x <- paste(spp, "references", "reference_genome", sep = "/")
  y <- paste(spp, "reference.fna", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}


AT_RIP_path <- function(spp){
  x <- paste(spp, "RIP_files", sep = "/")
  y = paste(spp, "AT_indices.tsv", sep = "_")
  paste(x, y, sep = "/")
}


###calling chromosomes
Ety_chrom <- function(chrom_n){
  paste("Ety_1756_", chrom_n, sep = "")
}

Ecl_chrom <- function(chrom_n){
  paste("Ecl_1605_22_", chrom_n, sep = "")
}

Etp_chrom <- function(chrom_n){
  paste("Etp76_", chrom_n, sep = "")
}
###########################################################
############################################################
#############################################################




TE_cols <- c("DNA" = "#272950", 
             "DNA/hAT" = "#3f4174",
             "DNA/MuLE" = "#6a6ba3",
             "DNA/PIF-Harbinger" = "#969ac9",
             "DNA/Tc1-Mariner" = "#a1c5df",
             
             "LINE" = "#5f7872",
             "LINE/Tad1" = "#abcdbf",
             
             "LTR" = "#9f6c71", 
             "LTR/Copia" = "#c39aa2",
             "LTR/Gypsy" = "#e4cbcf",
             
             
             "MITE" = "#e4d2ad",
             "Unknown" = "#888290")



lineage_cols <- c("Ecl" = "#a87d8a", "Ety" = "#87c6b1",  "Epo" = "#f7eabd")


group_all_tclass <- function(df){
  x <- df
  x$tclass[x$tclass == "Unknown.inc"] <- "Unknown"
  
  x$tclass[x$tclass == "DNA?"] <- "DNA"
  x$tclass[x$tclass == "DNA.inc"] <- "DNA"
  x$tclass[x$tclass == "DNA.inc?"] <- "DNA"
  
  x$tclass[x$tclass == "DNA/hAT?"] <- "DNA/hAT"
  x$tclass[x$tclass == "DNA/hAT.inc"] <- "DNA/hAT"
  x$tclass[x$tclass == "DNA/hAT.inc?"] <- "DNA/hAT"
  
  x$tclass[x$tclass == "DNA/MuLE?"] <- "DNA/MuLE"
  x$tclass[x$tclass == "DNA/MuLE.inc"] <- "DNA/MuLE"
  x$tclass[x$tclass == "DNA/MuLE.inc?"] <- "DNA/MuLE"
  x$tclass[x$tclass == "DNA/MULE"] <- "DNA/MuLE"
  x$tclass[x$tclass == "DNA/MULE?"] <- "DNA/MuLE"
  x$tclass[x$tclass == "DNA/MULE.inc"] <- "DNA/MuLE"
  x$tclass[x$tclass == "DNA/MULE.inc?"] <- "DNA/MuLE"
  
  x$tclass[x$tclass == "DNA/PIF-Harbinger?"] <- "DNA/PIF-Harbinger"
  x$tclass[x$tclass == "DNA/PIF-Harbinger.inc"] <- "DNA/PIF-Harbinger"
  x$tclass[x$tclass == "DNA/PIF-Harbinger.inc?"] <- "DNA/PIF-Harbinger"
  
  x$tclass[x$tclass == "DNA/PiggyBac?"] <- "DNA/PiggyBac"
  x$tclass[x$tclass == "DNA/PiggyBac.inc"] <- "DNA/PiggyBac"
  x$tclass[x$tclass == "DNA/PiggyBac.inc?"] <- "DNA/PiggyBac"
  
  x$tclass[x$tclass == "DNA/Tc1-Mariner?"] <- "DNA/Tc1-Mariner"
  x$tclass[x$tclass == "DNA/Tc1-Mariner.inc"] <- "DNA/Tc1-Mariner"
  x$tclass[x$tclass == "DNA/Tc1-Mariner.inc?"] <- "DNA/Tc1-Mariner"
  
  x$tclass[x$tclass == "LINE/Tad1?"] <- "LINE/Tad1"
  x$tclass[x$tclass == "LINE?"] <- "LINE"
  x$tclass[x$tclass == "LINE.inc"] <- "LINE"
  
  
  x$tclass[x$tclass == "LTR?"] <- "LTR"
  x$tclass[x$tclass == "LTR.inc"] <- "LTR"
  x$tclass[x$tclass == "LTR.inc?"] <- "LTR"
  
  x$tclass[x$tclass == "LTR/Gypsy?"] <- "LTR/Gypsy"
  x$tclass[x$tclass == "LTR/Gypsy.inc"] <- "LTR/Gypsy"
  x$tclass[x$tclass == "LTR/Gypsy.inc?"] <- "LTR/Gypsy"
  
  x$tclass[x$tclass == "LTR/Copia?"] <- "LTR/Copia"
  x$tclass[x$tclass == "LTR/Copia.inc"] <- "LTR/Copia"
  x$tclass[x$tclass == "LTR/Copia.inc?"] <- "LTR/Copia"
  
  x$tclass[x$tclass == "MITE?"] <- "MITE"
  df <- df
  df$group <- x$tclass
  df
  
}
