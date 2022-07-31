#######Functions for manual mosaic plot############

source("workflow_functions2.R")
AT_genes_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_intersect", sep = "/")
  y <- paste(spp, "AT_genes.gff3", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}

GC_genes_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_intersect", sep = "/")
  y <- paste(spp, "GC_genes.gff3", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}



AT_order_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_intersect", sep = "/")
  y <- paste(spp, "AT_order.gff3", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}

GC_order_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_intersect", sep = "/")
  y <- paste(spp, "GC_order.gff3", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}

AT_TE_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_intersect", sep = "/")
  y <- paste(spp, "AT_tclass.gff3", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}

GC_TE_path <- function(spp){
  x <- paste(spp, "results", "BEDTOOLS_intersect", sep = "/")
  y <- paste(spp, "GC_tclass.gff3", sep = "_")
  z <- paste(x, y, sep = "/")
  z
}

repeat_bed_name <- function(spp, order){

  y = paste(spp, order, sep = "_")
  z = paste(y, ".bed", sep = "")
 
  z
}



################################################################################
make_intersect_df <- function(AT_genes, GC_genes, AT_order, GC_order){
  AT_gene <- read.table(AT_genes, sep = "\t")
  GC_gene <- read.table(GC_genes, sep = "\t")
  AT_order <- read.table(AT_order, sep = "\t")
  GC_order <- read.table(GC_order, sep = "\t")
  
  gene_intersect <- rbind.data.frame(AT_gene, GC_gene)
  gene_intersect <- data.frame("compartment" = gene_intersect$V3,
                               "start" = gene_intersect$V4, 
                               "end" = gene_intersect$V5,
                               "component" = gene_intersect$V11,
                               "tstart" = gene_intersect$V12,
                               "tend" = gene_intersect$V13)
  order_intersect <- rbind.data.frame(AT_order, GC_order)
  order_intersect <- data.frame("compartment" = order_intersect$V3,
                                "start" = order_intersect$V4, 
                                "end" = order_intersect$V5,
                                "tname" = order_intersect$V13,
                                "tstart" = order_intersect$V11,
                                "tend" = order_intersect$V12)
  consensi <- consensi_names()
  order_intersect <- left_join(order_intersect, consensi, by = "tname")
  order_intersect <- group_all_tclass(order_intersect) 
  order_intersect <- mutate(order_intersect, 
                            order = sapply(str_split(order_intersect$group, "/"), "[[", 1))  
  order_intersect <- order_intersect[,c("compartment", "start", "end", "component", "tstart", "tend")]

  intersect <- rbind.data.frame(order_intersect, gene_intersect)
  intersect <- intersect %>% mutate(compartment_len = end - start, 
                                    component_len = tend - tstart)
  
  component_len <- intersect %>% group_by(compartment, compartment_len, component) %>% summarise(component_len_total = sum(component_len))
  component_count <- intersect %>% group_by(compartment, compartment_len, component) %>% count(component)
  new_df <- left_join(component_len, component_count)
  new_df
}


################################################################################

component_len_summary <- function(AT_genes, GC_genes, AT_order, GC_order){
  df <- make_intersect_df(AT_genes, GC_genes, AT_order, GC_order)
  split_df <- split(df, df$compartment) #split by GC or AT
  AT <- split_df[[1]] 
  GC <- split_df[[2]]
  
  
  AT_rich <- AT %>% group_by(compartment, component) %>% 
    summarise(sum_len = sum(component_len_total))
  GC_rich <- GC %>% group_by(compartment, component) %>% 
    summarise(sum_len = sum(component_len_total))
  df <- list("AT_rich" = AT_rich, "GC_rich" = GC_rich)
  df
} #works til here




################################################################################

component_len_raw <- function(spp, AT_genes, GC_genes, AT_order, GC_order){
  df <- component_len_summary(AT_genes, GC_genes, AT_order, GC_order)
  AT <- df[[1]] 
  GC <- df[[2]] 
  
  df <- data.frame("Unknown_AT" = subset(AT, 
                                         AT$component == "Unknown")[["sum_len"]],
                   "LINE_AT" = subset(AT, 
                                      AT$component == "LINE")[["sum_len"]],
                   "LTR_AT" = subset(AT, 
                                     AT$component == "LTR")[["sum_len"]],
                   "DNA_AT" = subset(AT, 
                                     AT$component == "DNA")[["sum_len"]],
                   "MITE_AT" = subset(AT, 
                                      AT$component == "MITE")[["sum_len"]],
                   "Gene_AT" = subset(AT, 
                                      AT$component == "gene")[["sum_len"]],
                   
                   
                   ##GC
                   "Unknown_GC" = subset(GC, 
                                         GC$component == "Unknown")[["sum_len"]],
                   "LINE_GC" = subset(GC, 
                                      GC$component == "LINE")[["sum_len"]],
                   "LTR_GC" = subset(GC, 
                                     GC$component == "LTR")[["sum_len"]],
                   "DNA_GC" = subset(GC, 
                                     GC$component ==  "DNA")[["sum_len"]],
                   "MITE_GC" = subset(GC, 
                                      GC$component == "MITE")[["sum_len"]],
                   "Gene_GC" = subset(GC, GC$component == "gene")[["sum_len"]])
  
  
  
  total_AT <- AT_rich(AT_path(spp), ref_path(spp))[["AT_bp"]]
  total_GC <- ref_length(ref_path(spp)) - total_AT
  
  df2 <- data.frame("total_AT" = total_AT,
                    "anno_AT" = sum(AT$sum_len),
                    "unanno_AT" = total_AT - sum(AT$sum_len), 
                    
                    "total_GC" = total_GC,
                    "anno_GC" =  sum(GC$sum_len), 
                    "unanno_GC" = total_GC - sum(GC$sum_len))
  df <- data.frame(df, df2)
  df
}



################################################################################
mosaic_proportions <- function(spp, AT_genes, GC_genes, AT_order, GC_order){
  component_df <- component_len_raw(spp, AT_genes, GC_genes, AT_order, GC_order)
  t_df <- as.data.frame(t(component_df))
  colnames(t_df) <- "sum_len"
  t_df <- rownames_to_column(t_df, var = "component")
  components <- t_df[1:12,] 
  AT_comp <- components[1:6,]
  GC_comp <- components[7:12,]
  
  AT_props <- AT_comp %>% mutate(prop = (sum_len/sum(sum_len)*100))
  GC_props <- GC_comp %>% mutate(prop = (sum_len/sum(sum_len)*100))
  base <- t_df[13:18,]
  base <- base %>% mutate(prop = 
                            (sum_len)/(subset(base, base$component=="total_AT")[["sum_len"]] + 
                                         subset(base, base$component=="total_GC")[["sum_len"]])*100) 
  proportions <- rbind.data.frame(AT_props, GC_props, base)
  rownames(proportions) <- NULL
  proportions <- column_to_rownames(proportions[,c(1,3)], var = "component")
  colnames(proportions) <- NULL
  props <- as.data.frame(t(proportions))
  props
}

TE_proportions <- function(spp, AT_TE, GC_TE){

  AT_TE <- read.table(AT_TE, sep = "\t")
  GC_TE <- read.table(GC_TE, sep = "\t")
  
  AT_int <- data.frame("compartment" = AT_TE$V3,
                       "start" = AT_TE$V4, 
                       "end" = AT_TE$V5,
                       "tname" = AT_TE$V12,
                       "tstart" = AT_TE$V10,
                       "tend" = AT_TE$V11) 
  
  GC_int <- data.frame("compartment" = GC_TE$V3,
                       "start" = GC_TE$V4, 
                       "end" = GC_TE$V5,
                       "tname" = GC_TE$V12,
                       "tstart" = GC_TE$V10,
                       "tend" = GC_TE$V11) 
  
  
  consensi <- consensi_names()
  AT_int <- left_join(AT_int, consensi, by = "tname")
  GC_int <- left_join(GC_int, consensi, by = "tname")
  AT_int <- group_all_tclass(AT_int)
  GC_int <- group_all_tclass(GC_int)

  
  
  GC_int$group <- paste(GC_int$group, "GC", sep = "_") #unique ID from AT
  
  GC_int$group <- paste(GC_int$group, "class", sep = "_") #need to discern it from order entries later on. 
  AT_int$group <- paste(AT_int$group, "class", sep = "_")
  
  

  
  int <- rbind.data.frame(AT_int, GC_int)
  int <- rename(int, component = group)
  int <- int[,c("compartment", "start", "end", "component", "tstart", "tend")]

  
  
  TE <- int %>% mutate(compartment_len = end - start, 
                       component_len = tend - tstart)
  TE_component_len <- TE %>% group_by(compartment, compartment_len, component) %>% summarise(component_len_total = sum(component_len))
  TE_component_count <- TE %>% group_by(compartment, compartment_len, component) %>% count(component)
  TE_summary <- left_join(TE_component_len, TE_component_count)
  TE_lens <- TE_summary %>% group_by(compartment, component) %>% 
    summarise(sum_len = sum(component_len_total)) 
  
  

  AT <- TE_lens %>% filter(str_detect(compartment, "AT"))
  GC <- TE_lens %>% filter(!str_detect(compartment, "AT"))

  
  DNA <- filter(AT, str_detect(AT$component, "DNA"))
  LINE <- filter(AT, str_detect(AT$component, "LINE"))
  LTR <- filter(AT, str_detect(AT$component, "LTR"))
  Unknown <- filter(AT, str_detect(AT$component, "^Unknown")) #NOTE: Must have '^' or will pick up on things like "LTR_unknown and double up. 
  MITE <- filter(AT, str_detect(AT$component, "MITE"))
  
  
  DNAg <- filter(GC, str_detect(GC$component, "DNA"))
  LINEg <- filter(GC, str_detect(GC$component, "LINE"))
  LTRg <- filter(GC, str_detect(GC$component, "LTR"))
  Unknowng <- filter(GC, str_detect(GC$component, "^Unknown"))
  MITEg <- filter(GC, str_detect(GC$component, "MITE"))
  
  a1 <- DNA %>% mutate(prop = (sum_len/sum(sum_len)))
  a2 <- LINE %>% mutate(prop = (sum_len/sum(sum_len)))
  a3 <- LTR %>% mutate(prop = (sum_len/sum(sum_len)))
  a4 <- Unknown %>% mutate(prop = (sum_len/sum(sum_len)))
  a5 <- MITE %>% mutate(prop = (sum_len/sum(sum_len)))
  
  g1 <- DNAg %>% mutate(prop = (sum_len/sum(sum_len)))
  g2 <- LINEg %>% mutate(prop = (sum_len/sum(sum_len)))
  g3 <- LTRg %>% mutate(prop = (sum_len/sum(sum_len)))
  g4 <- Unknowng %>% mutate(prop = (sum_len/sum(sum_len)))
  g5 <- MITEg %>% mutate(prop = (sum_len/sum(sum_len)))
  

  TE_p <- rbind.data.frame(a1, a2, a3, a4, a5,
                           g1, g2, g3, g4, g5)
  TE_p
  
  
  
  #issue starts here
  TE_p <- TE_p[,c(2,4)]
  TE_p <- column_to_rownames(TE_p, var = "component")
  colnames(TE_p) <- NULL
  TE_p <- as.data.frame(t(TE_p))
  colnames(TE_p) <- str_replace(colnames(TE_p), "/", ".") 
  colnames(TE_p) <- str_replace(colnames(TE_p), "-", ".")
  TE_p
  
}  

mosaic_base <- function(df){
  
  base <- ggplot() + geom_rect(data = df, 
                               aes(xmin = -0.1,
                                   xmax = 100.2,
                                   ymin = -0.1,
                                   ymax = 100.2),
                               colour = "white", fill = "white", size = 0.1) +
    theme_minimal() + theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                            axis.text.y = element_blank()) #keeping border black for now
  
  #adding all AT
  #adding all AT
  b <- geom_rect(data = df,
                 aes(xmin = 0,
                     xmax = total_AT,
                     ymin = 0,
                     ymax = 100),
                 colour = "white", 
                 fill = "white", size = .1)
  
  #adding GC
  d <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = 100,
                     ymin = 0,
                     ymax = 100),
                 colour = "white", 
                 fill = "white", size = .1)
  
  #adding annotated AT
  e <- geom_rect(data = df,
                 aes(xmin = 0,
                     xmax = anno_AT,
                     ymin = 0,
                     ymax = 100),
                 colour = "white", 
                 fill = "white", size = .1)
  
  #adding unannotated AT 
  
  f <- geom_rect(data = df,
                 aes(xmin = anno_AT,
                     xmax = anno_AT + unanno_AT,
                     ymin = 0,
                     ymax = 100),
                 colour = "white", 
                 fill = "#ebebeb", size = .1)

  
  #adding unnano GC
  g <- geom_rect(data = df,
                 aes(xmin = 100 - unanno_GC,
                     xmax = 100,
                     ymin = 0,
                     ymax = 100),
                 colour = "white", 
                 fill = "#ebebeb", size = .1)
  
  #adding annotated GC
  h <- geom_rect(data = df,
                 aes(xmin = total_AT,#adding 0.5 to keep black margin thicc
                     xmax = total_AT+anno_GC,
                     ymin = 0,
                     ymax = 100),
                 colour = "white", 
                 fill = "white", size = .1)

  
  base <- base+b+d+e+f+g+h
  base
  
  
  #Unknown
  i <- geom_rect(data = df,
                 aes(xmin = 0.1,
                     xmax = anno_AT,
                     ymin = 0,
                     ymax = Unknown_AT),
                 colour = "white", 
                 fill = "#888290", size = .1)
  
  #MITE
  j <- geom_rect(data = df,
                 aes(xmin = 0.1,
                     xmax = anno_AT,
                     ymin = Unknown_AT,
                     ymax = MITE_AT + Unknown_AT),
                 colour = "white", 
                 fill = "#e4d2ad", size = .1)
  
  #LTR
  k <- geom_rect(data = df,
                 aes(xmin = 0.1,
                     xmax = anno_AT,
                     ymin = MITE_AT  + Unknown_AT,
                     ymax = MITE_AT +  Unknown_AT + LTR_AT),
                 colour = "white", 
                 fill = "#CCBBCD", size = .1)
  
  
  #LINE
  l <- geom_rect(data = df,
                 aes(xmin = 0.1,
                     xmax = anno_AT,
                     ymin = MITE_AT  + Unknown_AT + LTR_AT ,
                     ymax = MITE_AT  + Unknown_AT + LTR_AT  + LINE_AT),
                 colour = "white", 
                 fill = "#5f7872", size = .1)
  
  
  #DNA
  m <- geom_rect(data = df,
                 aes(xmin = 0.1,
                     xmax = anno_AT -0.1,
                     ymin = MITE_AT  + Unknown_AT + LTR_AT  + LINE_AT  ,
                     ymax = MITE_AT  + Unknown_AT + LTR_AT  + LINE_AT  + DNA_AT),
                 colour = "#647588", 
                 fill = "#647588", size = .1)
  
  #genes
  n <- geom_rect(data = df,
                 aes(xmin = 0.1, 
                     xmax = anno_AT - 0.1,
                     ymin = MITE_AT  + Unknown_AT + LTR_AT  + LINE_AT  + DNA_AT ,
                     ymax = MITE_AT  + Unknown_AT + LTR_AT  + LINE_AT  + DNA_AT  + Gene_AT),
                 colour = "#A6B0BB", 
                 fill = "#bbbdb5", size = .1)
  

  
  
  
  #############GC
  
  #Unknown
  o <- geom_rect(data = df,
                 aes(xmin = total_AT,#adding 0.5 to keep black margin
                     xmax = total_AT + anno_GC,
                     ymin = 0,
                     ymax = Unknown_GC),
                 colour = "white", 
                 fill = "#888290", size = .1)
  
  #MITE
  p <- geom_rect(data = df,
                 aes(xmin = total_AT,#adding 0.5 to keep black margin 
                     xmax = total_AT + anno_GC,
                     ymin = Unknown_GC,
                     ymax = MITE_GC + Unknown_GC),
                 colour = "white", 
                 fill = "#e4d2ad", size = .1)
  
  #LTR
  q <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = total_AT + anno_GC,
                     ymin = MITE_GC  + Unknown_GC,
                     ymax = MITE_GC  + Unknown_GC + LTR_GC),
                 colour = "white", 
                 fill = "#CCBBCD", size = .1)
  
  
  #LINE
  r <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = total_AT + anno_GC,
                     ymin = MITE_GC  + Unknown_GC + LTR_GC ,
                     ymax = MITE_GC  + Unknown_GC + LTR_GC + LINE_GC),
                 colour = "white", 
                 fill = "#5f7872", size = .1)
  
  
  #DNA
  s <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = total_AT + anno_GC,
                     ymin = MITE_GC  + Unknown_GC + LTR_GC + LINE_GC,
                     ymax = MITE_GC  + Unknown_GC + LTR_GC + LINE_GC + DNA_GC),
                 colour = "white", 
                 fill = "#647588", size = .1)
  
  
  #genes
  t <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = total_AT + anno_GC,
                     ymin = MITE_GC  + Unknown_GC + LTR_GC + LINE_GC + DNA_GC ,
                     ymax = MITE_GC  + Unknown_GC + LTR_GC + LINE_GC + DNA_GC  + Gene_GC),
                 colour = "white", 
                 fill = "#bbbdb5", size = .1)
  

  base_order <- base + i + j + k + l + m + n + o + p + q + r + s + t
  base_order
}

mosaic_base <- function(df){
  base <- ggplot() + geom_rect(data = df, 
                               aes(xmin = -0.1,
                                   xmax = 100.2,
                                   ymin = -0.1,
                                   ymax = 100.2),
                               colour = "black", fill = "white", size = 0.1) +
    theme_minimal() + theme(panel.grid = element_blank(),axis.text.x = element_blank(),
                            axis.text.y = element_blank())

  #adding all AT
  b <- geom_rect(data = df,
                 aes(xmin = 0,
                     xmax = total_AT,
                     ymin = 0,
                     ymax = 100),
                 colour = "black", 
                 fill = "white", size = .1)
  
  #adding GC
  d <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = 100,
                     ymin = 0,
                     ymax = 100),
                 colour = "black", 
                 fill = "white", size = .1)
  
  #adding annotated AT
  e <- geom_rect(data = df,
                 aes(xmin = 0,
                     xmax = anno_AT,
                     ymin = 0,
                     ymax = 100),
                 colour = "black", 
                 fill = "white", size = .1)
  
  #adding unannotated AT 
  
  f <- geom_rect(data = df,
                 aes(xmin = anno_AT,
                     xmax = anno_AT + unanno_AT,
                     ymin = 0,
                     ymax = 100),
                 colour = "black", 
                 fill = "#ebebeb", size = .1)
  
  ##these parts are slightly redundant but it adds the lines
  
  #adding unnano GC
  g <- geom_rect(data = df,
                 aes(xmin = 100 - unanno_GC,
                     xmax = 100,
                     ymin = 0,
                     ymax = 100),
                 colour = "black", 
                 fill = "#ebebeb", size = .1)
  
  #adding annotated GC
  h <- geom_rect(data = df,
                 aes(xmin = total_AT,
                     xmax = total_AT+anno_GC,
                     ymin = 0,
                     ymax = 100),
                 colour = "black", 
                 fill = "white", size = .1)

  
  base <- base+b+d+e+f+g+h
  base
  
}

mosaic_plot <- function(df){
  mosaic <- mosaic_base(df) 
  mosaic_plot <- mosaic + 
    geom_rect(data = df,
              aes(xmin = 0,
                  xmax = (LTR_class*anno_AT),
                  ymin = Unknown_AT + MITE_AT,
                  ymax = Unknown_AT + MITE_AT + LTR_AT),
              color = "white",
              fill = "#9f6c71", size = 0.1) +
    
    geom_rect(data = df,
              aes(xmin = (LTR_class*anno_AT),
                  xmax = (LTR_class+LTR.Copia_class)*anno_AT,
                  ymin = Unknown_AT + MITE_AT,
                  ymax = Unknown_AT + MITE_AT + LTR_AT),
              colour = "white", 
              fill = "#c39aa2", size = .1) +

    
    geom_rect(data = df,
              aes(xmin = (LTR_class+LTR.Copia_class)*anno_AT,
                  xmax = (LTR_class+LTR.Copia_class+LTR.Gypsy_class)*anno_AT,
                  ymin = Unknown_AT  + MITE_AT,
                  ymax = Unknown_AT  + MITE_AT + LTR_AT),
              colour = "white", 
              fill = "#e4cbcf", size = .1) +
    ##GC_LTRs
    geom_rect(data = df,
              aes(xmin = total_AT,
                  xmax = total_AT + (anno_GC * LTR_GC_class),
                  ymin = Unknown_GC + MITE_GC,
                  ymax = Unknown_GC + MITE_GC + LTR_GC),
              colour = "white", 
              fill = "#9f6c71", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = total_AT + (anno_GC*LTR_GC_class),
                  xmax = ((LTR_GC_class+LTR.Copia_GC_class)*anno_GC) + total_AT,
                  ymin = Unknown_GC + MITE_GC,
                  ymax = Unknown_GC + LTR_GC + MITE_GC),
              colour = "white", 
              fill = "#c39aa2", size = .1)+
    
    geom_rect(data = df,
              aes(xmin = ((LTR_GC_class+LTR.Copia_GC_class)*anno_GC) + total_AT,
                  xmax = ((LTR_GC_class+LTR.Copia_GC_class+LTR.Gypsy_GC_class)*anno_GC) + total_AT,
                  ymin = Unknown_GC + MITE_GC,
                  ymax = Unknown_GC + MITE_GC + LTR_GC),
              colour = "white", 
              fill = "#e4cbcf", size = .1) +
    

    
    
    ##LINE AT
    geom_rect(data = df,
              aes(xmin = anno_AT,
                  xmax = LINE_class*anno_AT,
                  ymin = Unknown_AT + MITE_AT + LTR_AT,
                  ymax = Unknown_AT + MITE_AT + LTR_AT + LINE_AT),
              colour = "white", 
              fill = "#5f7872", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = LINE_class *anno_AT,
                  xmax = (LINE.Tad1_class + LINE_class)*anno_AT,
                  ymin = Unknown_AT + MITE_AT + LTR_AT,
                  ymax = Unknown_AT + MITE_AT + LTR_AT + LINE_AT),
              colour = "white", 
              fill = "#abcdbf", size = .1) +
    
    ##LINE GC
    
    geom_rect(data = df,
              aes(xmin = total_AT,
                  xmax = total_AT+ (LINE_GC_class*anno_GC),
                  ymin = Unknown_GC + MITE_GC + LTR_GC  ,
                  ymax = Unknown_GC + MITE_GC + LTR_GC + LINE_GC),
              colour = "white", 
              fill = "#5f7872", size = .1) +
    
    ##AT DNA
    geom_rect(data = df,
              aes(xmin = 0,
                  xmax = DNA_class * anno_AT,
                  ymin = Unknown_AT   + MITE_AT+ LTR_AT + LINE_AT ,
                  ymax = Unknown_AT   + MITE_AT+ LTR_AT + LINE_AT + DNA_AT),
              colour = "white",
              fill = "#272950", size = .1) +
    
    geom_rect(data = df,
              aes(xmin =  DNA_class *anno_AT,
                  xmax = (DNA_class + DNA.hAT_class)*anno_AT,
                  ymin = Unknown_AT   + MITE_AT+ LTR_AT + LINE_AT,
                  ymax = Unknown_AT   + MITE_AT+ LTR_AT + LINE_AT + DNA_AT),
              colour = "white", 
              fill = "#3f4174", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = (DNA_class + DNA.hAT_class)*anno_AT,
                  xmax = (DNA_class + DNA.hAT_class + DNA.MuLE_class )*anno_AT,
                  ymin = Unknown_AT   + MITE_AT + LTR_AT + LINE_AT,
                  ymax = Unknown_AT   + MITE_AT + LTR_AT + LINE_AT + DNA_AT),
              colour = "white", 
              fill = "#6a6ba3", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = (DNA_class + DNA.hAT_class + DNA.MuLE_class )*anno_AT,
                  xmax = (DNA_class + DNA.hAT_class + DNA.MuLE_class + DNA.PIF.Harbinger_class)*anno_AT,
                  ymin = Unknown_AT + MITE_AT + LTR_AT + LINE_AT,
                  ymax = Unknown_AT + MITE_AT + LTR_AT + LINE_AT + DNA_AT),
              colour = "white", 
              fill = "#969ac9", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = (DNA_class + DNA.hAT_class + DNA.MuLE_class + DNA.PIF.Harbinger_class)*anno_AT,
                  xmax = (DNA_class + DNA.hAT_class + DNA.MuLE_class + DNA.PIF.Harbinger_class + DNA.Tc1.Mariner_class)*anno_AT,
                  ymin = Unknown_AT + MITE_AT + LTR_AT + LINE_AT,
                  ymax = Unknown_AT + MITE_AT + LTR_AT + LINE_AT + DNA_AT),
              colour = "white", 
              fill = "#a1c5df", size = .1) +
    
    #GC DNA
    
    geom_rect(data = df,
              aes(xmin = total_AT,
                  xmax = total_AT+ (DNA_GC_class*anno_GC),
                  ymin = Unknown_GC + MITE_GC + LTR_GC + LINE_GC ,
                  ymax = Unknown_GC + MITE_GC + LTR_GC + LINE_GC + DNA_GC),
              colour = "white", 
              fill = "#272950", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = total_AT+ (DNA_GC_class*anno_GC),
                  xmax = total_AT+ (DNA_GC_class + DNA.hAT_GC_class)*anno_GC,
                  ymin = Unknown_GC + MITE_GC + LTR_GC + LINE_GC ,
                  ymax = Unknown_GC + MITE_GC + LTR_GC + LINE_GC + DNA_GC),
              colour = "white", 
              fill = "#3f4174", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = total_AT+ (DNA_GC_class + DNA.hAT_GC_class)*anno_GC,
                  xmax = total_AT+ (DNA_GC_class + DNA.hAT_GC_class + DNA.MuLE_GC_class)*anno_GC,
                  ymin = Unknown_GC + MITE_GC + LTR_GC + LINE_GC ,
                  ymax = Unknown_GC + MITE_GC + LTR_GC + LINE_GC + DNA_GC),
              colour = "white", 
              fill = "#6a6ba3", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = total_AT+(DNA_GC_class + DNA.hAT_GC_class + DNA.MuLE_GC_class)*anno_GC,
                  xmax = total_AT+(DNA_GC_class + DNA.hAT_GC_class + DNA.MuLE_GC_class + DNA.PIF.Harbinger_GC_class)*anno_GC,
                  ymin = Unknown_GC + MITE_GC + LTR_GC + LINE_GC ,
                  ymax = Unknown_GC + MITE_GC + LTR_GC + LINE_GC + DNA_GC),
              colour = "white", 
              fill = "#969ac9", size = .1) +
    
    geom_rect(data = df,
              aes(xmin = total_AT+(DNA_GC_class + DNA.hAT_GC_class + DNA.MuLE_GC_class + DNA.PIF.Harbinger_GC_class)*anno_GC,
                  xmax = total_AT+(DNA_GC_class + DNA.hAT_GC_class + DNA.MuLE_GC_class + DNA.PIF.Harbinger_GC_class + DNA.Tc1.Mariner_GC_class)*anno_GC,
                  ymin = Unknown_GC + MITE_GC + LTR_GC + LINE_GC ,
                  ymax = Unknown_GC + MITE_GC + LTR_GC + LINE_GC + DNA_GC),
              colour = "white", 
              fill = "#a1c5df", size = .1) 
  mosaic_plot
  
}

