setwd("D:\\Project\\MubberryTE")
library(tidyverse)

# Fig1 --------------------------------------------------------------------
library(ggtree)
library(ggplot2)
tree <- treeio::read.newick( "summary_run_tree.newick")

tree_annot <- read.table("summary_run_annot.txt", sep = "\t", header = TRUE)
temp <- fortify(tree)
edge <- merge(temp, tree_annot, by.x = "label", sort = FALSE,  all.x = TRUE)
tree2 <- treeio::as.treedata(edge)

p <- revts(ggtree(tree)) + 
  theme_tree2() +
  geom_nodelab()+
  scale_x_continuous(labels = abs, limits = c(-60, 40)) + 
  geom_tiplab(hjust = -.2, align = TRUE)
p
pie_data <- edge[, c(3, 10, 13)]

# ordered by the letters of name
pies <- nodepie(pie_data, cols = c(2,3))
pies <- lapply(pies, function(f) f + scale_fill_manual(values = c('#99d594', '#fc8d59')))

p1 <- p + geom_inset(pies, width = 0.1, height = 0.1)
p2 <- p1 %<+% edge +  geom_label2(aes(subset=`Sig.Contractions` != 0,
                                      x = branch, 
                                      label = paste0("+", `Sig.Expansions`, 
                                                     "/-", `Sig.Contractions`)), 
                                  size = 3,
                                  nudge_y = 0.25,
                                  nudge_x = -0.25)

# adding legend manually
p3 <- p2 + annotate("point", x= -40, y =c(6.5, 7), 
                    shape=15, 
                    color=c('#99d594', '#fc8d59'), 
                    size=3) + 
  annotate("text", x= -35, y =c(6.5, 7), 
           label= c("Expansions", "Contractions"))

p3 
#ggsave('tree2.pdf', width = 9, height = 6, plot = p3)




## Chromosome --------------------------------------------------------------
library(tidyverse)
Chrom <- read.table("chrom_length.tsv")
Chrom$idx <- as.numeric( gsub("Chr", "", Chrom$V2))

Chrom$V1 <- factor(Chrom$V1, levels = c(
  "Cannabis_sativa", "Artocarpus_heterophyllus", "Morus_yunnanensis",
  "Morus_notabilis", "Morus_mongolica", "Morus_indica",
  "Morus_alba", "Morus_atropurpurea"
  
  
))

idx <- Chrom$V1 %in%  c("Cannabis_sativa", "Artocarpus_heterophyllus", "Morus_indica")

Chrom <- Chrom[!idx, ]

ggplot(Chrom ) + 
  geom_col(aes(x = V1, y = V3/1000000, group = as.factor(idx)),
           colour = "white",fill = "#fc8d59", size = 1.2) +
  #geom_segment(aes(x = idx, y = 0, xend = idx, yend = V3/1000000), size = 5, lineend = "round") + 
  labs(x = '', y = 'Chromosomal length(Mb)')  + 
  theme_bw() + 
  theme(text = element_text(size = 14)) +
  coord_flip() + 
  theme(
        axis.ticks.y = element_blank())


ggsave('S1.pdf', width = 9, height = 6)










# GO ----------------------------------------------------------------------

gene_idx <- c("ZZB_Hap", "EVM0", "Mi", "S1Chr", "MnotChr", "FeDS")
OrgDb_name <- c( "Malba", "Matropurpurea", "Mindica", "Mmongolica", "Mnotabilis", "Myunnanensis")
Annot_file <- c("AnnotMalba", "AnnotMatropurpurea")

sp = 1
df <- read.table("var_gene.tsv")
idx <- grep(gene_idx[sp], df$V1)
gene_list = df[idx, ]
gene_list <- gsub("\\.[0-9]$", "", gene_list)


df <- read.table("TE2gene4GO.tsv")

go_rich <- clusterProfiler::enrichGO(gene = df$V1,
                                     OrgDb = paste0("org.", OrgDb_name[sp], ".eg.db"), 
                                     ont = 'CC', #MF,CC,BP,ALL
                                     keyType = "GID",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.01,
                                     qvalueCutoff = 0.01)

my_go <- go_rich@result



GID_K <- read.table(paste0(Annot_file[sp], "/GID_K.txt"), header = T)
GID_K$GID <- gsub("\\.[0-9]$", "", GID_K$GID)

K_pathway <- read.table(paste0(Annot_file[sp], "/K_pathway.txt"), header = T)
tmp <- merge(GID_K, K_pathway, by = "K")

pathway_des <- read.table(paste0(Annot_file[sp], "/pathway_Description.txt"), header = T, sep = "\t")
tmp <- merge(tmp, pathway_des, by = "PATHWAY")

K_pathway[K_pathway$PATHWAY == 'ko00942', ]
GID2ko <-  dplyr::select(tmp, c("PATHWAY", "GID"))
ko2Term <- dplyr::select(tmp, c("PATHWAY", "description"))

K_rich <- clusterProfiler::enricher(gene = df$V1,
                                    TERM2GENE = GID2ko,
                                    TERM2NAME = ko2Term,
                                    pAdjustMethod = 'fdr',
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.01)

my_k <- K_rich@result

grep("An", my_k$Description)

identical(df$geneID[1], df$geneID[6])










# WGD ---------------------------------------------------------------------
library(aggregatesR)
dupgenes <- read.table("Wgd_kaks.tsv", sep = "\t", 
                       header = FALSE, check.names = FALSE)
names(dupgenes) <- c("Type", "Ka", "Ks","Ka_Ks")
dupgenes = na.omit(dupgenes)
idx <- grep("Indica", dupgenes$Type)
dupgenes <- dupgenes[-idx, ]

dupgenes$Ka <- as.numeric(dupgenes$Ka, na.rm = TRUE)
dupgenes$Ks<- as.numeric(dupgenes$Ks, na.rm = TRUE)
dupgenes$Ka_Ks <- as.numeric(dupgenes$Ka_Ks, na.rm = TRUE)


ggplot(dupgenes) + 
  stat_density(aes(x = Ks, col =  Type, group = Type), 
                           geom = "line",  
                           position = "identity", 
                           linewidth = 0.8, 
                           alpha = 0.8) +
  scale_color_manual(
    values = c("#e41a1c", "#377eb8",  "#4daf4a", 
               "#984ea3", "#ff7f00", "#ffff33")) +
  theme_bw() + 
  labs(x = "Ks", y = 'Density')  +
  theme(axis.ticks.x = element_blank(),
        text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.6)) +
  scale_x_continuous(limits = c(0, 6))




# ggsave(filename =  "Fig3_B.pdf", width = 6,
#       height = 6, device = "pdf")


# TE ----------------------------------------------------------------------

## TE Distribution -------------------------------------------------------
df <- readxl::read_xlsx("meta_data.xlsx", sheet = 3)
df$Percent <- df$Percent * 100
df <- subset(df, Notes == "Re-annotated" )
colnames(df)[1] <- "name"

df$group <- paste0(df$Order, "/", df$Superfamily)
df$name <- factor(df$name, levels = c("Morus yunnanensis",
  "Morus notabilis", "Morus mongolica", "Morus indica",
  "Morus alba", "Morus atropurpurea"
  
  
))

ggplot(df) + geom_col(aes(y = name, x = Percent, fill = group)) + 
  theme_bw() + 
  labs(x = "Genome coverage of TEs (%)", y = "") +
  theme(axis.ticks.y = element_blank(),
        text = element_text(size = 24))

ggsave(filename = "Fig4_a.pdf", width = 9, height = 6)



## Insertion time ----------------------------------------------------------

tmp <- read.table('LTR_insertion_time.tsv')
tmp$V10 <- as.numeric(tmp$V10)/1000000

7451 - 7464

density_maximum <- function(data) {
  
  max_idx <- which.max(density(data)$y)
  max_x <- density(data)$x[max_idx]
  
}

9/11
tmp %>% group_by(V1) %>% 
  summarise(t = density_maximum(V10)  )

#"Copia_LTR_retrotransposon" 
#"Gypsy_LTR_retrotransposon"
#"LTR_retrotransposon"
ggplot(tmp) + stat_density(aes(x = V10, col =  V1, group = V1), 
                           geom = "line",  
                           position = "identity", 
                           linewidth = 0.8, 
                           alpha = 0.8) +
  scale_color_manual(
    values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
  theme_bw() + 
  labs(x = "Insertion time (Mya)", y = 'Density')  +
  theme(axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8))  

#ggsave(filename = "Fig4_b.pdf", width = 9, height = 6)
# 
# library(aggregatesR)
# df <- read.table('TEanno_chrom_len.tsv', header = FALSE)
# df$V2 <- gsub("[A-Za-z]*", "", df$V2)
# df$V3 <- as.numeric(gsub("MB", "", df$V3))
# max(df$V5 )
# mydfreq <- NULL
# for (i in unique(df$V6)){
#   tmp <- filter(df, V6 == i)
#   
#   freqtable <- tapply(tmp$V4, tmp$V6, function(x) binFreqTable(x, step = 10000))
#   
#   for (id in attributes(freqtable)$dimnames[[1]]){
#     freqtable[[id]]$group <- id 
#   }
#   freqtable <- bind_rows(freqtable)
#   freqtable$range <- as.numeric(gsub(' - .*', "",   freqtable$range))
#   mydfreq <- rbind(mydfreq, freqtable)
#   }
#   #var <- paste0('p', i) 
# 
# mydfreq <- separate(mydfreq, group, sep = "-", into = c("Sp", "Chr"))
# 
# ggplot(subset(df, V1 == 'Morus_atropurpurea' & V2 == "08"), 
#        aes(x = V3, y = V4/V5, group = V2)) + 
#     #geom_point(aes(col = Chr), alpha = 0.5) + 
#     geom_smooth(aes(col = V2), se = FALSE) + 
#     theme_bw() +
#     theme(axis.ticks.x = element_blank(),
#           #axis.text.x = element_text(angle = 90, size = 14, hjust = 1, vjust = 1),
#           legend.text = element_text(size = 14),
#           axis.text = element_text(size = 14),
#           axis.title = element_text(size = 14),
#           legend.title = element_blank(),
#           legend.position = 'top') +
#     guides(color = guide_legend(nrow = 1)) + 
#     labs(x = "", y = "Number of TEs")
# 
# ggplot(subset(df, V1 == 'Morus_notabilis' & V2 == "6"), 
#        aes(x = V3, y = V4/V5, group = V2)) + 
#   #geom_point(aes(col = Chr), alpha = 0.5) + 
#   geom_smooth(aes(col = V2), se = FALSE) + 
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         #axis.text.x = element_text(angle = 90, size = 14, hjust = 1, vjust = 1),
#         legend.text = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         legend.title = element_blank(),
#         legend.position = 'top') +
#   guides(color = guide_legend(nrow = 1)) + 
#   labs(x = "", y = "Number of TEs")
# 
# 
# #ggsave(p, filename = paste0("TEanno_", chrom, "_all.tiff"), width = 9, height = 6)

## Half life -------------------------------------------------------------
LTR <- read.table('LTR_insertion_time.tsv')

half_time <- data.frame()
for (i in unique(LTR$V1)){
  df <- filter(LTR, V1 == i)
  df %>% mutate(Insertion_BIN = (V10 %/% 1000000) / 2) %>% 
    group_by(Insertion_BIN) %>% 
    summarise(Count=n()) -> df
  
  log.model <- lm(log2(Count) ~ Insertion_BIN, df)
  log.model.df <- data.frame(x = df$Insertion_BIN, y = 2^(fitted(log.model)))
  log.model$coefficients
  
  a <- log.model$coefficients[[1]]
  b <- log.model$coefficients[[2]]
  
  fitted.curve<-function(y){return((log2(y)-a)/b)}
  
  half.y <- (2^a)/2
  half.x <- fitted.curve(half.y)
  
  half_tmp <- data.frame(V1 = i,
                         half_time =half.x)
  
  half_time <- rbind(half_time, half_tmp )
  
  fy <- max(df$Count) / 4 *3
  hy <- max(df$Count) / 2
  # formula <- TeX(paste0(
  #   '$y = ', round(2^coef(log.model)[1]), '*',
  #   '2^{-', abs(round(coef(log.model)[2], 4)), 't}$'
  # ))
  # half.formula = TeX(paste0(
  #   '$t^{\\frac{1}{2}} = ', round(half.x, 2), 'Mys$'
  # ))
  
  
  # ggplot(df, aes(x = Insertion_BIN, y = Count)) +
  #    geom_bar(stat='identity', fill='gray') +
  #    geom_line(data = log.model.df, aes(x, y), size = 0.8, color='red') +
  #   # geom_text(aes(x = half.x, y = (2^a)/2 * 1.5, label = round(half.x, 2))) +
  #    # annotate(geom='text', x=5, y=fy, label = formula, size=5) +
  #    #  annotate(geom='text', x=6, y=hy, label = half.formula, size=5) +
  #    annotate(
  #      geom = "segment", x = half.x, xend = half.x,
  #      y = (2^a)/2, yend = -Inf, lty = "dashed")+
  #    annotate(geom = "segment", x = half.x,
  #             xend = -Inf, y = (2^a)/2, yend = (2^a)/2,
  #             lty="dashed")+
  #    theme_bw() +
  #    labs(x = 'Insertion Time (Mys)')
  #  
  #  ggsave(p, file = paste0('half_time_', i, ".pdf"), width = 9, height = 6)
}



half_life <- NULL
for (i in unique(LTR$V1)){
  df <- dplyr::filter(LTR, V1 == i )
  geom1 <- MASS::fitdistr(df$V10, "geometric")
  p1 <- geom1$estimate["prob"]
  deathRate <- p1/(1 - p1)
  
  
  half.x = log(2) /deathRate/1000000
  
  half_tmp <- data.frame(V1 = i,
                         half_time = half.x)
  
  half_life <- rbind(half_life, half_tmp )
  
}


df <- readxl::read_xlsx("total_meta.xlsx")
df <- na.omit(df)
df$V1 <- gsub(".v.*", "", df$Version)
df$V1 <- fct_reorder(df$V1, df$ConsensusTEs)

half_time <- merge(df, half_life,  by = 'V1', sort = FALSE) 
half_time$`AssembleSize(bp)`
cor.test(half_time$AssembleSize, half_time$half_time)

## TE2gene -----------------------------------------------------------------
file_list <- list.files(pattern = "*.gend2TE.bed")
mydf <- NULL
for (fl in file_list){
  tmp <- data.table::fread(fl, sep = "\t")
  tmp$sp <- gsub("\\.gend2TE.bed", "", fl)
  mydf <- rbind(tmp, mydf)
  
}

plot_data <- mydf
plot_data$V11 <- ifelse(plot_data$V11 == 0, "TE", "non-TE")
plot_data <-dplyr::select(plot_data, c("sp", "V5", "V6", "V11"))
plot_data <- unique(plot_data)
head(plot_data)
plot_data %>% 
  group_by(sp, V6, V11) %>% 
  summarise(tb = table(V6)) %>% 
  ungroup() %>% 
  group_by(sp, V6) %>% 
  mutate(percent = tb/sum(tb)) -> plot_data 

plot_data$V6 <- factor(plot_data$V6, 
                       levels = c("gene", "mRNA", 
                                  "three_prime_UTR", "five_prime_UTR",
                                  "CDS", "exon"))
ggplot(subset(plot_data, V11 == "TE")) +
  geom_boxplot(aes(x = V6, y = percent)) + 
  theme_bw() + 
  labs(x = "", y = "Percentage of TE insertions") +
  theme(text = element_text(size = 20))
#ggsave(filename = "Fig4_c2.pdf", width = 9, height = 6)

library(agricolae)
tmp <- as.data.frame(subset(plot_data, V11 == "TE"))
tmp$V6 <- as.character(tmp$V6)
tmp$percent <- as.numeric(tmp$percent)

model<- aov(percent~V6, data = tmp)
out <- HSD.test(model, "V6")

out$groups

####
library(aggregatesR)
plot_data <- subset(mydf, V6 == 'gene')
plot_data %>% 
  group_by(sp, V5) %>% 
  summarise(dis = sum(V11)) -> plot_data

#plot_data <- subset(plot_data, dis <= 5000)
freqtable <- binFreqTable(plot_data$dis, step = 100)
freqtable$range <- as.numeric(gsub(' -.*', "",   freqtable$range))

ggplot(freqtable) + 
  geom_col(aes(x = range, y = frequency)) + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(y = "Number of genes", x = "Distance between TEs and their nearest genes")

#ggsave( filename = "S2.pdf", width = 9, height = 6)

plot_data <- subset(mydf, V6 == 'gene')
plot_data %>% 
  group_by(sp, V5) %>% 
  summarise(dis = sum(V11)) -> plot_data

plot_data$dis <- ifelse(plot_data$dis == 0, 'TE', 'non-TE')


plot_data %>% 
  group_by(sp, dis) %>% 
  summarise(tb = n()) %>% 
  ungroup() %>% 
  group_by(sp) %>% 
  mutate(percent = tb/sum(tb))-> plot_data
plot_data %>% group_by(dis) %>% 
  summarise(mean(percent))

ggplot(plot_data) + geom_col(aes(y = sp, x = percent, fill = dis), position = "stack") + 
  theme_bw() + 
  labs(x = "Gene ratio (%)", y = "") +
  theme(axis.ticks.y = element_blank(),
        text = element_text(size = 14))

#ggsave( filename = "S3.pdf", width = 9, height = 6)

plot_data <- subset(mydf, V6 == 'gene')
plot_data <- select(plot_data, c("V5", "V11", "sp"))

core_gene <- read.table("gene_Core.tsv")
core_gene$type <- "core"
dispensable_gene <- read.table("gene_Dispensable.tsv")
dispensable_gene$type <- "dispensable"
private_gene <- read.table("gene_Private.tsv")
private_gene$type <- "private"
softcore_gene <- read.table("gene_SoftCore.tsv")
softcore_gene$type <- "softcore"

gene_type <- rbind(core_gene, dispensable_gene, private_gene, softcore_gene)
gene_type <- select(gene_type, c("V9", "type"))
colnames(gene_type)[1] <- "V5"
plot_data <- merge(plot_data, gene_type, by = "V5")
head(plot_data)
plot_data %>% 
  group_by(V5, sp, type) %>% 
  summarise(dis = sum(V11)) -> plot_data


plot_data$dis <- ifelse(plot_data$dis == 0, 'TE', 'non-TE')


plot_data %>% 
  group_by(sp, type, dis) %>% 
  summarise(tb = n())  %>% 
  ungroup() %>% 
  group_by(sp, type) %>% 
  mutate(percent = tb/sum(tb)) -> plot_data 

plot_data$sp <- factor(plot_data$sp, levels = c("Morus_yunnanensis",
  "Morus_notabilis", "Morus_mongolica", "Morus_indica",
  "Morus_alba", "Morus_atropurpurea"))

plot_data$type <- factor(plot_data$type, levels = c("core", "softcore", "dispensable", "private"))
  
ggplot(plot_data) + geom_col(aes(y = sp, x = percent, fill = dis), position = "stack") + 
  facet_grid(.~type) +
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  theme_bw() + 
  labs(x = "Gene ratio (%)", y = "") +
  theme(axis.ticks.y = element_blank(),
        text = element_text(size = 20),
        legend.position = "bottom") 

#ggsave( filename = "Fig4_d.pdf", width = 9, height = 6)


# Pan-genome --------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(agricolae)
df <- read.table("gene_stats.tsv")
df$V1 <- factor(df$V1, levels = c("Core", "SoftCore", "Dispensable", "Private"))

model<- aov(V2~V1, data = df)
out <- HSD.test(model, "V1")

out$groups


p1 <- ggplot(df) + geom_boxplot(aes(x = V1, y = V2, fill = V1), show.legend = FALSE) + 
  scale_y_continuous(limits = c(0, 12000), position = "right") +
  theme_bw() + 
  scale_fill_manual(values = c("Core" = "#e64b35", 
                               "Dispensable" = "#00a087", 
                               'Private' = "#3c5488",
                               "SoftCore" = "#4dbbd5")) +
  labs(y = "Length of genes", x = "") +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p1
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">


df <- data.frame(name = c("Core", "Dispensable", "Private", "SoftCore"),
                 value = c( 120311/6, 21461/6, 3031/6,  20585/6))
#(120311 + 20585)/(23964  + 25642 + 27435 + 21657 + 27414 +24689)
df$name <- factor(df$name, levels = c("Core", "SoftCore", "Dispensable", "Private"))
p2 <- ggplot(df) + geom_col(aes(x = name , y = value, fill = name), show.legend = FALSE) + 
  theme_bw() + 
  scale_fill_manual(values = c("Core" = "#e64b35", 
                               "Dispensable" = "#00a087", 
                               'Private' = "#3c5488",
                               "SoftCore" = "#4dbbd5")) +
  labs(y = "Number of genes", x = "")

p2 + p1
ggsave(p2 + p1, filename = "Fig2_de.pdf", width = 9, height = 6)

19423232  -  19420795
length(unique(df[df$V10 <= 0, ]$var))
length(unique(df$var))
max(df$V10)

df <- read.table()

## panGenome_CPC2 ------------------------------------------------------
df <- read.table("panGenome_CPC2.tsv", sep = "\t")

df %>% 
  group_by(V9, V8) %>% 
  summarise(tb = table(V8)) ->  plot_data

plot_data$V9 <- factor(plot_data$V9, levels = c("Core", "SoftCore", "Dispensable", "Private"))
plot_data$V8 <- factor(plot_data$V8, levels = c( "noncoding", "coding"))

p2 <- ggplot( plot_data) + 
  geom_col(aes(x = V9, y = tb, fill = V8), position = 'fill', show.legend = FALSE) +
  theme_bw() + 
  geom_hline(yintercept = 2100/(2100 + 3499)) +
  labs(y = "Coding and non-coding gene proportion", x = "") +
  theme(text = element_text(size = 18),
        #legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p2
ggplot(df) + 
  geom_boxplot(aes(x = V9, y = V7, fill = V8), show.legend = FALSE) 
  scale_y_continuous(limits = c(0, 2000)) +
  theme_bw() + 
  # scale_fill_manual(values = c("Core" = "#e64b35", 
  #                              "Dispensable" = "#00a087", 
  #                              'Private' = "#3c5488",
  #                              "SoftCore" = "#4dbbd5")) +
  labs(y = "Length of genes", x = "") +
  theme(text = element_text(size = 14),
        legend.position = "bottom")

# SV2gene ---------------------------------------------------------------
library(aggregatesR)
df <- read.table("var2gene.tsv", sep = "\t")
df <- subset(df, V6 == "gene")
df$var <-  paste0(df$V7, "_", df$V8, "_", df$V9)
df <- df[df$V10 < 5000, ]
tmp <- df[df$V10 ==0, ]
length(unique(tmp$var))
tmp <- df[df$V10 < 5000, ]
length(unique(tmp$var))

freqtable <- binFreqTable(df$V10, step = 100)
freqtable$range <- as.numeric(gsub(' -.*', "",   freqtable$range))
ggplot(freqtable) + 
  geom_col(aes(x = range, y = frequency)) + 
  theme_bw() +
  theme(text = element_text(size = 20))+
  labs(y = "Number of SVs", x = "The distance from the closest Genes")


#ggsave(filename = "Fig2_g.pdf", width = 9, height = 6)

df <- read.table("var2gene.tsv", sep = "\t")
df <- subset(df, V10 == 0)
df <- unique(df)
plot_data <- as.data.frame(table(df$V6))
plot_data
ggplot(plot_data, aes(x = "", y = Freq, fill = Var1)) + 
  geom_bar(width = 1, stat = "identity") +
  ggrepel::geom_text_repel(aes(x = 1.5, y = cumsum(Freq)-Freq/2, 
                label = paste0(Var1,"\n(", Freq, ")")),
                size = 6, 
                nudge_x = 0.5) +  
 # guides(fill = "none") +
  coord_polar("y") +
  theme_minimal() + 
  labs(x = "", y = "") + 
  theme(axis.text = element_blank())
#ggsave(filename = "Fig2_g2.pdf", width = 6, height = 6)


library(ggvenn)

df <- read.table("svim.tsv", header = TRUE)
unique(df$Morus_alba)
unique(df$Morus_atropurpurea)
unique(df$Morus_notabilis)
unique(df$Morus_yunnanensis)

df$Morus_alba <- ifelse(df$Morus_alba == "./.", FALSE, TRUE)
df$Morus_atropurpurea <- ifelse(df$Morus_atropurpurea == "./.", FALSE, TRUE)
df$Morus_notabilis <- ifelse(df$Morus_notabilis == "./.", FALSE, TRUE)
df$Morus_yunnanensis <- ifelse(df$Morus_yunnanensis == "./.", FALSE, TRUE)

df <- as_tibble(df)
ggplot(df) +
  geom_venn(aes(A = Morus_atropurpurea, 
                B =  Morus_alba,
                C = Morus_notabilis,
                D = Morus_yunnanensis),
                fill_color = c("#e64b35", "#00a087", "#3c5488", "#4dbbd5"),
                set_name_color =c("#e64b35", "#00a087", "#3c5488", "#4dbbd5")) +
  coord_fixed() +
  theme_void()

#ggsave(filename = "Fig2_f.pdf", width = 9, height = 6)




# Sv2gene structure -------------------------------------------------------
library(gggenes)

df <- read.table("atrop2alba_anthocyanin.gene.tsv")
df$V9 <- gsub(";Accession=.*", "", 
              gsub("\\.[1-9].*", "", 
                   gsub("ID=", "", df$V9)))
df$V4 <- ifelse(df$V7 == "+", 1*df$V4, -1 * df$V4)
df$V5 <- ifelse(df$V7 == "+", 1*df$V5, -1 * df$V5)


df2 <- read.table("atrop2alba_anthocyanin.tsv")
df2$V9 <- gsub(";Accession=.*", "", 
               gsub("\\.[1-9].*", "", 
                    gsub("ID=", "", df2$V9)))

df2$V11 <- ifelse(df2$V7 == "+", 1*df2$V11, -1 * df2$V11)
df2$V12 <- ifelse(df2$V7 == "+", 1*df2$V12, -1 * df2$V12)

df3 <-  read.table("atrop2alba_anthocyanin.TE.tsv")
df3$V5 <- gsub(";Accession=.*", "", 
              gsub("\\.[1-9].*", "", 
                   gsub("ID=", "", df3$V5)))
df3$V8 <- ifelse(df3$V4 == "+", 1*df3$V8, -1 * df3$V8)
df3$V9 <- ifelse(df3$V4 == "+", 1*df3$V9, -1 * df3$V9)

df3 <- select(df3, c("V5", "V8", "V9"))
colnames(df3) <- c("V9", "V13", "V14")

ggplot() + geom_rect(data = subset(df, V3 == 'gene'), aes(xmin=V4, 
                           xmax = V5, 
                           ymin = -0.01, 
                           ymax = 0.06),
                           fill = '#1b9e77', alpha = 0.3 ) +
  geom_segment(data = subset(df, V3 == 'mRNA'), aes(x = V4, 
                                                    xend = V5, 
                                                    y =  0.01, 
                                                    yend = 0.01), 
               col = '#d95f02', size = 2 ) +
  geom_segment(data = subset(df, V3 == 'CDS'), aes(x = V4, 
                                                   xend = V5, 
                                                   y =  0.02, 
                                                   yend =0.02), 
               col = '#7570b3', size = 2 ) +
  geom_segment(data = subset(df, V3 == 'exon'),  aes(x = V4, 
                                                     xend = V5, 
                                                     y =  0.03, 
                                                     yend = 0.03), 
               col = '#e6ab02' , size = 2) +
  geom_segment(data = subset(df, V3 == 'three_prime_UTR'),  aes(x = V4, 
                                                     xend = V5, 
                                                     y =  0.04, 
                                                     yend = 0.04), 
               col = '#4daf4a' , size = 2) +
  geom_segment(data = subset(df, V3 == 'five_prime_UTR'),  aes(x = V4, 
                                                     xend = V5, 
                                                     y =  0.05, 
                                                     yend = 0.05), 
               col = '#377eb8' , size = 2) +
  geom_rect(data = df2, aes(xmin = V11,
                            xmax = V12,
                            ymin = -0.01,
                            ymax = 0.03), fill = "black", alpha = 0.2) +
  geom_rect(data = df3, aes(xmin = V13,
                            xmax = V14,
                            ymin = 0.03,
                            ymax = 0.06), fill = "red", alpha = 0.2) +
  facet_wrap(.~V9, ncol = 1,  scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 1))  + 
  labs(x = "", y = "") 


# Anthocyanin content -----------------------------------------------------
library(agricolae)
df <- readxl::read_xlsx("meta_data.xlsx", sheet = 2)
model<- aov(Anthocyanin~Species, data = df)
out <- HSD.test(model, "Species")

out$groups
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(df, varname="Anthocyanin", 
                    groupnames=c("Species"))

df2
ggplot(df2) + geom_col(aes(x = Species, y = Anthocyanin), width = 0.2) +
  geom_errorbar(aes(x = Species,  ymin = Anthocyanin -sd, ymax = Anthocyanin + sd),
                width=.02,
                position=position_dodge(.9)) +
  theme_bw() +
  labs(x = "", y = "Anthocyanin content (mg/g fresh weight)") +
  theme(text = element_text(size = 14)) 

#ggsave(filename = "S4b.pdf", width = 3, height = 6)


# Syngap DEG -----------------------------------------------------------------
alba <- read.csv("alba_gene_count_matrix.csv")
rownames(alba) <- alba$gene_id

atrop <- read.csv("atrop_gene_count_matrix.csv")
rownames(atrop) <- atrop$gene_id
meta_data <- readxl::read_xlsx("meta_data.xlsx")

tmp <- subset(meta_data, BiologicalCondition == "DS")
atrop <- atrop[, match(tmp$ID, colnames(atrop)), ]
colnames(atrop) <- paste0("S", tmp$Time, "_",  tmp$Replicate)

tmp <- subset(meta_data, BiologicalCondition == "BS")
alba <-alba[, match(tmp$ID, colnames(alba)), ]
colnames(alba) <- paste0("S", tmp$Time, "_",  tmp$Replicate)

alba_cpm <- data.frame(edgeR::cpm(alba))
alba_cpm$GeneID <- rownames(alba)
atrop_cpm <- data.frame(edgeR::cpm(atrop))
atrop_cpm$GeneID <- rownames(atrop)

# write.table(atrop_cpm, file = "atrop_cpm.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(alba_cpm, file = "alba_cpm.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
gene_pair <- read.table("alba.atropurpurea.final.genepair", sep = "\t")
colnames(gene_pair) <- c("gene_id1", "gene_id2")
gene_pair$gene_id1 <- gsub("\\.[0-9]*$", "", gene_pair$gene_id1)
gene_pair$gene_id2 <- gsub("\\.[0-9]*", "", gene_pair$gene_id2)
gene_pair <- unique(gene_pair[c(1,2)])

# alba_cpm <- rename(alba_cpm, c( "gene_id1" = "GeneID" ))
# atrop_cpm <- rename(atrop_cpm, c( "gene_id2" = "GeneID" ))
# alba_cpm <- merge(gene_pair, alba_cpm, by = "gene_id1")
# merged_cpm <- merge(alba_cpm, atrop_cpm, by = "gene_id2")
364
1300/480 * 354
alba_cpm <- read.table( "alba_cpm.tsv", header = TRUE)
colnames(alba_cpm)[-10]  <- paste0("alba_", colnames(alba_cpm)[-10] )

atrop_cpm <- read.table( "atrop_cpm.tsv", header = TRUE)
colnames(atrop_cpm)[-10]  <- paste0("atrop_", colnames(atrop_cpm)[-10] )
syngap_alba <- read.table("syngap_alba_DEG.idx")
idx <- alba_cpm$gene_id1 %in% syngap_alba$V1
alba_cpm <- alba_cpm[idx, ]

idx <- gene_pair$gene_id1 %in% syngap_alba$V1
idx <- atrop_cpm$gene_id2 %in% gene_pair[idx, ]$gene_id2
atrop_cpm <- atrop_cpm[idx, ]

gene_pair <- unique(gene_pair[c(1,2)])
alba_cpm <- merge(alba_cpm, gene_pair, by = "gene_id1")
merged_cpm <- merge(alba_cpm, atrop_cpm, by = "gene_id2")
merged_cpm$ID <- paste0(merged_cpm$gene_id1, "=", merged_cpm$gene_id2)
merged_cpm %>% pivot_longer(cols = -c(, "gene_id2", "gene_id1", "ID")) -> plot_data
plot_data$name <- factor(plot_data$name, levels = c("alba_S4_1", "alba_S4_2", "alba_S4_3", 
                                                    "alba_S6_1", "alba_S6_2", "alba_S6_3", 
                                                    "alba_S10_1", "alba_S10_2", "alba_S10_3",
                                                    "atrop_S4_3", "atrop_S4_2", "atrop_S4_1", 
                                                    "atrop_S6_3", "atrop_S6_2", "atrop_S6_1",
                                                    "atrop_S10_3", "atrop_S10_2", "atrop_S10_1" ))

ggplot(plot_data) + geom_tile(aes(x =  name, y = ID, fill = log10(value))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradient2(midpoint = 2, 
                       low = '#313695',
                       mid = "white",
                       high = '#d73027') +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "top",
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "", y = "")
ggsave(filename = "test.pdf", width = 9, height = 9)

### GO Enrich ---------------------------------------------------------------
gene_idx <- c("ZZB_Hap", "EVM0", "Mi", "S1Chr", "MnotChr", "FeDS")
OrgDb_name <- c( "Malba", "Matropurpurea", "Mindica", "Mmongolica", "Mnotabilis", "Myunnanensis")
Annot_file <- c("AnnotMalba", "AnnotMatropurpurea")

sp = 1
gene_list <- atrop_cpm[idx, ]$gene_id1

top_num = 30
my_go <- NULL
keys <- c("MF", "CC", "BP")
for (key in keys){
  enrich <- clusterProfiler::enrichGO(gene = gene_list,
                                       OrgDb = paste0("org.", OrgDb_name[sp], ".eg.db"), 
                                       ont = key, #MF,CC,BP,ALL
                                       keyType = "GID",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.01,
                                       qvalueCutoff = 0.01)
  
  
  tmp = na.omit(enrich@result)
  tmp<- head(tmp, n = top_num)
  tmp$ONTOLOGY <- key
  my_go <- rbind(my_go, tmp)
  
}

write.csv(my_go, file = "my_GO4FC5.csv", row.names = FALSE, quote = FALSE)

go_rich  %>%
  group_by(ONTOLOGY) %>%
  slice_min(n = 8, order_by = p.adjust) -> mydf
mydf$Description

mydf$Description <- fct_reorder(mydf$Description, mydf$p.adjust)
ggplot(mydf) +
  geom_point(aes(x = -log10(p.adjust),
                 y = Description,
                 size = Count)) +
  theme_bw() +
  labs(x =  bquote(-log[10](p.adjust)), y = "") 
# #guides(size = "none") +
# scale_color_manual(values = col_platte) +
# theme(axis.text.y = element_text(colour = '#33a02c'))


#ggsave(filename = "Figures/S10.pdf", width = 9, height = 6)


col_platte <- c("BP" = '#33a02c',
                'MF' = '#e31a1c',
                'CC' = '#1f78b4')

col_text <- rep(c( '#33a02c', '#1f78b4', '#e31a1c'), each = top_num)


enrich_dat$Description <- gsub("oxidoreductase activity.*", 
                               "oxidoreductase activity",
                               enrich_dat$Description)
enrich_dat$Description <- fct_reorder(enrich_dat$Description, enrich_dat$ONTOLOGY)

ggplot(enrich_dat, 
       aes(x = Description, y = Count, fill = ONTOLOGY)) +
  geom_col() +
  #scale_fill_hue() +
  scale_fill_manual(values = col_platte) +
  theme_bw() +
  labs(y = "The number of genes") +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, 
                               colour = col_text, 
                               hjust = 1,
                               vjust = 0.5),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.position = c(0.8,0.9)) 

#ggsave(filename = "Figures/Fig5d.pdf", width = 9, height = 6)




### KEGG -----------------------------------------------------------------

GID_K <- read.table(paste0(Annot_file[sp], "/GID_K.txt"), header = T)
GID_K$GID <- gsub("\\.[0-9]$", "", GID_K$GID)

K_pathway <- read.table(paste0(Annot_file[sp], "/K_pathway.txt"), header = T)
tmp <- merge(GID_K, K_pathway, by = "K")

pathway_des <- read.table(paste0(Annot_file[sp], "/pathway_Description.txt"), header = T, sep = "\t")
tmp <- merge(tmp, pathway_des, by = "PATHWAY")

K_pathway[K_pathway$PATHWAY == 'ko00942', ]
GID2ko <-  dplyr::select(tmp, c("PATHWAY", "GID"))
ko2Term <- dplyr::select(tmp, c("PATHWAY", "description"))

df <- readxl::read_xlsx("alba.atropurpurea.final.clean.genepair.EVI.xlsx")
df$Rank <- as.numeric(rownames(df))
tmp <- head(df, 100)

library(tidyverse)

ggplot(df) + geom_point(aes(x = Rank, y = EVI)) +
  geom_point(data = tmp, aes(x = Rank, y = EVI), col = "red") +
  theme_bw() +
  theme(
    text = element_text(size = 24),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    #legend.direction = "horizontal",
    legend.position = c(0.8, 0.2)) +
  labs(x = "Samples ranked by EVI (high to low)", y = "expression variation index (EVI)")

#ggsave(filename = "Fig5b.pdf", width = 9, height = 6)

tmp <- head(df, 120)
K_rich <- clusterProfiler::enricher(gene = tmp$GeneID1,
                                    TERM2GENE = GID2ko,
                                    TERM2NAME = ko2Term,
                                    pAdjustMethod = 'fdr',
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.01)

my_k <- K_rich@result

my_k  <- head(my_k, n = top_num)

# write.csv(my_k, file = "syngap_top100.csv", row.names = FALSE, quote = FALSE)

my_k <- head(read.csv("syngap_top100.csv", sep = ","), n = 10)
my_k$pvalue <- as.numeric(my_k$pvalue)
my_k$Description <- fct_reorder(my_k$Description, desc(my_k$pvalue))
ggplot(my_k,
       aes(y = Description,
           x = -log2(pvalue),
           size = Count
           #fill = ONTOLOGY
       )) +
  geom_point() +
  labs(y = "") +
  #scale_y_continuous(limits = c(0, max(enrich_dat$Count * 1.1))) +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    #legend.direction = "horizontal",
    legend.position = c(0.8, 0.2)) +
  labs(x = "-log2(pvalue)")


#ggsave(filename = "Fig5c.pdf", width = 9, height = 6)
# GenePair -------------------------------------------------------------

#alba	23890 3950
#atropurpurea	25004 5064

df1 <- data.frame(   gene_num = c(953, 18987, 3950),
              type = c("2wayblast", "Synteny", "Unpaired"))

df1 <- arrange(df1, desc(type))

df2 <- data.frame(   gene_num = c(953, 18987, 5064),
                     type = c("2wayblast", "Synteny", "Unpaired"))

df2 <- arrange(df2, desc(type))


p1 <- ggplot(df1, aes(x = "", y = gene_num, fill = type))+
  geom_bar(width = 1.5, stat = "identity", show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(x = 1.5, y = cumsum(gene_num)-gene_num/2,
                               label = paste0(round(gene_num/sum(gene_num) * 100, 2), '%')),
                           size = 3,
                           nudge_x = 0.5,
  ) +
  coord_polar("y") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text = element_blank(),
        text = element_text(size = 14))
p2 <- ggplot(df2, aes(x = "", y = gene_num, fill = type))+
  geom_bar(width = 1.5, stat = "identity", show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(x = 1.5, y = cumsum(gene_num)-gene_num/2,
                               label = paste0(round(gene_num/sum(gene_num) * 100, 2), '%')),
                           size = 3,
                           nudge_x = 0.5,
  ) +
  coord_polar("y") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text = element_blank(),
        text = element_text(size = 14))

p1 + p2

#ggsave(filename = "Fig4a.pdf", width = 9, height = 6)



### PCP2Genepairs --------------------------------------------------------

atrop_gff <- read.table("Morus_atropurpurea.Chrom.gff")
atrop_gff$V9 <- gsub("\\..*", "", 
  gsub("ID=", "",
  gsub(";Parent=.*", "", 
  gsub(";Accession=.*", "", atrop_gff$V9))))


alba_gff <- read.table("Morus_alba.Chrom.gff")
alba_gff$V9 <- gsub("\\.[1-9].*", "", 
                     gsub("ID=", "",
                          gsub(";Parent=.*", "", 
                               gsub(";Accession=.*", "", alba_gff$V9))))

alba_gff$var_len <- alba_gff$V5 -  alba_gff$V4
atrop_gff$var_len <- atrop_gff$V5 -  atrop_gff$V4


alba_pcp <- read.delim2("Morus_alba.CPC2.txt", header = TRUE)
alba_pcp$ID <- gsub("\\.[1-9].*", "", alba_pcp$ID)
alba_pcp$coding_probability <- as.numeric(alba_pcp$coding_probability)

atrop_pcp <-  read.table("Morus_atropurpurea.CPC2.txt", header = TRUE)
atrop_pcp$ID <- gsub("\\.[1-9].*", "", atrop_pcp$ID)
atrop_pcp$coding_probability <- as.numeric(atrop_pcp$coding_probability)


alba_unpaired <- read.table("Unpaired_alba_gene.idx", header = TRUE)
atrop_unpaired <- read.table("Unpaired_atropurpurea_gene.idx", header = TRUE)

idx <-alba_pcp$ID  %in% alba_unpaired$gene_id
alba_unpaired_pcp <- alba_pcp[idx, ]
alba_paired_pcp <- alba_pcp[!idx, ]


idx <- atrop_pcp$ID  %in% atrop_unpaired$gene_id
atrop_unpaired_pcp <- atrop_pcp[idx, ]
atrop_paired_pcp <- atrop_pcp[!idx, ]





ggplot() + 
  geom_violin(data = alba_paired_pcp,  
               aes(y = peptide_length, x = "A"), 
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_violin(data = alba_unpaired_pcp,  
            aes(y = peptide_length, x = "B"), 
            draw_quantiles = c(0.25, 0.5, 0.75)) +
 geom_violin(data = atrop_paired_pcp,  
            aes(y = peptide_length, x = "C"), 
            draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_violin(data = atrop_unpaired_pcp,  
              aes(y = peptide_length, x = "D"), 
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  theme_bw() +
  theme(text = element_text(size = 20))


#ggsave(filename = "S6_2.pdf", width = 9, height = 6)

## qPCR -----------------------------------------------------------------
df1 <- read.csv("PCR/Alba.csv")
ref_g1 <- subset(df1, PreID == "R1")
mean_ct1 <- mean(ref_g1$Cq)
df1$d_ct1 <- df1$Cq - mean_ct1



df2 <- read.csv("PCR/Atropurpurea.csv")
ref_g2 <- subset(df2, PreID == "R1")
mean_ct2 <- mean(ref_g2$Cq)
df2$d_ct <- df2$Cq - mean_ct2


df <- data.frame(df1[, 3])
colnames(df) <- "PreID"

df$d_ct_ct <- df2$d_ct - df1$d_ct
df$fc <- 2^( -1 * df$d_ct_ct)

pre_df <- readxl::read_xlsx("PCR/primer.xlsx")
df <- merge(df, pre_df, by = "PreID")

df %>% group_by(PreID, GeneID) %>% 
  summarise(fc_m = paste0(mean(fc), "+" , sd(fc))) %>% 
  arrange( desc(GeneID)) 

my_qPCR <- df
my_qPCR$fc_m <- round(my_qPCR$fc_m, 2)
unique(my_qPCR$GeneID)

## DEG2SV ------------------------------------------------------------------

my_gene_ids <- unique(my_qPCR$GeneID)[!unique(my_qPCR$GeneID) %in% c( "Action")]



for (gene_id in unique(my_gene_ids)) {
## gene pair
gene_pair <- read.table("alba.atropurpurea.final.genepair", sep = "\t")
colnames(gene_pair) <- c("gene_id1", "gene_id2")
gene_pair$gene_id1 <- gsub("\\.[0-9]*$", "", gene_pair$gene_id1)
gene_pair$gene_id2 <- gsub("\\.[0-9]*", "", gene_pair$gene_id2)
gene_pair <- unique(gene_pair[c(1,2)])
gene_pair$comb_gene <- paste0(gene_pair$gene_id1, "_", gene_pair$gene_id2)

tmp <- read.table("atropurpurea.mongolica.final.genepair", sep = "\t")
colnames(tmp) <- c("gene_id2", "gene_id3")
tmp$gene_id2 <- gsub("\\.[0-9]*", "", tmp$gene_id2)
tmp$gene_id3 <- gsub("\\.t[1-9]*", "", tmp$gene_id3)
gene_pair <- merge(gene_pair, tmp, all.x = TRUE, by = 'gene_id2')

gene_pos <- grep(gene_id, gene_pair$gene_id1)

gene2sv <- read.table("qPCR_gene2variants.tsv")
gene2sv$V5 <- gsub(";.*", "", gsub("ID=", "", gene2sv$V5))
gene2sv <- subset(gene2sv, V11 != -1)
gene2sv$V8 <- gene2sv$V8 - gene2sv$V2
gene2sv$V9 <- gene2sv$V9 - gene2sv$V2

alba_gene2TE <- read.table("qPCR_alba2TE.tsv")
alba_gene2TE$V5 <- gsub(";.*", "", gsub("ID=", "", alba_gene2TE$V5))
alba_gene2TE$V10 <- gsub("ID=", "", alba_gene2TE$V10)
alba_gene2TE$V8 <- alba_gene2TE$V8 - alba_gene2TE$V2
alba_gene2TE$V9 <-alba_gene2TE$V9 - alba_gene2TE$V2


atrop_gene2TE <- read.table("qPCR_atrop2TE.tsv")
atrop_gene2TE$V5 <- gsub(";.*", "", gsub("ID=", "", atrop_gene2TE$V5))
atrop_gene2TE$V10 <- gsub("ID=", "", atrop_gene2TE$V10)
atrop_gene2TE$V8 <- atrop_gene2TE$V8 - atrop_gene2TE$V2
atrop_gene2TE$V9 <- atrop_gene2TE$V9 - atrop_gene2TE$V2


mongolica_gene2TE <- read.table("qPCR_mongolica2TE.tsv")
mongolica_gene2TE$V5 <- gsub(";.*", "", gsub("ID=", "", mongolica_gene2TE$V5))
mongolica_gene2TE$V10 <- gsub("ID=", "", mongolica_gene2TE$V10)
mongolica_gene2TE$V8 <- mongolica_gene2TE$V8 - mongolica_gene2TE$V2
mongolica_gene2TE$V9 <- mongolica_gene2TE$V9 - mongolica_gene2TE$V2


#### blat data
blat_dat <- read.table("syngap_alba2atropurpurea_blast.tsv")

colnames(blat_dat) <- c("Query", "Subject", "identity", "length", "mismatches", "gap", 
                        "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score")

#blat_dat <- subset(blat_dat, identity > 90)
tmp1 <- separate(blat_dat, col = "Query", into = c("Chrom1", "s_e", "gene_id2"), sep = "_")
tmp2 <- separate(tmp1 , col = "s_e", into = c("atrop_start", "atrop_end"), sep = "-")
tmp1 <- separate(tmp2 , col = "Subject", into = c("s_e", "gene_id1"), sep = "_ZZB_")
tmp1$gene_id1 <- paste0("ZZB_", tmp1$gene_id1)
blat_dat <- separate(tmp1 , col = "s_e", into = c("Chrom2", "alba_start", "alba_end"))
blat_dat$comb_gene <- paste0(blat_dat$gene_id1, "_", blat_dat$gene_id2)


idx <- blat_dat$comb_gene %in% gene_pair$comb_gene
blat_dat <- blat_dat[idx, ]
blat_dat$atrop_start <- as.numeric(blat_dat$atrop_start)
blat_dat$atrop_end <- as.numeric(blat_dat$atrop_end)
blat_dat$alba_start <- as.numeric(blat_dat$alba_start)
blat_dat$alba_end <- as.numeric(blat_dat$alba_end)

blat_dat$atrop_end <- blat_dat$atrop_end - blat_dat$atrop_start + 5000
blat_dat$atrop_start <- blat_dat$atrop_start - blat_dat$atrop_start - 5000

blat_dat$alba_end <- blat_dat$alba_end - blat_dat$alba_start + 5000
blat_dat$alba_start <- blat_dat$alba_start - blat_dat$alba_start - 5000

blat_dat$q_end <- blat_dat$q_end - 5000
blat_dat$q_start <-  blat_dat$q_start - 5000

blat_dat$s_end <- blat_dat$s_end - 5000
blat_dat$s_start <-  blat_dat$s_start - 5000


#### gff 
alba_gff <- read.table("syngap_alba_DEG.gff")
atrop_gff <- read.table("syngap_atropurpurea_DEG.gff")
mongolica_gff <- read.table("syngap_mongolica_DEG.gff")


gene_pair$gene_id1[gene_pos ]
#### plot 
y_pos = 0.5
idx <- grep(gene_pair$gene_id1[gene_pos ], alba_gff$V5)
alba_gff <- alba_gff[idx, ]

idx <- grep(gene_pair$gene_id2[gene_pos ], atrop_gff$V5)
atrop_gff <- atrop_gff[idx, ]

idx  <- grep(gene_pair$gene_id1[gene_pos ], blat_dat$gene_id1) 
blat_dat <- blat_dat[idx, ]


blat_dat$atrop_start <- as.numeric(blat_dat$atrop_start) - as.numeric(filter(atrop_gff, V4 == "gene")$V2)
blat_dat$atrop_end <- as.numeric(blat_dat$atrop_end) - as.numeric(filter(atrop_gff, V4 == "gene")$V2)

blat_dat$alba_start <- as.numeric(blat_dat$alba_start) - as.numeric(filter(alba_gff, V4 == "gene")$V2)
blat_dat$alba_end <- as.numeric(blat_dat$alba_end) - as.numeric(filter(alba_gff, V4 == "gene")$V2)


tmp <- select(blat_dat, c("q_start", "q_end", "s_start", "s_end"))
tmp$group <- paste0("g", rownames(tmp))
polygon_dat <- pivot_longer(tmp, cols = c(q_start, q_end, s_start, s_end))

polygon_dat$yaxis <- ifelse(polygon_dat$name == 's_start', y_pos, 
                            ifelse(polygon_dat$name== "s_end", y_pos, 0))


### gene info
idx <- grep(gene_pair$gene_id1[gene_pos ], alba_gff$V5)
alba_gff <- alba_gff[idx, ]
alba_gff$V3 <- alba_gff$V3 - subset(alba_gff, V4 == "gene")$V2
alba_gff$V2 <- alba_gff$V2 - subset(alba_gff, V4 == "gene")$V2


idx <- grep(gene_pair$gene_id2[gene_pos ], atrop_gff$V5)
atrop_gff<- atrop_gff[idx, ]
atrop_gff$V3 <- atrop_gff$V3 - subset(atrop_gff, V4 == "gene")$V2
atrop_gff$V2 <- atrop_gff$V2 - subset(atrop_gff, V4 == "gene")$V2


# idx <- grep(gene_pair$gene_id3[gene_pos ], mongolica_gff$V5)
# mongolica_gff<- mongolica_gff[idx, ]
# mongolica_gff$V3 <- mongolica_gff$V3 - subset(mongolica_gff, V4 == "gene")$V2
# mongolica_gff$V2 <- mongolica_gff$V2 - subset(mongolica_gff, V4 == "gene")$V2

### TE insertion
idx <- grep(gene_pair$gene_id1[gene_pos ], alba_gene2TE$V5)
alba_gene2TE <- alba_gene2TE[idx, ]

idx <- grep(gene_pair$gene_id2[gene_pos ], atrop_gene2TE$V5)
atrop_gene2TE <- atrop_gene2TE[idx, ]

idx <- grep(gene_pair$gene_id3[gene_pos ],mongolica_gene2TE$V5)
mongolica_gene2TE <- mongolica_gene2TE[idx, ]

### gene2sv
idx <- grep(gene_pair$gene_id1[gene_pos ], gene2sv$V5)
gene2sv <- gene2sv[idx, ]
gene2sv <- subset(gene2sv, V6 == "gene")

col_pat <- c("exon" = "#d73027",  "CDS" = "#fdae61",
             "five_prime_UTR" = "#006837", "three_prime_UTR" = "#006837")

## gene2sv
idx <- grep(gene_pair$gene_id1[gene_pos ], gene2sv$V5)
gene2sv <- gene2sv[idx, ]



## plot
ggplot() + 
  ## blast
  geom_polygon(data = polygon_dat, aes(x = value,
                                       y = yaxis,
                                       group = group,
                                       fill = group),
               alpha = 0.5, show.legend = FALSE) +
  geom_segment(data = blat_dat, aes(x = s_start, xend = s_end, 
                                    y = y_pos, yend = y_pos), size = 2) +
  geom_segment(data = blat_dat, aes(x = q_start, xend = q_end, 
                                    y = 0, yend = 0), size = 2) + 
  ## gene body
  geom_segment(data = subset(atrop_gff, V4 == "gene"), 
               aes( x = V2, xend = V3, 
                    y = 0, yend = 0), size =12, col = "grey60") +
  geom_segment(data = subset(atrop_gff, V4 == "exon"), 
               aes( x = V2, xend = V3, 
                    y = 0, yend =0, col = V4), size =8) +
  geom_segment(data = subset(atrop_gff, V4 == "CDS"), 
               aes( x = V2, xend = V3, 
                    y = 0, yend = 0, col = V4), size =6) +
  geom_segment(data = subset(atrop_gff, V4 == "five_prime_UTR"), 
               aes( x = V2, xend = V3, 
                    y = 0, yend = 0, col = V4), size = 4) +
  geom_segment(data = subset(atrop_gff, V4 == "three_prime_UTR"), 
               aes( x = V2, xend = V3, 
                    y = 0, yend = 0, col = V4), size = 4) +
  
  
  geom_segment(data = subset(alba_gff, V4 == "gene"), 
               aes( x = V2, xend = V3, 
                    y = y_pos, yend = y_pos), size =12, col = "grey60") +
  geom_segment(data = subset(alba_gff, V4 == "exon"), 
               aes( x = V2, xend = V3, 
                    y = y_pos, yend = y_pos, col = V4), size =8) +
  geom_segment(data = subset(alba_gff, V4 == "CDS"), 
               aes( x = V2, xend = V3, 
                    y = y_pos, yend = y_pos, col = V4), size =6) +
  geom_segment(data = subset(alba_gff, V4 == "five_prime_UTR"), 
               aes( x = V2, xend = V3, 
                    y = y_pos, yend = y_pos, col = V4), size = 4) +
  geom_segment(data = subset(alba_gff, V4 == "three_prime_UTR"), 
               aes( x = V2, xend = V3, 
                    y = y_pos, yend = y_pos, col = V4), size = 4) +
  
  # geom_segment(data = subset(mongolica_gff, V4 == "gene"), 
  #              aes( x = V2, xend = V3, 
  #                   y = -1, yend = -1), size =12, col = "grey60") +
  # geom_segment(data = subset(mongolica_gff, V4 == "exon"), 
  #              aes( x = V2, xend = V3, 
  #                   y = -1, yend = -1, col = V4), size =8) +
  # geom_segment(data = subset(mongolica_gff, V4 == "CDS"), 
  #              aes( x = V2, xend = V3, 
  #                   y = -1, yend = -1, col = V4), size =6) +
  # geom_segment(data = subset(mongolica_gff, V4 == "five_prime_UTR"), 
  #              aes( x = V2, xend = V3, 
  #                   y = -1, yend = -1, col = V4), size = 4) +
  # geom_segment(data = subset(mongolica_gff, V4 == "three_prime_UTR"), 
  #              aes( x = V2, xend = V3, 
  #                   y = -1, yend = -1, col = V4), size = 4) +
  
 
  ## TE
  geom_segment(data =  alba_gene2TE,
               aes( x = V8, xend = V9,  y = y_pos + 0.1, 
                    yend = y_pos + 0.1), size =4) +
  ggrepel::geom_text_repel(data = alba_gene2TE, 
                           aes( x = (V8 + (V9 - V8)/2), y = y_pos + 0.1, label = V11) ) +
  geom_segment(data =  atrop_gene2TE,
               aes( x = V8, xend = V9,  y = -0.1, 
                    yend = -0.1), size =4) +
  ggrepel::geom_text_repel(data = atrop_gene2TE, 
                             aes( x = (V8 + (V9 - V8)/2), y = -0.1, label = V11) ) +
  # geom_segment(data =  mongolica_gene2TE,
  #              aes( x = V8, xend = V9,  y = -1.1, 
  #                   yend = -1.1), size =4) +
  
  ## SV
  geom_segment(data =  gene2sv,
               aes( x = V8, xend = V9,  y = -0.2, 
                    yend = -0.2), size =4) +
  geom_text(aes(x = -5000, y = -0.2, label = "SV")) +
  geom_text(aes(x = -5000, y = -0.1, label = "TE")) +
  # geom_text(aes(x = -5000, y = 1.2, label = "SV")) +
  geom_text(aes(x = -5000, y = y_pos + 0.1, label = "TE")) +
  # geom_text(aes(x = -5000, y = -1.1, label = "TE")) +
  theme_bw() + 
  scale_color_manual(values = col_pat) +
  ylab(label = "") + 
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.ticks.y = element_blank()) + 
  guides(fill = "none", col = "none") + 
  scale_y_continuous(breaks = c(0, y_pos), 
                     labels = c("Atrop", "Alba"))


#ggsave(filename = paste0("Fig5_", gene_id, ".pdf"), width = 9, height = 6)
}


