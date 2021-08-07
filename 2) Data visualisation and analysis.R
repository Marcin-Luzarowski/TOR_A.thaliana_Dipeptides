####GENERATE BASIC STATISTICAL INFORMATION ABOUT THE SAMPLES#######

library(tidyverse)

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
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


setwd("your working directory")

#Prepare SD data

df_sd <- read.delim("normalized_annotated_short_day.txt")

#read sample list file
order_sd <- read.delim("order_sd.txt")

df_sd_long <- gather(df_sd, "Sample", "Intensity", 2:ncol(df_sd))

df_sd_long <- separate(df_sd_long, Sample, sep = "_", into = c("Mutant", "Time", "Replica"))
df_sd_long$Mutant <- factor(df_sd_long$Mutant, levels = c("WT", "rb5"))
df_sd_long$Time <- factor(df_sd_long$Time, levels = c("0h", "1h", "2h", "4h", "7h", "8h", "12h", "16h", "20h", "23h", "24h"))

#Prepare LD data
df_ld <- read.delim("normalized_annotated_long_day.txt")
order_ld <- read.delim("order_ld.txt")

df_ld_long <- gather(df_ld, "Sample", "Intensity", 2:ncol(df_ld))

df_ld_long <- separate(df_ld_long, Sample, sep = "_", into = c("Mutant", "Time", "Replica"))
df_ld_long$Mutant <- factor(df_ld_long$Mutant, levels = c("WT", "rb5"))
df_ld_long$Time <- factor(df_ld_long$Time, levels = c("0h", "1h", "2h", "4h", "8h", "12h", "15h", "16h", "20h", "23h", "24h"))

df_sd_long_mean_sd <- data_summary(df_sd_long, varname = "Intensity", groupnames = c("Metabolite", "Mutant", "Time"))
write.table(df_sd_long_mean_sd, file = "short_day_mean_sd.txt", sep = "\t")

df_ld_long_mean_sd <- data_summary(df_ld_long, varname = "Intensity", groupnames = c("Metabolite", "Mutant", "Time"))
write.table(df_ld_long_mean_sd, file = "long_day_mean_sd.txt", sep = "\t")

####GENERATE LINE PLOTS SHOWING SINGLE DATA POINTS FOR CHOSEN DIPEPTIDES####
library(tidyverse)
library(ggsci)
library(reshape2)
library(ggdendro)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(scales)

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
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

setwd("your working directory")
#Prepare SD data

df_sd <- read.delim("normalized_annotated_short_day.txt")
order_sd <- read.delim("order_sd.txt")

df_sd_long <- gather(df_sd, "Sample", "Intensity", 2:ncol(df_sd))

df_sd_long <- separate(df_sd_long, Sample, sep = "_", into = c("Mutant", "Time", "Replica"))
df_sd_long$Mutant <- factor(df_sd_long$Mutant, levels = c("WT", "rb5"))
df_sd_long$Time <- factor(df_sd_long$Time, levels = c("0h", "1h", "2h", "4h", "7h", "8h", "12h", "16h", "20h", "23h", "24h"))

df_sd_long$Intensity <- log(df_sd_long$Intensity, base = 2)


#Prepare LD data
df_ld <- read.delim("normalized_annotated_long_day.txt")
order_ld <- read.delim("order_ld.txt")

df_ld_long <- gather(df_ld, "Sample", "Intensity", 2:ncol(df_ld))

df_ld_long <- separate(df_ld_long, Sample, sep = "_", into = c("Mutant", "Time", "Replica"))
df_ld_long$Mutant <- factor(df_ld_long$Mutant, levels = c("WT", "rb5"))
df_ld_long$Time <- factor(df_ld_long$Time, levels = c("0h", "1h", "2h", "4h", "8h", "12h", "15h", "16h", "20h", "23h", "24h"))

df_ld_long$Intensity <- log(df_ld_long$Intensity, base = 2)



#Prepare a list of metabolites 

dipeptides_to_plot <- read.delim("Dipeptides_to_plot_list.txt")
dipeptides_to_plot <- as.character(dipeptides_to_plot$Metabolite)

met_list <- dipeptides_to_plot # here you must include list of common metabolites!
i <- 1
plot_i <- 1
my_ggplot_list_barplot_sd <- vector('list')
my_ggplot_list_barplot_ld <- vector('list')

myColors <- c("Black", "Red")

while(i <= length(met_list)){
  df_sd_long_tmp <- filter(df_sd_long, Metabolite == met_list[i])
  df_sd_long_summary <- data_summary(df_sd_long_tmp, varname = "Intensity", groupnames = c("Mutant", "Time"))
  df_sd_long_summary$x_axis <- rep(c(0, 1, 2, 4, 7, 8, 12, 16, 20, 23, 24),2)
  
  df_ld_long_tmp <- filter(df_ld_long, Metabolite == met_list[i])
  df_ld_long_summary <- data_summary(df_ld_long_tmp, varname = "Intensity", groupnames = c("Mutant", "Time"))
  df_ld_long_summary$x_axis <- rep(c(0, 1, 2, 4, 8, 12, 15, 16, 20, 23, 24),2)
  
  #PREPARE TMP - SINGLE POINTS ON GRAPH SHORT DAY
  df_sd_long_tmp <- separate(data = df_sd_long_tmp, col = Time, into = "x_axis", sep = "h", remove = FALSE)
  df_sd_long_tmp$x_axis <- as.numeric(df_sd_long_tmp$x_axis)
  
  #PREPARE TMP - SINGLE POINTS ON GRAPH LONG DAY
  df_ld_long_tmp <- separate(data = df_ld_long_tmp, col = Time, into = "x_axis", sep = "h", remove = FALSE)
  df_ld_long_tmp$x_axis <- as.numeric(df_ld_long_tmp$x_axis)
  
  rects_sd <- data.frame(xstart = c(0,8), xend = c(8,24), Conditions = c("light", "darkness"))
  rects_sd$conditions <- factor(rects_sd$Conditions, levels = c("light", "darkness"))
  
  rects_ld <- data.frame(xstart = c(0,16), xend = c(16,24), Conditions = c("light", "darkness"))
  rects_ld$conditions <- factor(rects_sd$Conditions, levels = c("light", "darkness"))
  
  my_breaks_sd <- floor(pretty(seq(min(df_sd_long_summary$Intensity, na.rm = TRUE), (max(df_sd_long_summary$Intensity, na.rm = TRUE)+1)*1.1)))
  
  p <- ggplot()+
    geom_line(data = df_sd_long_summary, aes(x = x_axis, y = Intensity, group = Mutant, color = Mutant))+
    geom_point(data = df_sd_long_tmp, aes(x = x_axis, y = Intensity, group = Mutant, color = Mutant), size = 1)+
    scale_y_continuous(breaks = my_breaks_sd)+
    theme_classic()+
    scale_color_manual(values = myColors)+
    labs(y = "Intensity [log2]", x = "Time [h]", title = met_list[i])+
    theme(plot.title = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.position = "None",
          legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)))+
    geom_rect(data = rects_sd, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = Conditions), alpha = 0.2)+
    scale_fill_manual(values = c("grey80", "khaki1"))
  
  my_breaks_ld <- floor(pretty(seq(min(df_ld_long_summary$Intensity, na.rm = TRUE), (max(df_ld_long_summary$Intensity, na.rm = TRUE)+1)*1.1)))
  
  p2 <- ggplot()+
    geom_line(data = df_ld_long_summary, aes(x = x_axis, y = Intensity, group = Mutant, color = Mutant))+
    geom_point(data = df_ld_long_tmp, aes(x = x_axis, y = Intensity, group = Mutant, color = Mutant), size = 1)+
    scale_y_continuous(breaks = my_breaks_sd)+
    theme_classic()+
    scale_color_manual(values = myColors)+
    labs(y = "Intensity [log2]", x = "Time [h]", title = met_list[i])+
    theme(plot.title = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.position = "None",
          legend.direction = "horizontal",
          legend.background = element_rect(fill=alpha(0.4)))+
    geom_rect(data = rects_ld, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = Conditions), alpha = 0.2)+
    scale_fill_manual(values = c("grey80", "khaki1"))
  
  my_ggplot_list_barplot_sd[[plot_i]] <- p
  my_ggplot_list_barplot_ld[[plot_i]] <- p2
  
  i <- i+1
  plot_i <- plot_i+1
  
}

names(my_ggplot_list_barplot_sd) <- dipeptides_to_plot

sd_plot <- ggarrange(plotlist = my_ggplot_list_barplot_sd,
                     ncol = 1,
                     nrow = 3)

pdf("sd_line_plots_single_points_list.pdf", onefile = TRUE)
sd_plot
dev.off()

names(my_ggplot_list_barplot_ld) <- dipeptides_to_plot

ld_plot <- ggarrange(plotlist = my_ggplot_list_barplot_ld,
                     ncol = 1,
                     nrow = 3)

pdf("ld_line_plots_single_points_list.pdf", onefile = TRUE)
ld_plot
dev.off()

####PERFORM INTERACTOME ANALYSIS BASED ON PROMIS DATASET FROM A.THALIANA####
library(tidyverse)
library(ggsci)
library(reshape2)
library(ggdendro)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(pheatmap)

setwd("your working directory")

#read table from PROMIS dataset
cor_table <- read.delim("PROMIS_cor_table.txt", check.names = FALSE)

#read clustering results
wt_ld <- read.delim("wt_ld_kmeans_clustering_list.txt")
wt_ld <- rownames_to_column(wt_ld, "Metabolite")

wt_sd <- read.delim("wt_sd_kmeans_clustering_list.txt")
wt_sd <- rownames_to_column(wt_sd, "Metabolite")

tor_ld <- read.delim("tor_ld_kmeans_clustering_list.txt")
tor_ld <- rownames_to_column(tor_ld, "Metabolite")

tor_sd <- read.delim("tor_sd_kmeans_clustering_list.txt")
tor_sd <- rownames_to_column(tor_sd, "Metabolite")

jn_tables <- full_join(wt_ld, wt_sd, by = "Metabolite")
jn_tables <- full_join(jn_tables, tor_ld, by = "Metabolite")

jn_tables <- full_join(jn_tables, tor_sd, by = "Metabolite")

colnames(jn_tables) <- c("Metabolite", "WT_LongDay", "WT_ShortDay", "TOR_LongDay", "TOR_ShortDay")

affected_metabolites <- data.frame(jn_tables$Metabolite)
colnames(affected_metabolites) <- "Metabolite"
affected_dipeptides <- filter(affected_metabolites, grepl("-", Metabolite))
promis_metabolites <- data.frame(colnames(cor_table))
colnames(promis_metabolites) <- "Metabolite"
promis_metabolites <- separate(promis_metabolites, Metabolite, sep = "_", into = c("Metabolite", "Peak"))

common_metabolites <- inner_join(promis_metabolites, affected_dipeptides, by = "Metabolite")
number_of_dipeptides_in_both_datasets <- unique(common_metabolites$Metabolite)
common_metabolites <- paste(common_metabolites$Metabolite, common_metabolites$Peak, sep = "_")

interactions <- data.frame(TAIR=character(),
                           Dipeptides=character(),
                           stringsAsFactors=FALSE) 

i <- 1

while(i <= nrow(cor_table)){
  cor_table_tmp <- select(cor_table, all_of(common_metabolites))
  cor_table_tmp <- data.frame(t(cor_table_tmp[i,]))
  colnames(cor_table_tmp) <- "Correlation"
  cor_table_tmp <- cor_table_tmp %>%
    filter(Correlation >= 0.7)
  if(nrow(cor_table_tmp) >= 3){ #here one can change with how many dipeptides protein should co-elute
    cor_table_tmp <- rownames_to_column(cor_table_tmp, var = "Metabolite")
    cor_table_tmp <- separate(cor_table_tmp, Metabolite, sep = "_", into = c("Metabolite", "Peak"))
    dipeptides <- paste(cor_table_tmp$Metabolite, collapse = " | ")
    interactions_tmp <- data.frame(cor_table$TAIR[i], cor_table$Peak_Maximum[i], dipeptides)
    colnames(interactions_tmp) <- c("TAIR", "Peak_maximum", "Dipeptides")
    interactions <- rbind(interactions, interactions_tmp)
  }
  i <- i+1
}

unique_proteins <- unique(interactions$TAIR)
write.table(unique_proteins, "proteins_interacting_with_3dipeptides.txt", sep = "\t")
write.table(interactions, "putative_interactions_3dipeptides.txt", sep = "\t")



####GENERATION OF HEATMAP AND K-MEANS CLUSTERING####
library(tidyverse)
library(ggsci)
library(reshape2)
library(ggdendro)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(pheatmap)


setwd("your working directory")
######PREPARE SD DATA
df_sd <- read.delim("normalized_annotated_short_day.txt")
order_sd <- read.delim("order_sd.txt")

########Prepare LD data
df_ld <- read.delim("normalized_annotated_long_day.txt")
order_ld <- read.delim("order_ld.txt")


##############Prepare a list of metabolites identified in both data sets
list_ld <- data.frame(df_ld$Metabolite)
colnames(list_ld) <- "Metabolite"

list_sd <- data.frame(df_sd$Metabolite)
colnames(list_sd) <- "Metabolite"

list_of_metabolites <- inner_join(list_ld, list_sd, by = "Metabolite")

#####Short day data - just metabolites found in both datasets
sd_heatmap_df <- inner_join(list_of_metabolites, df_sd, by = "Metabolite")

wt_order_sd <- order_sd[1:51,]
sample_list <- wt_order_sd$Sample_complex

wt_sd_heatmap_df <- sd_heatmap_df[,sample_list]
rownames(wt_sd_heatmap_df) <- sd_heatmap_df$Metabolite

samp <- wt_order_sd
groups <- factor(samp$Sample, levels = unique(samp$Sample))

tmp_anova <- wt_sd_heatmap_df
tmp_anova[is.na(tmp_anova)] <- 0
tmp_anova <- log2(tmp_anova)
tmp_anova[tmp_anova == -Inf] <- 0
tmp_anova[tmp_anova == 0] <- NA

wt_sd_heatmap_df_initial <- wt_sd_heatmap_df
wt_sd_heatmap_df <- tmp_anova

#### P-values Select metabolites significantly affected during the course of experiment - anova was applied on log2 values
anova.2 <- function(x,y) anova(lm(x ~ y))$Pr[1]
p.values <- apply(wt_sd_heatmap_df, 1, anova.2, groups)
p.values.adj <- p.adjust(p.values, "bonferroni")

wt_sd_heatmap_df$p_val <- p.values.adj
wt_sd_heatmap_df <- dplyr::filter(wt_sd_heatmap_df, p_val <= 0.05)

save_p_vals <- data.frame("Metabolite" = rownames(wt_sd_heatmap_df), 
                          "p.val" = wt_sd_heatmap_df$p_val)

wt_sd_heatmap_df <- wt_sd_heatmap_df[,sample_list]

affected_metabolites <- data.frame(rownames(wt_sd_heatmap_df))
colnames(affected_metabolites) <- "Metabolite"
wt_sd_heatmap_df_initial <- rownames_to_column(wt_sd_heatmap_df_initial, "Metabolite")
wt_sd_heatmap_df <- inner_join(wt_sd_heatmap_df_initial, affected_metabolites, by = "Metabolite")
rownames(wt_sd_heatmap_df) <- wt_sd_heatmap_df$Metabolite
wt_sd_heatmap_df <- wt_sd_heatmap_df[,-1]

write.table(x = save_p_vals, file = "ANOVA_sign_met_WT_SD.txt", sep = "\t")

##########CALCULATE MEAN VALUES
wt_sd_heatmap_df[is.na(wt_sd_heatmap_df)] <- 0
tmpR <- apply(wt_sd_heatmap_df, 1, function(x) {tapply(unlist(x), groups, mean, na.rm=TRUE)})
tmpR <- as.data.frame(t(tmpR))

tmpR_scaled <- t(scale(t(tmpR)))

#Prepare heatmap for combined data
paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(tmpR_scaled), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(tmpR_scaled)/paletteLength, max(tmpR_scaled), length.out=floor(paletteLength/2)))

class_all <- read.delim("pheatmap_class_sd_wt.txt", row.names = 1)

#colours <- list(Treatment = c("WT_LD" = "grey20", 
#                             "rb5_LD" = "indianred",
#                            "WT_SD" = "grey40",
#                           "rb5_SD" = "indianred1"))

colours <- list(Conditions = c("WT_SD_light" = "yellow",
                               "WT_SD_darkness" = "black"))


#class_all <- list(c(rep("Control", 10), rep("Cold", 10), rep("Heat", 10)))
#General heatmap for all dipeptides
set.seed(1)
pheatmap(tmpR_scaled,
         annotation_col = class_all,
         annotation_legend = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize_col = 6,
         fontsize = 8,
         main = "",
         #width = 6,
         filename = "Heatmap_WT_Short_day.jpeg",
         #height = 8,
         cellheight = 8,
         cellwidth = 12,
         na_col = "grey",
         #breaks = breaksList,
         #color = colorRampPalette(colorpanel(256, "blue", "white", "red"))(length(breaksList)),
         breaks = myBreaks,
         color = myColor,
         treeheight_row = 0,
         annotation_colors = colours,
         kmeans_k = NA,
         cutree_rows = NA
)

#K-means included
set.seed(1)
heatmap_obj_kmeans <- pheatmap(tmpR_scaled,
                               annotation_col = class_all,
                               annotation_legend = FALSE,
                               cluster_rows = TRUE,
                               cluster_cols = FALSE,
                               clustering_distance_rows = "euclidean",
                               clustering_method = "complete",
                               show_rownames = TRUE,
                               show_colnames = FALSE,
                               fontsize_row = 8,
                               fontsize_col = 8,
                               fontsize = 8,
                               main = "",
                               #width = 6,
                               filename = "Heatmap_WT_Short_day_3k.jpeg",
                               #height = 8,
                               cellheight = 48,
                               cellwidth = 12,
                               na_col = "grey",
                               #breaks = breaksList,
                               #color = colorRampPalette(colorpanel(256, "blue", "white", "red"))(length(breaksList)),
                               breaks = myBreaks,
                               color = myColor,
                               treeheight_row = 0,
                               annotation_colors = colours,
                               kmeans_k = 3,
                               cutree_rows = NA
)

wt_sd_kmeans_clustering_list <- data.frame(heatmap_obj_kmeans$kmeans$cluster)
write.table(wt_sd_kmeans_clustering_list, "wt_sd_kmeans_clustering_list.txt", sep = "\t")

###############HERE STARTS ANALYSIS OF LONG DAY DATA

#####Long day data - just metabolites found in both datasets
ld_heatmap_df <- inner_join(list_of_metabolites, df_ld, by = "Metabolite")
wt_order_ld <- order_ld[1:55,]
sample_list <- wt_order_ld$Sample_complex

wt_ld_heatmap_df <- ld_heatmap_df[,sample_list]
rownames(wt_ld_heatmap_df) <- ld_heatmap_df$Metabolite

samp <- wt_order_ld
groups <- factor(samp$Sample, levels = unique(samp$Sample))

tmp_anova <- wt_ld_heatmap_df
tmp_anova[is.na(tmp_anova)] <- 0
tmp_anova <- log2(tmp_anova)
tmp_anova[tmp_anova == -Inf] <- 0
tmp_anova[tmp_anova == 0] <- NA

wt_ld_heatmap_df_initial <- wt_ld_heatmap_df
wt_ld_heatmap_df <- tmp_anova

#### P-values Select metabolites significantly affected during the course of experiment
anova.2 <- function(x,y) anova(lm(x ~ y))$Pr[1]
p.values <- apply(wt_ld_heatmap_df, 1, anova.2, groups)
p.values.adj <- p.adjust(p.values, "bonferroni")

wt_ld_heatmap_df$p_val <- p.values.adj
wt_ld_heatmap_df <- dplyr::filter(wt_ld_heatmap_df, p_val <= 0.05)

save_p_vals <- data.frame("Metabolite" = rownames(wt_ld_heatmap_df), 
                          "p.val" = wt_ld_heatmap_df$p_val)

wt_ld_heatmap_df <- wt_ld_heatmap_df[,sample_list]

affected_metabolites <- data.frame(rownames(wt_ld_heatmap_df))
colnames(affected_metabolites) <- "Metabolite"
wt_ld_heatmap_df_initial <- rownames_to_column(wt_ld_heatmap_df_initial, "Metabolite")
wt_ld_heatmap_df <- inner_join(wt_ld_heatmap_df_initial, affected_metabolites, by = "Metabolite")
rownames(wt_ld_heatmap_df) <- wt_ld_heatmap_df$Metabolite
wt_ld_heatmap_df <- wt_ld_heatmap_df[,-1]

write.table(x = save_p_vals, file = "ANOVA_sign_met_WT_LD.txt", sep = "\t")

##########calculate mean on initial data
wt_ld_heatmap_df[is.na(wt_ld_heatmap_df)] <- 0
tmpR <- apply(wt_ld_heatmap_df, 1, function(x) {tapply(unlist(x), groups, mean, na.rm=TRUE)})
tmpR <- as.data.frame(t(tmpR))

tmpR_scaled <- t(scale(t(tmpR)))

#Prepare heatmap for combined data
paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(tmpR_scaled), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(tmpR_scaled)/paletteLength, max(tmpR_scaled), length.out=floor(paletteLength/2)))

class_all <- read.delim("pheatmap_class_ld_wt.txt", row.names = 1)

#colours <- list(Treatment = c("WT_LD" = "grey20", 
#                             "rb5_LD" = "indianred",
#                            "WT_SD" = "grey40",
#                           "rb5_SD" = "indianred1"))

colours <- list(Conditions = c("WT_LD_light" = "yellow",
                               "WT_LD_darkness" = "black"))


#class_all <- list(c(rep("Control", 10), rep("Cold", 10), rep("Heat", 10)))
#General heatmap for all dipeptides
set.seed(1)
pheatmap(tmpR_scaled,
         annotation_col = class_all,
         annotation_legend = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize_col = 6,
         fontsize = 8,
         main = "",
         #width = 6,
         filename = "Heatmap_WT_Long_day.jpeg",
         #height = 8,
         cellheight = 8,
         cellwidth = 12,
         na_col = "grey",
         #breaks = breaksList,
         #color = colorRampPalette(colorpanel(256, "blue", "white", "red"))(length(breaksList)),
         breaks = myBreaks,
         color = myColor,
         treeheight_row = 0,
         annotation_colors = colours,
         kmeans_k = NA,
         cutree_rows = NA
)

#K-means included
set.seed(1)
heatmap_obj_kmeans <- pheatmap(tmpR_scaled,
                               annotation_col = class_all,
                               annotation_legend = FALSE,
                               cluster_rows = TRUE,
                               cluster_cols = FALSE,
                               clustering_distance_rows = "euclidean",
                               clustering_method = "complete",
                               show_rownames = TRUE,
                               show_colnames = FALSE,
                               fontsize_row = 8,
                               fontsize_col = 8,
                               fontsize = 8,
                               main = "",
                               #width = 6,
                               filename = "Heatmap_WT_Long_day_3k.jpeg",
                               #height = 8,
                               cellheight = 48,
                               cellwidth = 12,
                               na_col = "grey",
                               #breaks = breaksList,
                               #color = colorRampPalette(colorpanel(256, "blue", "white", "red"))(length(breaksList)),
                               breaks = myBreaks,
                               color = myColor,
                               treeheight_row = 0,
                               annotation_colors = colours,
                               kmeans_k = 3,
                               cutree_rows = NA
)


wt_ld_kmeans_clustering_list <- data.frame(heatmap_obj_kmeans$kmeans$cluster)
write.table(wt_ld_kmeans_clustering_list, "wt_ld_kmeans_clustering_list.txt", sep = "\t")
