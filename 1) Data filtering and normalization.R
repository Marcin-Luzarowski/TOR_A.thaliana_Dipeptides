#set working directory
setwd("your working directory")

#load libraries
library(tidyverse)

### normalize all samples (columns) to its median factor (median/median of all)
medianNormalize <- function(x)
{
  # dataNew <- x
  x[x == 0] <- NA
  dataNew <- sweep(x, 2, apply(x, 2, median, na.rm = TRUE)/median(as.matrix(x), na.rm = TRUE), "/") ### Take data set and DIVIDE each value in particular column by medFactor calculated above
  dataNew <- as.data.frame(dataNew) ### Save dataset and a data.frame
}

###########################LOAD DATA#####################

d1 <- read.delim("pp_peaks.txt", check.names = FALSE)
annotation <- read.delim("annotation.txt")

rownames(d1) <- d1$Peak.ID

d1_order <- as.data.frame(read.delim("order.txt"))

for_order <- as.character(d1_order$Name)

d1_ordered <- d1[for_order]


###################################################################################
### Filtering

#filter intensity 10000 and save only features which appear in 31 samples
freq_d1 <- rowSums(d1_ordered > 10000)

df1 <- cbind(d1_ordered, freq_d1)

df1 <- df1[df1$freq_d1 > 30,]
df1 <- df1[,1:105]

###################################################################################
### Normalization

#Normalization to median
df1_n <- medianNormalize(df1)

##################################################################################
### Normalization

#Pick annotated metabolites

df1_a <- rownames_to_column(df1, "Peak.ID")
df1_a <- inner_join(df1_a, annotation, by = "Peak.ID")

df1_n_a <- rownames_to_column(df1_n, "Peak.ID")
df1_n_a <- inner_join(df1_n_a, annotation, by = "Peak.ID")

rownames(df1_a) <- df1_a$Entry_ID
rownames(df1_n_a) <- df1_n_a$Entry_ID

df1_a <- df1_a[,2:106]
df1_n_a <- df1_n_a[,2:106]

colnames(df1_n_a) <- as.character(d1_order$Sample_complex)

#save normalized data
write.table(df1_n_a, "normalized_annotated.txt", sep = "\t")


