#library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DEqMS)
library(fgsea)

#Yeast Triple-Knockout TMT Standard

#Set the working directory for the project
setwd("/home/trainee/Documents/TKO_R/")

#If the table has been exported as a tab-separated text file
all_proteins <- read.table("MS3_FAIMS_Proteins_PD24.txt", sep = '\t', header = TRUE)
#If the table is in Excel format
#all_proteins <- read_excel("MS3_FAIMS_Proteins_PD24.xlsx")

colnames(all_proteins)
dim(all_proteins)
#Filter in the Master proteins and filter out contaminants
#Remove the unused columns
all_proteins <- all_proteins %>%
  filter(!grepl("cont_",Accession)) %>%
  filter(Master == "IsMasterProtein") %>%
  select(3:38)
#Add gene names from the description
all_proteins$Gene <- apply(
  all_proteins, 1, function(x) {
    descr <- x["Description"]
    #OBS: format-dependent operations
    gn <- strsplit(descr, "GN=")[[1]][[2]]
    strsplit(gn, " ")[[1]][[1]]
  }
)
dim(all_proteins)

#Construct a function to rename the ratios
#and simultaneously put the ratio columns in a list
#OBS: pay attention to the particular format of  the column names

rename_ratios <- function(df) {
  #Initialize the vector of ratio columns
  ratio_columns <- character()
  for ( c in colnames(df) ) {
    #check the columns one by one
    if ( grepl("Abundance.Ratio", c) ) {
      new_colname <- strsplit(c, "\\.")[[1]][[3]]
      names(df)[names(df) == c] <- new_colname
      #all_proteins <- rename(all_proteins, new_colname = c)
      ratio_columns <- c(ratio_columns, new_colname)
    }
  }
  list(df, ratio_columns)
}

res <- rename_ratios(all_proteins)
all_proteins <- res[[1]]
quan_columns <- res[[2]]
quan_columns

dfWide <- all_proteins %>%
  filter(!grepl("cont_",Accession)) %>%
  subset (select=c("Accession", quan_columns) ) %>%
  na.omit()

rownames(dfWide) <- dfWide$Accession
dfWide$Accession <- NULL
dfWide <- log2(dfWide)
#Look at the distribution of quan values
summary(dfWide)
#Box Plot
boxplot(
  Log2_Abund~Sample, data = gather(dfWide, Sample, Log2_Abund),
  main = "Original Log2 Ratios"
  )
#Check how the distributions change after normalization on median
#For each column, subtract the median of the column from each of it's values
dfNorm <- mapply('-', dfWide, apply(dfWide,2,median))
#Transform into a dataframe
dfNorm <- as.data.frame(dfNorm, row.names = row.names(dfWide))
boxplot(
  Log2_Abund~Sample, data = gather(dfNorm, Sample, Log2_Abund),
  main = "Normalized Log2 Ratios"
  )
#If desired, use the normalized data for the downstream operations
dfWide <- dfNorm
dfWide$Accession <- row.names(dfWide)
dfLong <- gather(dfWide, Sample, Log2_Abund, -Accession)
dfWide$Accession <- NULL
head(dfWide)

#Add group names based on the samples
dfLong$Group <- apply(
  dfLong, 1, function(x) {
    sampleName <- x["Sample"]
    #OBS: format-dependent operations
    strsplit(sampleName, "_")[[1]][[1]]
  }
)
head(dfLong)

#PLOTS

#Histograms of abundance by sample
ggplot( dfLong, aes(x = Log2_Abund) ) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(title="Histograms on Proteins", x="Log2 TMT Ratio") +
  facet_wrap( ~Sample, scales = "fixed") +
  geom_histogram(bins=100) +
  theme_bw()

#Calculate and plot principal components (on samples)
# Calculate principal components (for PCA)
dfWide.t <-  t(dfWide)
dfWide.pca <- prcomp(dfWide.t, center = TRUE, scale. = FALSE)
summary(dfWide.pca)
dfWide.pca <- as.data.frame(dfWide.pca$x)
dfWide.pca$Group <- sapply(
  as.character( row.names(dfWide.pca) ),
  function(x) {
    strsplit(x, "_")[[1]][[1]]
  }
)
#Principal components 1 and 2
ggplot(
  dfWide.pca,
  aes(x = PC1, y = PC2, colour = Group )
  ) +
  geom_point(shape=19, size=4, alpha = 0.7)+
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  ggtitle("PCA On Proteins") +
  theme_classic()
#Principal components 3 and 4
ggplot(
  dfWide.pca,
  aes(x = PC3, y = PC4, colour = Group )
) +
  geom_point(shape=19, size=4, alpha = 0.7)+
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  ggtitle("PCA On Proteins") +
  theme_classic()

#Pairwise scatter plots
#Consider this for the tables of the moderate size
#This is probably not suitable for datasets with many dozens of samples
#For those, use the correlation heatmap
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}
upper.panel<-function(x, y){
  points(x,y, pch = 19, cex = 0.4, col = alpha("blue", 0.2))
}
pairs(dfWide, lower.panel = panel.cor, upper.panel = upper.panel,
      main = "Correlations on Proteins")

#Correlation heatmap
correlation_mat <- as.data.frame(round(cor(dfWide),2))
correlation_mat$S1 <- row.names(correlation_mat)
correlation_mat <- gather(correlation_mat, S2, Pearson_Corr, -S1)
ggplot(
    data = correlation_mat,
    aes(x = S1, y = S2, fill = Pearson_Corr)
  ) + 
  geom_tile() +
  #Add customized coor gradient from -1 to 1
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "grey", 
    midpoint = 0, limit = c(-1,1),
    
  ) +
  labs(title="Pearson Correlations on Proteins") +
  theme(
    axis.title.x=element_blank(), axis.ticks.x=element_blank(),
    axis.title.y=element_blank(), axis.ticks.y=element_blank()
  )


#Calculate one-way ANOVA P-values
dfANOVA <- dfWide
cols_anova <- c("his4_1", "his4_2", "his4_3", "met6_1", "met6_2", "met6_3", "ura2_1", "ura2_2", "ura2_3")
#Create the named list with our groups
groups <- list()
for (i in cols_anova) {
  local_name <- strsplit(i, "_")[[1]][[1]]
  if ( local_name %in% names(groups) ) {
    groups[[local_name]] <- c(groups[[local_name]], i)
  } else {
    groups[[local_name]] <- c(i)
  }
}
groups

dfANOVA$anovaPval <- apply(dfANOVA, 1, function(x) {
  df <- as.data.frame(x)
  #Select the samples for ANOVA
  #Important if you need to exclude some of the samples from the calculation
  df$Sample <- rownames(df)
  df <- df[ df$Sample %in% cols_anova, ]
  
  #Define groups in sync with the selected columns
  #OBS: format-dependent
  df$Group <- as.vector(
    sapply(
      cols_anova,
      function(x) { strsplit(x, "_")[[1]][[1]] }
      )
    )
  anovaResults <- aov(x ~ Group, data = df)
  #This Very exciting expression is how to extract the p-value from the aov summary
  return(summary(anovaResults)[[1]]["Pr(>F)"][[1]][[1]])
})
#Benjamini-Hochberg correction for multiple testing
dfANOVA$adjPval <- p.adjust(dfANOVA$anovaPval, method = "BH")
#Add group averages
for ( i in names(groups) ) {
  dfANOVA[i] <- apply(
    dfANOVA, 1, function(x) {
      #print(x)
      #print(typeof(x))
      mean( x[ groups[[i]] ] )
    }
  )
}
#Add max LogFC
dfANOVA$MaxLog2FC <- apply(
  dfANOVA, 1, function(x) {
    log_means <- x[ names(groups) ]
    max(log_means) - min(log_means)
  }
)

#Filter on adjusted p-value and minimal fold-change
print("Proteins with adjusted ANOVA P <= 0.05 and FC > log2(1.3) [30%]")
dfANOVA.Sign <- dfANOVA %>%
  filter(adjPval <= 0.05 & MaxLog2FC >= log2(1.3) ) %>%
  select(cols_anova)
dim(dfANOVA.Sign)
heatmap(
  as.matrix( dfANOVA.Sign ), scale = "none",
  col = colorRampPalette(c("red", "grey", "blue"))(256)
)
#Export the significant proteins according to ANOVA
#Add the info from the original protein table before exporting
dfANOVA %>%
  filter(adjPval <= 0.05 & MaxLog2FC >= log2(1.3) ) %>%
  merge(all_proteins, by.x="row.names", by.y="Accession",  suffixes=c("", "_"), sort=FALSE) %>%
  write.table("TKO_ANOVA_Significant_fromR.txt", sep="\t", na="", row.names = FALSE)

#Look at pairwise comparisons and create volcano plots
calc_ttest <- function(df, groupping, gr1, gr2, maxAdjP, minFC) {
  df <- df[ c( groupping[[gr1]], groupping[[gr2]]  ) ]
  #Log2 fold change group2 - group1
  df$Log2FC <- apply(
    df, 1, function(x) {
      mean( x[ groupping[[gr2]] ] ) - mean( x[ groupping[[gr1]] ] )
    }
  )
  #T-test with equal variance
  df$T_Pval <- apply(
    df, 1, function(x) {
      res <- t.test(
        x[ groupping[[gr2]] ], x[ groupping[[gr1]] ],
        alternative = "two.sided", var.equal = TRUE
        )
      mean( x[ groupping[[gr2]] ] ) - mean( x[ groupping[[gr1]] ] )
      res$p.value
    }
  )
  #Benjamini-Hochberg correction for multiple testing
  df$adjPval <- p.adjust(df$T_Pval, method = "BH")
  df$Log10adjPval <- -1*log10(df$adjPval)
  #Add the categorical column
  df$Diff_Abund <- apply(
    df, 1, function(x) {
      if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] >= minFC) {
        return( paste("Up in", gr2) )
      } else if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] <= -1*minFC) {
        return( paste("Up in", gr1) )
      } else {
        return('Non-significant')
      }
    }
  )
  df
}
maxAdjP <- 0.05
minLog2FC <- round(log2(1.3), 3)
gr1 <- "met6"
gr2 <- "his4"
dfTtest <- calc_ttest(dfWide, groups, gr1, gr2, maxAdjP, minLog2FC )
#How many proteins have adj. P <= 0.05?
dim( filter(dfTtest, adjPval <= maxAdjP) )
#How many proteins are differentially abundant according to our criteria?
dim(dfTtest %>%
      filter(adjPval <= maxAdjP) %>%
      filter( Log2FC >= minLog2FC | Log2FC <= -1*minLog2FC ))
head(dfTtest)
#Add columns with the gene names and other info
dfTtest <- merge(
  dfTtest, all_proteins,
  by.x="row.names", by.y="Accession",
  suffixes=c("", "_"), sort=FALSE
  )

#Save the differentially abundant proteins
dfTtest %>%
  filter(adjPval <= maxAdjP) %>%
  filter( Log2FC >= minLog2FC | Log2FC <= -1*minLog2FC ) %>%
  write.table(
    paste("TKO_Ttest_", gr1, "_vs_", gr2, "_Significant_fromR.txt", sep=""),
    sep="\t", na="", row.names = FALSE
    )
#Volcano plot
ggplot(
  dfTtest,
  aes(x = Log2FC, y = Log10adjPval, colour = Diff_Abund )
) +
  geom_point(shape=19, size=2, alpha = 0.6)+
  geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_vline(xintercept = -1*minLog2FC, colour = "gray65") +
  geom_vline(xintercept = minLog2FC, colour = "gray65") +
  ggtitle(
    paste(
      "T-test ", gr1, " vs ", gr2,
      " Adjusted P-value<=", maxAdjP, " Log2 FC>=", minLog2FC,
      sep=""
      )
    ) +
  theme_classic() +
  theme(
    legend.title = element_blank(), legend.text = element_text(size=12),
    plot.title = element_text(size=16)
    ) +
  labs(x = paste("Log2 FC", gr2, "-", gr1), y = "-Log10 Adj. P-value" ) +
  geom_text(
    data = subset(dfTtest, Log2FC >=0.9 | Log2FC <= -0.8),
    aes( Log2FC, Log10adjPval, label = Gene),
    alpha = 0.6, hjust = 0.5, vjust = -0.6
    )

#Use DEqMS instead of the T test or ANOVA
dfD <- dfWide[cols_anova]
#Define the design vector
cond = as.factor(
  c("his4", "his4", "his4", "met6", "met6", "met6", "ura2", "ura2", "ura2")
)
design = model.matrix(~0+cond)
colnames(design) = gsub("cond","",colnames(design))
#Make contrasts
x <- c(
  "his4-met6", "his4-ura2", "ura2-met6" 
  )
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(dfD, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)
#Extract PSM count information
psm_count_table <- dfD %>%
  merge(
    all_proteins[c("Accession", "Number.of.PSMs")],
    by.x="row.names", by.y="Accession",  suffixes=c("", "_"), sort=FALSE
    )
row.names(psm_count_table) <- psm_count_table$Row.names
psm_count_table <- psm_count_table[c("Number.of.PSMs")]
fit3$count = psm_count_table[rownames(fit3$coefficients),"Number.of.PSMs"]
fit4 = spectraCounteBayes(fit3)
#Visualize the DEqMS calculations
VarianceBoxplot(
  fit4, n=30, main="TKO Variance according to DEqMS", xlab="PSM count"
  )
#Extract the results
head(fit4$coefficients)
#Let's look at the contrast ura2-met6, column 3
current_contrast <- colnames(fit4$coefficients)[[1]]
res_DEqMS = outputResult(fit4, coef_col = 1)
head(res_DEqMS)
#adj.P.Val in the results comes from LIMMA, and sca.P.Value is from DEqMS
res_DEqMS$log.adj.P.Val  = -log10(res_DEqMS$adj.P.Val )
res_DEqMS$log.sca.adj.pval = -log10(res_DEqMS$sca.adj.pval)

maxAdjP <- 0.05
minLog2FC <- round(log2(1.3), 3)

#How many proteins have adj. P <= 0.05 according to LIMMA?
dim(
  filter(res_DEqMS, adj.P.Val <= maxAdjP)
)
#How many proteins have adj. P <= 0.05 and Log2 FC >= log2(1.3) [30%]
#according to LIMMA?
dim(res_DEqMS %>%
      filter(adj.P.Val <= maxAdjP) %>%
      filter( logFC >= minLog2FC | logFC <= -1*minLog2FC ))

#How many proteins have adj. P <= 0.05 according to DEqMS?
dim(
  filter(res_DEqMS, sca.adj.pval <= maxAdjP)
)
#How many proteins have adj. P <= 0.05 and Log2 FC >= log2(1.3) [30%]?
dim(res_DEqMS %>%
      filter(sca.adj.pval <= maxAdjP) %>%
      filter( logFC >= minLog2FC | logFC <= -1*minLog2FC ))
#Export the filtered DEqMS results
res_DEqMS %>%
  filter(sca.adj.pval <= maxAdjP) %>%
  filter( logFC >= minLog2FC | logFC <= -1*minLog2FC ) %>%
  write.table(
    paste("TKO_DEqMS_", current_contrast, "_Significant_fromR.txt", sep=""),
    sep="\t", na="", row.names = FALSE
  )
#current_contrast == paste(gr2, "-", gr1, sep="")
gr1 <- strsplit(current_contrast, "-")[[1]][[2]]
gr2 <- strsplit(current_contrast, "-")[[1]][[1]]
#Add categorical annotation columns for LIMMA and DEqMS
res_DEqMS$Diff_LIMMA <- apply(
  res_DEqMS, 1, function(x) {
    local_p <- as.numeric(x[["adj.P.Val"]])
    local_fc <- as.numeric(x[["logFC"]])
    if (local_p <= maxAdjP & local_fc >= minLog2FC) {
      return( paste("Up in", gr2) )
    } else if (local_p <= maxAdjP & local_fc <= -1*minLog2FC) {
      return( paste("Up in", gr1) )
    } else {
      return('Non-significant')
    }
  }
)
res_DEqMS$Diff_DEqMS <- apply(
  res_DEqMS, 1, function(x) {
    local_p <- as.numeric(x[["sca.adj.pval"]])
    local_fc <- as.numeric(x[["logFC"]])
    if (local_p <= maxAdjP & local_fc >= minLog2FC) {
      return( paste("Up in", gr2) )
    } else if (local_p <= maxAdjP & local_fc <= -1*minLog2FC) {
      return( paste("Up in", gr1) )
    } else {
      return('Non-significant')
    }
  }
)
head(res_DEqMS)
res_DEqMS <- merge(
  res_DEqMS, all_proteins,
  by.x="row.names", by.y="Accession",
  suffixes=c("", "_"), sort=FALSE
)
dim(res_DEqMS)
#Export the Whole DEqMS table
write.table(
  res_DEqMS,
  paste("TKO_DEqMS_", current_contrast, "_AllProteins.txt", sep=""),
  sep="\t", na="", row.names = FALSE
)

#Volcano plot for LIMMA
ggplot(
  res_DEqMS,
  aes(x = logFC, y = log.adj.P.Val, colour = Diff_LIMMA )
) +
  geom_point(shape=19, size=2, alpha = 0.6)+
  geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_vline(xintercept = -1*minLog2FC, colour = "gray65") +
  geom_vline(xintercept = minLog2FC, colour = "gray65") +
  ggtitle(
    paste(
      "LIMMA ", current_contrast,
      " Adjusted P-value<=", maxAdjP, " Log2 FC>=", minLog2FC,
      sep=""
    )
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(), legend.text = element_text(size=12),
    plot.title = element_text(size=16)
  ) +
  labs(x = paste("Log2 FC",current_contrast), y = "-Log10 Adj. P-value" ) +
  geom_text(
    data = subset(res_DEqMS, logFC >=1.0 | logFC <= -0.8),
    aes( logFC, log.adj.P.Val, label = Gene),
    alpha = 0.6, hjust = 0.5, vjust = -0.6
  )

#Volcano plot for DEqMS
ggplot(
  res_DEqMS,
  aes(x = logFC, y = log.sca.adj.pval, colour = Diff_DEqMS )
) +
  geom_point(shape=19, size=2, alpha = 0.6)+
  geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_vline(xintercept = -1*minLog2FC, colour = "gray65") +
  geom_vline(xintercept = minLog2FC, colour = "gray65") +
  ggtitle(
    paste(
      "DEqMS ", current_contrast,
      " Adjusted P-value<=", maxAdjP, " Log2 FC>=", minLog2FC,
      sep=""
    )
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(), legend.text = element_text(size=12),
    plot.title = element_text(size=16)
  ) +
  labs(x = paste("Log2 FC",current_contrast), y = "-Log10 Adj. P-value" ) +
  geom_text(
    data = subset(res_DEqMS, logFC >=0.9 | logFC <= -0.8),
    aes( logFC, log.sca.adj.pval, label = Gene),
    alpha = 0.6, hjust = 0.5, vjust = -0.6
  )

#Functional annotation with fgsea
#Download the GO Molecular Function, GO Biological Process from Uniprot
#Available under Creative Commons license
#Converted to .gmt files that are common for GSEA
annotations <- gmtPathways("go-bp_gene-symbols_uniprot.gmt")

#Prepare the rank file for the his4-met6 contrast
#Replace the UNIPROT accession with the Gene name
dfGSEA <- res_DEqMS[c("Gene", "logFC", "log.sca.adj.pval")]
#row.names(dfGSEA) <- res_DEqMS$Gene
#Calculate the rank as a product Log2FC * -Log10(adj.P.val.SCA)
dfGSEA$rank <- apply(
  dfGSEA, 1, function(x) {
    as.numeric(x[[2]]) * as.numeric(x[[3]])
  }
)
head(dfGSEA)
#Created the named vector of ranks for GSEA
ranks <- as.numeric(dfGSEA$rank)
names(ranks) <- dfGSEA$Gene
ranks <- sort(ranks)
ranks[1:10]
#Enrichment analysis
res_gsea <- fgseaMultilevel(annotations, ranks, minSize=15, maxSize=500)
#res_gsea <- fgsea(pathways, ranks, nperm=1000, minSize=15, maxSize=500)
head(res_gsea)
#head(res_gsea[order(pval), ])
res_gsea[ES > 0][head(order(pval), n=10),]
res_gsea[ES < 0][head(order(pval), n=10),]
plotEnrichment(
  annotations[["methionine biosynthetic process [GO:0009086]"]],
  ranks
) + labs(title="methionine biosynthetic process [GO:0009086]")
annotations[["methionine biosynthetic process [GO:0009086]"]]
plotEnrichment(
  annotations[["one-carbon metabolic process [GO:0006730]"]],
  ranks
) + labs(title="one-carbon metabolic process [GO:0006730]")
annotations[["one-carbon metabolic process [GO:0006730]"]]
