getwd()

library(SNPRelate)
library(vcfR)
library(adegenet)
library(ape)
library(pegas)
library(hierfstat)
library(knitr)
library(StAMPP)
library(related)
library(ggplot2)
library(ggtree)

#first we'll look at the VCF file omit data with coverage outliers and high missing data
library(vcfR)
vcf <- read.vcfR("/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/populations.snps.vcf", verbose = FALSE )
vcf
head(vcf)
head(is.polymorphic(vcf, na.omit = TRUE))
head(is.biallelic(vcf))
#look at sequencing depth
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), xlab= "Sample", ylab="Number of reads per locus")
#filter out variants outside of 95% confidence interval
sums <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[1,])
dp[dp2 < 0] <- NA
dp2 <- sweep(dp, MARGIN=2, FUN = "-", sums[2,])
dp[dp2 > 0] <- NA
dp[dp < 4] <- NA
#plot data w/ extreme values omitted
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), xlab= "Sample", ylab="Number of reads per locus")
#current VCF
vcf
#set filtered values to NA
is.na(vcf@gt[,-1][ is.na(dp)]) <- TRUE
#new VCF;filtering by DP increases missing data
vcf 
#now look at genotype quality (GQ)
gq <- extract.gt(vcf, element='GQ', as.numeric=TRUE)
boxplot(gq, las=3, col=c("#C0C0C0", "#808080"), xlab= "Sample", ylab="Genotype quality")
#filter out variants with scores lower than 20
gq <- ifelse(gq<20, NA, gq)
boxplot(gq, las=3, col=c("#C0C0C0", "#808080"), xlab= "Sample", ylab="Genotype quality")
#current VCF
vcf
#set filtered values to NA
is.na(vcf@gt[,-1][ is.na(gq)]) <- TRUE
#new VCF;filtering by GQ further increases missing data
vcf 
#now look at missing data in samples and loci
#heatmap of all samples and first 1000 loci
heatmap.bp(dp[1:1000,], rlabels = FALSE)
#exclude samples with greater than 70% missing data
miss <- apply(dp, MARGIN = 2, function(x){sum(is.na(x))})
miss <- miss / nrow(dp)
vcf@gt <- vcf@gt[, c(TRUE, miss < 0.7)]
vcf
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp[1:1000,], rlabels = FALSE)
#exclude variants with greater than 20% missing data
miss2 <- apply(dp, MARGIN = 1, function(x){sum(is.na(x))})
miss2 <- miss2 / ncol(dp)
vcf <- vcf[miss2 < 0.2, ]
#in this case reduced number of variants by half
vcf
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp[1:1000,], rlabels = FALSE)
#write new VCF file
write.vcf(vcf, "/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/vcfR.filtered.vcf")

#use SNPRelate to perform LD-based SNP pruning 
library(SNPRelate)
vcf.fn <- "/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/vcfR.filtered.vcf"
#convert VCF file to GDS file for SNPRelate
snpgdsVCF2GDS(vcf.fn, "ab.gds", method="biallelic.only")
snpgdsSummary("ab.gds")
genofile <- snpgdsOpen("ab.gds")
#pop info - must make sure these are accurate w/ excluded data
pop_code <- rep(c("A", "B", "C", "D"),length.out=96)
table(pop_code)

#LD based SNP pruning
set.seed(1000)
#suggested to try diff thresholds; lower threshold is more conservative and will result in fewer makers
#run loop with different LD thresholds (loop not working so abandoned)
snpset0 <- snpgdsLDpruning(genofile, ld.threshold=0, autosome.only=FALSE)
snpset0.1 <- snpgdsLDpruning(genofile, ld.threshold=0.1, autosome.only=FALSE)
snpset0.2 <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only=FALSE)
snpset0.3 <- snpgdsLDpruning(genofile, ld.threshold=0.3, autosome.only=FALSE)
snpset0.4 <- snpgdsLDpruning(genofile, ld.threshold=0.4, autosome.only=FALSE)
snpset0.5 <- snpgdsLDpruning(genofile, ld.threshold=0.5, autosome.only=FALSE)
snpset0.6 <- snpgdsLDpruning(genofile, ld.threshold=0.6, autosome.only=FALSE)
snpset0.7 <- snpgdsLDpruning(genofile, ld.threshold=0.7, autosome.only=FALSE)
snpset0.8 <- snpgdsLDpruning(genofile, ld.threshold=0.8, autosome.only=FALSE)
snpset0.9 <- snpgdsLDpruning(genofile, ld.threshold=0.9, autosome.only=FALSE)
snpset1 <- snpgdsLDpruning(genofile, ld.threshold=1, autosome.only=FALSE)
numloci <- c(length(unlist(snpset0)), length(unlist(snpset0.1)), length(unlist(snpset0.2)), length(unlist(snpset0.3)),
             length(unlist(snpset0.4)), length(unlist(snpset0.5)), length(unlist(snpset0.6)), length(unlist(snpset0.7)),
             length(unlist(snpset0.8)), length(unlist(snpset0.9)), length(unlist(snpset1)))
thresh <- seq(0,1,0.1)
#plot retained loci vs LD threshold
plot(thresh, numloci, xlab = "LD threshold", ylab = "Retained SNP loci")
#chose ld threshold of 0.5
head(snpset0.5$chrScaffold_1)


# Get all selected snp id
snpset.id <- unlist(snpset0.5)
snpset.id2 <- as.data.frame(snpset.id)

#open VCF file to select included SNPs; exclude header
ab_vcf <- read.table("/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/populations.snps.vcf", skip = 15)
#select SNPs based on ID (index)
ab_vcf2 <- ab_vcf[rownames(ab_vcf) %in% snpset.id2$snpset.id,]
#export file
#write.table(ab_vcf2, "/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/populations.snps.select.vcf", 
#            quote = FALSE, row.names = FALSE, col.names = FALSE,sep="\t")

#### manually append header information to new VCF from old VCF in text editor

#convert VCF file to genlight
library(vcfR)
vcf <- read.vcfR("/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/populations.snps.select.vcf")
ab.genlight <- vcfR2genlight(vcf)
#use adgenet to add pop info (strata)
#create pop info df
#first, vector of sample ids
sample <- c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21","S22","S23","S24",
            "S25","S26","S27","S28","S29","S30","S31","S32","S33","S34","S35","S36","S37","S38","S39","S40","S41","S42","S43","S44","S45","S46","S47","S48",
            "S49","S50","S51","S52","S53","S54","S55","S56","S57","S58","S59","S60","S61","S62","S63","S64","S65","S66","S67","S68","S69","S70","S71","S72",
            "S73","S74","S75","S76","S77","S78","S79","S80","S81","S82","S83","S84","S85","S86","S87","S88","S89","S90","S91","S92","S93","S94","S95","S96")
#create dummy pop vector
pop <- rep(c("A", "B", "C", "D"),length.out=96)
info2 <- as.data.frame(cbind(sample, pop))
write.table(info2, file = "/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/pop.info", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
library(adegenet)
#can add pop info as strata (adegenet)
strata(ab.genlight) <- info2
strata(ab.genlight)
ab.genlight$strata
#or can add pop info as pop field 
pop(ab.genlight) <- paste(pop)
pop(ab.genlight)

#distribution of allele frequencies
myFreq <- glMean(ab.genlight)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba= TRUE, xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20) 
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=2)

#now convert VCF file for hierfstat
#load pegas hierfstat package
library(pegas)
library(hierfstat)
#need table with sample and pop info - note can include gender and other vars
vcf_select <- "/Users/jd/Library/Mobile\ Documents/com~apple~CloudDocs/PSRF/Abalone/SbfI-test-singleSNP/populations.snps.select.vcf"
table(info2$pop) #shows sample sizes per pop
x.1 <- VCFloci(vcf_select)
base <- c("A","T","G","C")
snps <- which(x.1$REF %in% base & x.1$ALT %in% base) 
x <- read.vcf(vcf_select,which.loci=snps)
vcf_hfstat <- genind2hierfstat(loci2genind(x),pop=info2$pop) 

hw <- hw.test(x)

#Basic statistics per population with hierfstat
basicstat <- basic.stats(vcf_hfstat, diploid = TRUE, digits = 2) 
names(basicstat)
basicstat_df <- as.data.frame(basicstat[c(3:5)])
basicstat_mean <- apply(basicstat_df, 2, FUN = mean, na.rm = TRUE)
popnames <- c("Pop A", "Pop B", "Pop C", "Pop D")
Ho <- basicstat_mean[1:4]
Hs <- basicstat_mean[5:8]
Fis <- basicstat_mean[9:12]
stat_table <- as.data.frame(cbind(popnames, Ho, Hs, Fis))
colnames(stat_table)[1] <- "Population"
stat_table$Ho <- as.numeric(as.character(stat_table$Ho))
stat_table$Hs <- as.numeric(as.character(stat_table$Hs))
stat_table$Fis <- as.numeric(as.character(stat_table$Fis))
#generate html table with 3 digits; later can add Ne estimates
library(knitr)
kable(stat_table, format = "html", digits = 3, row.names = FALSE)

#overall Fst and Fis
wc_stat <- wc(vcf_hfstat,diploid=TRUE)
#global Fst
wc_stat$FST
#pairwise W&C Fst
wc_pairwise <- pairwise.WCfst(vcf_hfstat,diploid=TRUE)
#html table
kable(wc_pairwise, format = "html", digits = 3, row.names = FALSE)


#pca (adegenet)
pca1 <- glPca(ab.genlight)
scatter(pca1, posi="topleft")

#pca scatter plot with pop groups (adegenet)
s.class(pca1$scores, as.factor(pop), col = c('#d7191c','#fdae61','#abdda4','#2b83ba'), 
        axesell = FALSE, grid = FALSE)
add.scatter.eig(pca1$eig,2,1,2)

#phylo tree
library(ape)
njtree <- nj(dist(as.matrix(ab.genlight)))
plot(njtree, typ="u", cex=0.7, show.tip.label = TRUE)
#make sure pop colors match
tiplabels(pch=20, col = c('#d7191c','#fdae61','#abdda4','#2b83ba'), cex=2)
segments(rep(0, 4), c(6,4,2,0), rep(2, 4), c(6,4,2,0), lwd = 5, col = c('#d7191c','#fdae61','#abdda4','#2b83ba'))
text(rep(2.5, 4), c(6,4,2,0), paste(c("A","B","C","D")), adj = 0)


#plot tree with ggtree
#library(ggplot2)
#library(ggtree)
#library(treeio)
#ggtree(njtree)
#figure out how to subset populations


############################################################################
#From Mac - relatedness using StaMPP package

head(ab.genlight)
# convert to allele frequencies in StAMPP
ab_freq <- stamppConvert(ab.genlight, type = "genlight")
(head(ab_freq))
#pop level Fst
stamppFst(ab_freq, nboots = 1000, percent = 95, nclusters = 4)


# calculate the genomic relationship matrix in StAMPP
pinto_Gmatrix <- stamppGmatrix(ab_freq)
head(pinto_Gmatrix)

#visualize gmatrix via heatmap
library(data.table)
library(pheatmap)
#need to modify the Gmatrix so the sample names are on both rows and columns before loading in to make the heatmap work
colnames(pinto_Gmatrix) <- rownames(pinto_Gmatrix)
pmap<-pheatmap(pinto_Gmatrix, cluster_rows=T, cluster_cols=T, clustering_distance_rows='euclidean', clustering_distance_cols='euclidean', clustering_method='ward', show_rownames=T, show_colnames=F)
x2 <- dist(pinto_Gmatrix)
x3 <- hclust(x2)
plot(x3)

# ##########################################################################
# #Relatedness using related package
# 
# #convert to format with alleles separated into columns
# rel_df <- loci2genind(x,pop=info2$pop, oneColPerAll = TRUE) 
# test <- as.loci(x)
# test2 <- genind2df(test, pop=info2$pop, oneColPerAll = TRUE)
# test3 <- genind2df(loci2genind(x), pop=info2$pop, oneColPerAll = TRUE)
# 
# 
# #create populations vector and merge with sample names
# pops <- rep(c("AA", "BB", "CC", "DD"),each = 24)
# pops2 <- paste(pops,row.names(rel_df),sep="_")
# 
# #bind pop/sample names to dataframe
# rel_df[,2:length(rel_df)] <- lapply(rel_df[,2:length(rel_df)], as.numeric)
# rel_df2 <- cbind(pops2, rel_df[,2:length(rel_df)])
# rel_df2$pops2 <- as.character(rel_df2$pops2)
# 
# #read in data for related package
# input <- readgenotypedata(rel_df2)
# #simulate 100 parentoffspring pairs, 100 full-sib pairs, 100 half-sib pairs, and 100 unrelated pairs 
# simdata <- familysim(sample(input$freqs, size=100, replace=F), 100)
# #run four different relatedness estimators on simulated data
# output <- coancestry(simdata, quellergt =1, wang =1, ritland =1, lynchli =1)
# simrel <- cleanuprvals(output$relatedness, 100)
# #select the range of rows and columns that correspond to the appropriate
# #relatedness value and estimator.
# wangpo <- simrel [1:100 , 6]
# wangfs <- simrel [(100 + 1) : (2 * 100) , 6]
# wanghs <- simrel [((2 * 100) + 1) : (3 * 100) , 6]
# wangur <- simrel [((3 * 100) + 1) : (4 * 100) , 6]
# lynchlipo <- simrel [1:100 , 7]
# lynchlifs <- simrel [(100 + 1) : (2 * 100) , 7]
# lynchlihs <- simrel [((2 * 100) + 1) : (3 * 100) , 7]
# lynchliur <- simrel [((3 * 100) + 1) : (4 * 100) , 7]
# ritlandpo <- simrel [1:100 , 9]
# ritlandfs <- simrel [(100 + 1) : (2 * 100) , 9]
# ritlandhs <- simrel [((2 * 100) + 1) : (3 * 100) , 9]
# ritlandur <- simrel [((3 * 100) + 1) : (4 * 100) , 9]
# quellergtpo <- simrel [1:100 , 10]
# quellergtfs <- simrel [(100 + 1) : (2 * 100) , 10]
# quellergths <- simrel [((2 * 100) + 1) : (3 * 100) , 10]
# quellergtur <- simrel [((3 * 100) + 1) : (4 * 100) , 10]
# #create a list of labels for the different estimators, with each repeated the appropriate number
# #of times (100 in this case).
# wang <- rep ("W", 100)
# lynchli <- rep ("LL", 100)
# ritland <- rep ("R", 100)
# quellergt <- rep ("QG", 100)
# estimator2 <- c(wang, lynchli, ritland, quellergt)
# Estimator <- rep(estimator2,4)
# #create a list of labels for the different relatedness types
# po <- rep ("Parent - Offspring", (4 * 100))
# fs <- rep ("Full - Sibs", (4 * 100))
# hs <- rep ("Half - Sibs ", (4 * 100))
# ur <- rep ("Unrelated", (4 * 100))
# relationship <- c(po, fs, hs, ur)
# #combine the different values for each estimator based on relatedness type, as lists
# relatednesspo <- c(wangpo, lynchlipo, ritlandpo, quellergtpo)
# relatednessfs <- c(wangfs, lynchlifs, ritlandfs, quellergtfs)
# relatednesshs <- c(wanghs, lynchlihs, ritlandhs, quellergths)
# relatednessur <- c(wangur, lynchliur, ritlandur, quellergtur)
# Relatedness_Value <- c(relatednesspo, relatednessfs, relatednesshs, relatednessur)
# #combine the data
# combineddata <- as.data.frame(cbind(Estimator, relationship, Relatedness_Value))
# combineddata$Relatedness_Value <- as.numeric(as.character(combineddata$Relatedness_Value))
# #plot
# ggplot(combineddata, aes(x = Estimator, y = Relatedness_Value) , ylim = c(-0.5, 1.0)) +
#   geom_boxplot () +
#   facet_wrap (~ relationship )
# #chose the Wang estimator based on low variance
# relvalues <- simrel [, 6]
# label1 <- rep ("PO", 100)
# label2 <- rep (" Full ", 100)
# label3 <- rep (" Half ", 100)
# label4 <- rep (" Unrelated ", 100)
# Relationship <- c(label1, label2, label3, label4 )
# newdata <- as.data.frame(cbind(Relationship, relvalues))
# newdata$relvalues <- as.numeric(as.character(newdata$relvalues))
# #density plot
# qplot(relvalues, ..density.. , data = newdata, geom ="density", colour = Relationship, 
#       xlab ="Relatedness Value", ylab ="Density")
# library("grid")
# grid.locator()
# 
# #run Wang estimator on full dataset
# related_wang <- coancestry(rel_df2, wang =1)
# library("reshape2")
# test <- as.data.frame(related_wang$relatedness)
# test2 <- acast(test, ind1.id ~ ind2.id ~ wang)
# heatmap(test2)
# 
# ########################################################################
# #relatedness using SNPrelate package
# set.seed(100)
# ibs <- snpgdsIBS(genofile, snp.id=snpset.id, autosome.only=FALSE, num.thread=2)
# pop.idx <- order(pop_code)
# image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))
# set.seed(100)
# ibs.hc <- snpgdsHCluster(ibs)
# rv <- snpgdsCutTree(ibs.hc)
# plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
# rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))
# plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
# 
