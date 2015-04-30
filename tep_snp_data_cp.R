##tep SNP analysis

#Set up
library(plyr)
library(ggplot2)
library(caTools)

setwd("~/tep_project/")
tep_data<-read.table("tep1_merge_cp.vcf", header = FALSE, stringsAsFactors=FALSE)
summary(tep_data)
head(tep_data)
colnames(tep_data) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","tep_calls","M82_calls")

#get rid of alleles with more than one alternate allele
tep_data <- tep_data[!grepl(",",tep_data$ALT),]

#fix empty data entry
tep_data$tep_calls[tep_data$tep_calls=="."] <- "NA:NA:NA:NA:NA:NA:NA,NA,NA"
tep_data$M82_calls[tep_data$M82_calls=="."] <- "NA:NA:NA:NA:NA:NA:NA,NA,NA"


tep_tmp <- ldply(strsplit(tep_data$tep_calls, split = ":"))
colnames(tep_tmp) <- c("tep.gt","tep.read.depth","tep.ref.depth","tep.ref.qual","tep.alt.depth","tep.alt.qual","tep.lik")
head(tep_tmp)

M82_tmp <- ldply(strsplit(tep_data$M82_calls, split = ":"))
colnames(M82_tmp) <- c("M82.gt","M82.read.depth","M82.ref.depth","M82.ref.qual","M82.alt.depth","M82.qual","M82.lik")
head(M82_tmp)


tep_tmp[,c(-1,-7)] <- apply(tep_tmp[,c(-1,-7)],2,as.numeric)

M82_tmp[,c(-1,-7)] <- apply(M82_tmp[,c(-1,-7)],2,as.numeric)

tep_data<- cbind(tep_data, tep_tmp, M82_tmp)
names(tep_data)

#next steps
# Filter to keep rows where M82 = "1/1"
head(tep_data$M82.gt, 15)
hmz_M82<-subset(tep_data, M82.gt == "1/1")


M82_snp_data<-hmz_M82[grep("TYPE=snp", hmz_M82$INFO),]


# Filtering on depth or quality for both M82 and tep
small_data<-sample(M82_snp_data[1:1000,])
hist(small_data$tep.read.depth[small_data$tep.read.depth < 100])
hist(small_data$M82.read.depth[small_data$M82.read.depth < 100])
hist(small_data$tep.alt.qual[small_data$tep.alt.qual < 400])
hist(small_data$M82.qual[small_data$M82.qual < 6000]) 

#Based on the histograms, filter by read depth
# 2<tep read depth<50 and 8<M82 read depth<70

hist(small_data$tep.read.depth[small_data$tep.read.depth < 50 & small_data$tep.read.depth > 2])
hist(small_data$M82.read.depth[small_data$M82.read.depth < 70 & small_data$M82.read.depth > 8])

#testing
small_tepdata<-small_data[small_data$tep.read.depth<50 & small_data$M82.read.depth <70,]
small_tepdata<-small_tepdata[small_tepdata$tep.read.depth> 2 & small_tepdata$M82.read.depth > 8,]

#Before filtering for read depth, the number of observations for M82_snp_data = 1,041,919

M82_snp_data<-M82_snp_data[ M82_snp_data$tep.read.depth < 50 & M82_snp_data$M82.read.depth < 70,]

#number of obs for M82_snp_data is 1,036,563

M82_snp_data<-M82_snp_data[ M82_snp_data$tep.read.depth > 2 & M82_snp_data$M82.read.depth > 8,]

#After filtering for read depth, number of obs is 878,566

#Check to make sure filtering worked:
hist(M82_snp_data$tep.read.depth)
hist(M82_snp_data$M82.read.depth)

# Calculate for tep %M82
M82_snp_data$tep_percent <- (M82_snp_data$tep.read.depth / M82_snp_data$tep.ref.depth)*100
#test1<-ggplot(M82_snp_data, aes(POS, tep_percent))+geom_point(alpha=0.5)
#test1+facet_wrap(~CHROM, ncol = 2)+geom_smooth()

M82_snp_data$M82_percent <- M82_snp_data$M82.read.depth / M82_snp_data$M82.ref.depth*100





# plot tep %M82
# repeat as running average
