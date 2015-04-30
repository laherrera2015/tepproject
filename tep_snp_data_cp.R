##tep SNP analysis

#Set up
library(plyr)
library(ggplot2)
library(caTools)

setwd("~/tep_project/")
#setwd("~/git/tepproject") #for julin

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
colnames(M82_tmp) <- c("M82.gt","M82.read.depth","M82.ref.depth","M82.ref.qual","M82.alt.depth","M82.alt.qual","M82.lik")
head(M82_tmp)

#convert values to numeric
tep_tmp[,c(-1,-7)] <- apply(tep_tmp[,c(-1,-7)],2,as.numeric)

M82_tmp[,c(-1,-7)] <- apply(M82_tmp[,c(-1,-7)],2,as.numeric)

summary(tep_tmp)

summary(M82_tmp)

tep_data<- cbind(tep_data, tep_tmp, M82_tmp)

names(tep_data)

#next steps
# Filter to keep rows where M82 = "1/1"
head(tep_data$M82.gt, 15)
hmz_M82<-subset(tep_data, M82.gt == "1/1")

M82_snp_data<-hmz_M82[grep("TYPE=snp", hmz_M82$INFO),]

M82_snp_data <- M82_snp_data[!is.na(M82_snp_data$tep.read.depth),]
M82_snp_data <- M82_snp_data[!is.na(M82_snp_data$M82.read.depth),]

# Filtering on depth or quality for both M82 and tep

#mistake below in commented out line
#small_data<-sample(M82_snp_data[1:10000,])

small_data <- M82_snp_data[sample(nrow(M82_snp_data),10000),]

hist(small_data$tep.read.depth[small_data$tep.read.depth < 100],breaks=50)
hist(small_data$M82.read.depth[small_data$M82.read.depth < 100],breaks=50)
hist(small_data$tep.alt.qual[small_data$tep.alt.qual < 400],breaks=40)
hist(small_data$M82.alt.qual[small_data$M82.alt.qual < 6000],breaks=60) 

#Based on the histograms, filter by read depth
# 10<tep read depth<50 and 10<M82 read depth<70

hist(small_data$tep.read.depth[small_data$tep.read.depth < 50 & small_data$tep.read.depth > 3])
hist(small_data$M82.read.depth[small_data$M82.read.depth < 70 & small_data$M82.read.depth > 09])

#testing
small_tepdata<-small_data[small_data$tep.read.depth<50 & small_data$M82.read.depth <70,]
small_tepdata<-small_tepdata[small_tepdata$tep.read.depth> 9 & small_tepdata$M82.read.depth > 9,]

#Before filtering for read depth, the number of observations for M82_snp_data = 988,194

M82_snp_data<-M82_snp_data[ M82_snp_data$tep.read.depth < 50 & M82_snp_data$M82.read.depth < 70,]

#number of obs for M82_snp_data is 982,857

M82_snp_data<-M82_snp_data[ M82_snp_data$tep.read.depth > 9 & M82_snp_data$M82.read.depth > 9,]

#After filtering for read depth, number of obs is 824,860

#Also filter on quality

M82_snp_data<-M82_snp_data[ M82_snp_data$M82.alt.qual > 800,]

#now number of obs is 792,533

M82_snp_data<-M82_snp_data[ M82_snp_data$tep.alt.qual > 100 | M82_snp_data$tep.ref.qual > 100,]

#now 777,325

#Check to make sure filtering worked:
hist(M82_snp_data$tep.read.depth)
hist(M82_snp_data$M82.read.depth)
hist(M82_snp_data$tep.alt.qual,breaks=40)
hist(M82_snp_data$M82.alt.qual,breaks=40)

# Calculate for tep %M82
M82_snp_data$tep_percent <- (M82_snp_data$tep.alt.depth / M82_snp_data$tep.read.depth)*100

M82_snp_data$tep.runmean.15 <- unlist(tapply(M82_snp_data$tep_percent,M82_snp_data$CHROM,runmean,k=15,endrule="constant"))

#plot each chrosome and save to a pdf

pdf("tep1_runmean15_dep10.pdf",width=12,height=8)
for(chrom in unique(M82_snp_data$CHROM)) {
  pl <- ggplot(M82_snp_data[M82_snp_data$CHROM==chrom,],aes(x=POS))
  pl <- pl + geom_point(alpha=0.5,aes(y=tep_percent))
  pl <- pl + geom_line(color="blue",lwd=1,aes(y=tep.runmean.15))
  print(pl + ggtitle(chrom))
}
dev.off()

summary(M82_snp_data)

pl <- ggplot(M82_snp_data,aes(x=POS,y=tep_percent))
pl <- pl + geom_point(alpha=0.5)
pl <- pl + geom_smooth()
pl <- pl  + facet_wrap( ~ CHROM, ncol=2)
pl

#test1<-ggplot(M82_snp_data, aes(POS, tep_percent))+geom_point(alpha=0.5)
#test1+facet_wrap(~CHROM, ncol = 2)+geom_smooth()


pl <- ggplot(M82_snp_data[M82_snp_data$CHROM!="SL2.50ch00",],aes(x=POS))
pl <- pl + geom_point(alpha=0.5,aes(y=tep_percent))
pl <- pl + geom_line(color="blue",lwd=1,aes(y=tep.runmean.15))
pl <- pl  + facet_wrap( ~ CHROM, ncol=2)
pl

#M82_snp_data$M82_percent <- M82_snp_data$M82.read.depth / M82_snp_data$M82.ref.depth*100

#plot a single chromosome

pl <- ggplot(M82_snp_data[M82_snp_data$CHROM=="SL2.50ch01",],aes(x=POS))
pl <- pl + geom_point(alpha=0.5,aes(y=tep_percent))
pl <- pl + geom_line(color="blue",lwd=1,aes(y=tep.runmean.15))
pl




# plot tep %M82
# repeat as running average
