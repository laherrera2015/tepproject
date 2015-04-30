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
tep_snp_data <-subset(tep_data, tep.gt != "0/0")

tep_snp_data<-tep_snp_data[grep("TYPE=snp", tep_snp_data$INFO),]

tep_snp_data <- tep_snp_data[!is.na(tep_snp_data$tep.read.depth),]

# Filtering on depth or quality for both M82 and tep

#mistake below in commented out line
#small_data<-sample(tep_snp_data[1:10000,])

small_data <- tep_snp_data[sample(nrow(tep_snp_data),10000),]

hist(small_data$tep.read.depth[small_data$tep.read.depth < 100],breaks=50)
hist(small_data$tep.alt.qual[small_data$tep.alt.qual < 400],breaks=40)

tep_snp_data<-tep_snp_data[ tep_snp_data$tep.read.depth < 60,]

tep_snp_data<-tep_snp_data[ tep_snp_data$tep.read.depth > 9,]

#Also filter on quality

#tep_snp_data<-tep_snp_data[ tep_snp_data$tep.alt.qual > 100 | tep_snp_data$tep.ref.qual > 100,]

#Check to make sure filtering worked:
hist(tep_snp_data$tep.read.depth)
hist(tep_snp_data$tep.alt.qual,breaks=40)

# Calculate for tep %M82
tep_snp_data$tep_percent <- (tep_snp_data$tep.alt.depth / tep_snp_data$tep.read.depth)*100

tep_snp_data$tep.runmean.15 <- unlist(tapply(tep_snp_data$tep_percent,tep_snp_data$CHROM,runmean,k=15,endrule="constant"))
tep_snp_data$tep.runmean.20 <- unlist(tapply(tep_snp_data$tep_percent,tep_snp_data$CHROM,runmean,k=20,endrule="constant"))


#plot each chrosome and save to a pdf

pdf("tep_only_tep1_runmean20_dep10.pdf",width=12,height=8)
for(chrom in unique(tep_snp_data$CHROM)) {
  pl <- ggplot(tep_snp_data[tep_snp_data$CHROM==chrom,],aes(x=POS))
  pl <- pl + geom_point(alpha=0.5,aes(y=tep_percent))
  pl <- pl + geom_line(color="blue",lwd=1,aes(y=tep.runmean.20))
  print(pl + ggtitle(chrom))
}
dev.off()

summary(tep_snp_data)

pl <- ggplot(tep_snp_data,aes(x=POS,y=tep_percent))
pl <- pl + geom_point(alpha=0.5)
pl <- pl + geom_smooth()
pl <- pl  + facet_wrap( ~ CHROM, ncol=2)
pl

#test1<-ggplot(tep_snp_data, aes(POS, tep_percent))+geom_point(alpha=0.5)
#test1+facet_wrap(~CHROM, ncol = 2)+geom_smooth()


pl <- ggplot(tep_snp_data[tep_snp_data$CHROM!="SL2.50ch00",],aes(x=POS))
pl <- pl + geom_point(alpha=0.5,aes(y=tep_percent))
pl <- pl + geom_line(color="blue",lwd=1,aes(y=tep.runmean.15))
pl <- pl  + facet_wrap( ~ CHROM, ncol=2)
pl

#tep_snp_data$M82_percent <- tep_snp_data$M82.read.depth / tep_snp_data$M82.ref.depth*100

#plot a single chromosome

pl <- ggplot(tep_snp_data[tep_snp_data$CHROM=="SL2.50ch01",],aes(x=POS))
pl <- pl + geom_point(alpha=0.5,aes(y=tep_percent))
pl <- pl + geom_line(color="blue",lwd=1,aes(y=tep.runmean.15))
pl




# plot tep %M82
# repeat as running average

