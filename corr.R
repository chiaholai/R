setwd('C:/Users/Marco/Desktop/RefERP_analysis_2019/corr_20190604')
library(psych)

data <- read.csv('MeanAmp_YA34&OA34.csv', header=TRUE)
data.YA <- data[c(1:2)]
data.OA <- data[c(3:4)]

npsydata <- read.csv('SubInfo_Neuropsy.csv', header=TRUE)
npsy.YA <- subset(npsydata,Group == "YA")
npsy.OA <- subset(npsydata,Group == "OA")
npsy.YA.filt <- npsy.YA[c(9:11,13)]
npsy.OA.filt <- npsy.OA[c(9:11,13)]

corr.data.YA <- cbind.data.frame(data.YA,npsy.YA.filt)
corr.data.OA <- cbind.data.frame(data.OA,npsy.OA.filt)

corr.test.YA <- corr.test(corr.data.YA, corr.data.YA, method="pearson", adjust="fdr", alpha=.05)
corr.test.OA <- corr.test(corr.data.OA, corr.data.OA, method="pearson", adjust="fdr", alpha=.05)