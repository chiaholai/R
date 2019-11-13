#chant <- c("FP1","FPZ","FP2","AF3","AF4",
#         "F7","F5","F3","F1","FZ","F2","F4","F6","F8",
#         "FT7","FC5","FC3","FC1","FCZ","FC2","FC4","FC6","FT8",
#         "C5","C3","C1","CZ","C2","C4","C6")

#chpos <- c("TP7","CP5","CP3","CP1","CPZ","CP2","CP4","CP6","TP8",
#         "P7","P5","P3","P1","PZ","P2","P4","P6","P8",
#         "PO7","PO5","PO3","POZ","PO4","PO6","PO8",
#         "O9","O1","OZ","O2","O10")

setwd('C:/Users/Marco/Desktop/RefERP_analysis_2019/ANOVA_afex/data')
YAdata <- read.delim('YA34_RefAmb_200to600_lop.txt',sep = '\t')
OAdata <- read.delim('OA34_RefAmb_200to600_lop.txt',sep = '\t')

YAdata$anteriority <- "NA"
YAdata$anteriority[grepl("P|O",YAdata$chlabel)] <- "pos"
YAdata$anteriority[grepl("F|C5$|C3$|C1$|CZ$|C2$|C4$|C6$",YAdata$chlabel)] <- "ant"

YAdata$condition <- "NA"
YAdata$condition[grepl("5",YAdata$bini)] <- "Unamb"
YAdata$condition[grepl("6",YAdata$bini)] <- "Amb"

electrode <- rep(1:30,136)
YAdata <- data.frame(YAdata,electrode)
YAdata$age <- "YA"

YAdata$subgroup <- "NA"
YAdata$subgroup[grepl("^subj001$|^subj002$|^subj006$|^subj007$|^subj008$|^subj012$|^subj013$|^subj015$|^subj018$|^subj019$|^subj022$|^subj023$|^subj024$|^subj025$|^subj026$|^subj028$|^subj029$|^YAsubj038$|^YAsubj041$|^YAsubj042$|^YAsubj043$",
                      YAdata$ERPset)] <- "N"
YAdata$subgroup[grepl("^subj003$|^subj005$|^subj009$|^subj016$|^subj020$|^subj021$|^subj027$|^YAsubj030$|^YAsubj031$|^YAsubj032$|^YAsubj033$|^YAsubj036$|^YAsubj037$",
                      YAdata$ERPset)] <- "P"

OAdata$anteriority <- "NA"
OAdata$anteriority[grepl("P|O",OAdata$chlabel)] <- "pos"
OAdata$anteriority[grepl("F|C5$|C3$|C1$|CZ$|C2$|C4$|C6$",OAdata$chlabel)] <- "ant"

OAdata$condition <- "NA"
OAdata$condition[grepl("5",OAdata$bini)] <- "Unamb"
OAdata$condition[grepl("6",OAdata$bini)] <- "Amb"

electrode <- rep(1:30,136)
OAdata <- data.frame(OAdata,electrode)
OAdata$age <- "OA"

OAdata$subgroup <- "NA"
OAdata$subgroup[grepl("^OAsubj002$|^OAsubj004$|^OAsubj007$|^OAsubj008$|^OAsubj009$|^OAsubj012$|^OAsubj013$|^OAsubj015$|^OAsubj016$|^OAsubj017$|^OAsubj021$|^OAsubj023$|^OAsubj024$|^OAsubj030$|^OAsubj032$|^OAsubj038$|^OAsubj040$",
                      OAdata$ERPset)] <- "N"
OAdata$subgroup[grepl("^OAsubj001$|^OAsubj003$|^OAsubj005$|^OAsubj006$|^OAsubj010$|^OAsubj014$|^OAsubj020$|^OAsubj025$|^OAsubj026$|^OAsubj027$|^OAsubj028$|^OAsubj029$|^OAsubj033$|^OAsubj034$|^OAsubj036$|^OAsubj037$|^OAsubj039$",
                      OAdata$ERPset)] <- "P"

data <- rbind(YAdata,OAdata)

data$age <- factor(data$age)
data$condition <- factor(data$condition)
data$anteriority <- factor(data$anteriority)
data$electrode <- factor(data$electrode)


require(afex)
require(emmeans)

stats.age <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "anteriority", "electrode"), between = "age", data = data)

##YA group
YAdata$condition <- factor(YAdata$condition)
YAdata$anteriority <- factor(YAdata$anteriority)
YAdata$electrode <- factor(YAdata$electrode)

stats.YA <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "anteriority", "electrode"), data = YAdata)
em.YA <- emmeans(stats.YA, spec="condition", by="anteriority")
simple.YA <- test(pairs(em.YA), joint = TRUE)

##OA group
OAdata$condition <- factor(OAdata$condition)
OAdata$anteriority <- factor(OAdata$anteriority)
OAdata$electrode <- factor(OAdata$electrode)

stats.OA <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "anteriority", "electrode"), data = OAdata)
em.OA <- emmeans(stats.OA, spec="condition", by="anteriority")
simple.OA <- test(pairs(em.OA), joint = TRUE)

YAdata_ant <- subset(YAdata,anteriority == "ant")
YAdata_pos <- subset(YAdata,anteriority == "pos")
OAdata_ant <- subset(OAdata,anteriority == "ant")
OAdata_pos <- subset(OAdata,anteriority == "pos")

stats.YA.ant <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "electrode"), data = YAdata_ant)
stats.YA.pos <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "electrode"), data = YAdata_pos)
stats.OA.ant <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "electrode"), data = OAdata_ant)
stats.OA.pos <- aov_ez(dv = "value", id = "ERPset", within = c("condition", "electrode"), data = OAdata_pos)

save(stats.age, file = "RefAmb200to600_btYA34&OA34.RData")
save(stats.YA, em.YA, simple.YA, stats.YA.ant, stats.YA.pos, file = "RefAmb200to600_wtYA34.RData")
save(stats.OA, em.OA, simple.OA, stats.OA.ant, stats.OA.pos, file = "RefAmb200to600_wtOA34.RData")