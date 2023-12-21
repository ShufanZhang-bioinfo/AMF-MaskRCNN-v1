## raw segmentation data and confidence score
library(stringr)
library(dplyr)
library(tidyr)

##read in data
#datatab=read.table("~/Downloads/segmentation (2).txt",header=TRUE,sep=",")
tmp.data = read.table("~/Dropbox (Edison_Lab@UGA)/AMF/AMF Imaging 2021/3_analysis/inference_on_14_RILS_accessions/inference using Georgia model/segmentation (2).txt",header=TRUE,sep=",")
tmp.data2 = tmp.data[,-6]

##split out accession, region slide and file name info
tmp.splittab <- tmp.data2[,"filename"] %>% str_split(string=.,pattern="\\/",simplify=TRUE) 
as.data.frame(tmp.splittab[,c(13:16)]) -> tmp.data3
tmp.data3$V1 %>% str_split( string=.,pattern = "-",simplify=TRUE) -> tmp.data4
tmp.data4[,1] -> tmp.data3$V1
cbind(tmp.data3, tmp.data2) -> seg
colnames(seg) <- c("accession", "region", "slide","imagename", "filename", "fileid", 
                   "height", "width", "class", "pixel.area", "confidence.score")

#clear tmp files
rm(list=ls(pattern="^tmp."))

#data cleaning: remove accessions used in training
seg[!seg$accession %in% c("E46", "EZY"),] -> seg

##create a key that can group the segments by image
seg[order(seg$filename),] -> seg
as.vector(table(seg$filename)) -> freq # table() sorts the element by alphabetic order

id <- c()
for (i in freq){
  id <- c(id,seq_len(i))
}

seg$segid <- id
seg$key <- paste(seg$accession, seg$region, seg$slide, seg$fileid, seg$segid, sep = "-")
seg$key2 <- paste(seg$accession, seg$region, seg$slide, seg$fileid, sep = "-")
seg$key3 <- paste(seg$accession, seg$region, seg$slide, sep = "-")
seg[,c(13:15,1:3,6,12, 9:11)] -> seg2

# data cleaning: remove segments with zero pixels
table(seg2$pixel.area == 0)
seg2[!seg2$pixel.area == 0,] -> seg2

# data cleaning: remove dumpster category "others"
seg2[!seg2$class == "others", ] -> seg3

# separate amf catergories from roots
seg.amf <- seg3[!seg3$class == "root",]
seg.root <- seg3[seg3$class == "root",]

# it is observed that some images don't have roots but AMF; to verify whether that's abnormal
length(intersect(unique(seg.amf$fileid), unique(seg.root$fileid))) # there are 396 images that don't have roots but contain amf
length(unique(seg.amf$fileid)) - length(intersect(unique(seg.amf$fileid), unique(seg.root$fileid)))
intersect(unique(seg.amf$fileid), unique(seg.root$fileid)) -> withroot
setdiff(unique(seg.amf$fileid), unique(withroot)) -> noroot
noroot[sample(1:394,10)] -> checkid
seg[seg$fileid %in% checkid,] -> checklist # 80% images marked as no root don't have roots 

# calculate overall percent colonization and count density for each sample/slide
# the amount of root area per slide will be used in denominator
## calculate root pixel area per slide
seg.root %>%  group_by(accession, region, slide) %>% summarise(root.pixel.area = sum(pixel.area)) -> seg.root.sum # the amount of root area per slide will be used to calculate %colonization
seg.root.sum$key <- paste(seg.root.sum$accession, seg.root.sum$region, seg.root.sum$slide, sep = "-")

seg.amf %>% group_by(accession, region, slide) %>% summarise(amf.pixel.area = sum(pixel.area),
                                                             amf.count = n()) -> seg.amf.sum
seg.amf.sum$key <- paste(seg.amf.sum$accession, seg.amf.sum$region, seg.amf.sum$slide, sep = "-")
merge(seg.amf.sum, seg.root.sum[,4:5], by = "key") -> col.overall

col.overall$percent_colonization <- col.overall$amf.pixel.area / col.overall$root.pixel.area
col.overall$count_density <- col.overall$amf.count / col.overall$root.pixel.area

# calculate accessions orders by overall percent colonization and overall count density
seg.amf %>% group_by(accession) %>% summarise(amf.pixel.area = sum(pixel.area),
                                              amf.count = n()) -> seg.amf.sum.accession
seg.root %>% group_by(accession) %>% summarise(root.pixel.area = sum(pixel.area)) -> seg.root.sum.accession

merge(seg.amf.sum.accession,seg.root.sum.accession, by = "accession") -> col.accession
col.accession$total.pc <- col.accession$amf.pixel.area / col.accession$root.pixel.area
col.accession$total.cd <- col.accession$amf.count / col.accession$root.pixel.area
col.accession[order(col.accession$total.pc, decreasing = TRUE),]$accession -> accession_order_pc
col.accession[order(col.accession$total.cd, decreasing = TRUE),]$accession -> accession_order_cd

# order accession in col.overall by accession_order_pc
factor(col.overall$accession, levels = accession_order_pc) -> col.overall$accession

# calculate count and avg segmentation size for each AMF class
seg.amf %>% group_by(accession, region, slide, class) %>% summarise(count = n(),
                                                                    total.size = sum(pixel.area)) -> seg.amf.sum.class 
seg.amf.sum.class$avg.size <- seg.amf.sum.class$total.size / seg.amf.sum.class$count

# convert data from long to wide
seg.amf.sum.class$key <- paste(seg.amf.sum.class$accession, seg.amf.sum.class$region, seg.amf.sum.class$slide, sep = "-")

col.class.count.spread <- seg.amf.sum.class %>% subset(select = c('key', 'class', 'count')) %>% 
  spread(key = class, value = count)
colnames(col.class.count.spread) <- c("key", "arb_count", "exH_count", "inH_count", "sp_count", "ves_count")

col.class.size.spread <- seg.amf.sum.class %>% subset(select = c('key', 'class', 'avg.size')) %>% 
  spread(key = class, value = avg.size)
colnames(col.class.size.spread) <- c("key", "arb_size", "exH_size", "inH_size", "sp_size", "ves_size")

col.overall %>% subset(select = c("key", "accession", "region","slide", "percent_colonization", "count_density"))-> tmp

merge(tmp, col.class.count.spread, by = "key", all.x = TRUE) -> tmp2
merge(tmp2, col.class.size.spread, by = "key") -> col.overall.class.predictors

# convert NAs in size and count to 0
col.overall.class.predictors[is.na(col.overall.class.predictors)] <- 0

# scale all class level predictors
cbind(tmp, scale(col.overall.class.predictors[,7:16])) -> col.overall.class.predictors.scaled
colnames(col.overall.class.predictors.scaled)[7:16] <- c("arb_count_scaled", "exH_count_scaled", "inH_count_scaled", "sp_count_scaled", "ves_count_scaled",
                                                         "arb_size_scaled", "exH_size_scaled", "inH_size_scaled", "sp_size_scaled", "ves_size_scaled")

# order region
factor(col.overall.class.predictors.scaled$region, levels = c("MID","TOP", "BOT")) -> col.overall.class.predictors.scaled$region

# reformat
gather(col.overall.class.predictors.scaled[,c(1:4,7:16)], key = "traits", value = "value",-key,-accession,-region,-slide) -> col.overall.class.predictors.scaled.long
#===========================================================
# calculate 2nd version col.overall.class.predictor.scaled for class count densities
#===========================================================
merge(col.overall.class.predictors.scaled[,c(1:4,6:11)], col.overall[,c(1,7)]) -> col.overall.class.predictors.scaled2
for (i in c(6:10)){
  col.overall.class.predictors.scaled2[,i] <- col.overall.class.predictors.scaled2[,i] / col.overall.class.predictors.scaled2[,11]
}
colnames(col.overall.class.predictors.scaled2)[6:10] <- c("arb_cd_scaled", "exH_cd_scaled", "inH_cd_scaled", "sp_cd_scaled",  "ves_cd_scaled")

#===========================================================
# for plotting model pc and cd
#===========================================================
col.overall.class.predictors.scaled[,c(1:5,7,12)] -> pc
gather(pc, key = "predictor.type", value = "predictor.value", -key,-accession,-region,-slide,-percent_colonization)-> pc.long

col.overall.class.predictors.scaled[,c(1:4,6,7,8,10,11)] -> cd
gather(cd, key = "predictor.type", value = "predictor.value", -key,-accession,-region,-slide,-count_density)-> cd.long


#===========================================================
# for plotting count density ratio between class
#===========================================================
## combine slides, calculate count for each AMF class 
seg.amf %>% group_by(accession, region, class) %>% summarise(count = n()) -> seg.amf.sum.class2
## combine slides, calculate root area per accession
seg.root %>% group_by(accession, region) %>% summarise(root.pixel.area = sum(pixel.area)) -> seg.root.sum2

seg.amf.sum.class2$key = paste(seg.amf.sum.class2$accession, seg.amf.sum.class2$region, sep = "-") 
seg.root.sum2$key = paste(seg.root.sum2$accession, seg.root.sum2$region, sep = "-") 
merge(seg.amf.sum.class2, seg.root.sum2[,-c(1,2)], by = "key", all.x = TRUE) -> col.class.count2
col.class.count2$count_density <- col.class.count2$count / col.class.count2$root.pixel.area

col.class.count2[!col.class.count2$class == "AMF internal hypha",] -> col.class.count2.sub

## combine slides and regions, calculate count for each AMF class 
seg.amf %>% group_by(accession, class) %>% summarise(count = n()) -> seg.amf.sum.class3
## combine slides and regions, calculate root area per accession
seg.root %>% group_by(accession) %>% summarise(root.pixel.area = sum(pixel.area)) -> seg.root.sum3
## calculate count density
merge(seg.amf.sum.class3, seg.root.sum3, by = "accession", all.x = TRUE) -> col.class.count3
col.class.count3$count_density <- col.class.count3$count / col.class.count3$root.pixel.area
## exclude internal hypha as it is indicated as insignificant
col.class.count3[!col.class.count3$class == "AMF internal hypha",] -> col.class.count3.sub

## rename amf classes
class_levels <- c("AMF arbuscule", "AMF external hypha", "AMF internal hypha", "AMF spore", "AMF vesicle")
class_rename <- list("AMF arbuscule" = "arb", 
                     "AMF external hypha" = "exH",
                     "AMF internal hypha" = "inH", 
                     "AMF spore" = "sp", 
                     "AMF vesicle" = "ves")
factor(col.class.count3.sub$class, levels = class_levels) -> col.class.count3.sub$class
factor(col.class.count3$class, levels = class_levels) -> col.class.count3$class

factor(col.class.count2.sub$class, levels = class_levels) -> col.class.count2.sub$class
factor(col.class.count2$class, levels = class_levels) -> col.class.count2$class

col.class.count3.sub$class2 <- col.class.count3.sub$class
col.class.count3$class2 <- col.class.count3$class
col.class.count2.sub$class2 <- col.class.count2.sub$class
col.class.count2$class2 <- col.class.count2$class

levels(col.class.count3.sub$class2) <- c("arb", "exH","inH", "sp", "ves")
levels(col.class.count3$class2) <- c("arb", "exH","inH", "sp", "ves")
levels(col.class.count2.sub$class2) <- c("arb", "exH","inH", "sp", "ves")
levels(col.class.count2$class2) <- c("arb", "exH","inH", "sp", "ves")

#===========================================================
# for plotting class size
#===========================================================
factor(seg.amf.sum.class$class, levels = class_levels) -> seg.amf.sum.class$class
levels(seg.amf.sum.class$class) <- c("arb", "exH","inH", "sp", "ves")

#save.image("workspace_modeltable.rds")

# drafts to be organized
col.overall[,c(1:4,8,9)] -> tmp
tmp$percent_colonization^(1/3) -> tmp$trans_pc
tmp$count_density^(1/2) -> tmp$trans_cd
