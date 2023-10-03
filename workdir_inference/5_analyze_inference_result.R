library(dplyr)
library(stringr)

getwd()
setwd("/Users/lovely_shufan/Dropbox (Edison_Lab@UGA)/AMF/AMF Imaging 2022/0_inference_using_MaskRCNN_2021/2_infer_result/GA_GWAS_2022/")
list.files()
data_blk2 = read.csv("block2_segmentation.csv", header = TRUE, stringsAsFactors = FALSE) 

table(data_blk2$filename, data_blk2$annotations)

# get accession name for block2 data
filenames <- strsplit(data_blk2$filename, split = "/")
plantid = c()
for (i in 1:length(filenames)){
  tmp = filenames[[i]][7]
  id = substr(tmp, 1, nchar(tmp)-4)
  plantid = c(plantid, id)
}
data_blk2$plantid <- plantid
unique(plantid)

# fix plant ids
data_blk2[data_blk2$plantid == "e2.1_PI656065_C0", ]$plantid <- "PI656065"

# remove duplicates of a plant
data_blk2_cleaned <- data_blk2
rmlist <- c("PI329541_C10_R09-1", "PI329541_C10_R09-2", "PI330168_C11_R09-1","PI330168_C11_R09-2",
            "PI521280_C06_R06-1", "PI569418_C14_R19-1", "PI569418_C14_R19-2")
data_blk2_cleaned[!data_blk2_cleaned$plantid %in% rmlist,] -> data_blk2_cleaned

# get accession name 
accession <- substr(data_blk2_cleaned$plantid, 1,8)
accession[accession == "PI92270_"] <- "PI92270"
accession[accession == "PI19770_"] <- "PI19770"

data_blk2_cleaned$accession <- accession

# calculate fungal structure frequency by accession
data_blk2_cleaned %>% group_by(accession, annotations) %>% summarise(freq = n()) -> table_frq

# calculate total root area per accession
data_blk2_cleaned[data_blk2_cleaned$annotations == "root",] %>% group_by(accession) %>% summarise(rootArea = sum(area)) -> rootarea

# calculate count density
merge(table_frq, rootarea, by = "accession", all.x = TRUE) -> count_density
count_density[!count_density$annotations == "root",] -> count_density
count_density[!count_density$annotations == "others",] -> count_density
count_density$countDensity <- count_density$freq / count_density$rootArea

# calculate percent colonization
data_blk2_cleaned[!data_blk2_cleaned$annotations == "root" & !data_blk2_cleaned$annotations == "others",] %>% 
  group_by(accession, annotations) %>% 
  summarise(sumarea = sum(area)) -> sumfungalarea

merge(sumfungalarea, rootarea, by = "accession", all.x = TRUE) -> percent_col
percent_col$percentCol <- percent_col$sumarea / percent_col$rootArea

# calculate overall percent colonization
data_blk2_cleaned[!data_blk2_cleaned$annotations == "root" & !data_blk2_cleaned$annotations == "others",] %>% 
  group_by(accession) %>% 
  summarise(sumarea = sum(area)) -> totalfungalarea

merge(totalfungalarea, rootarea, by = "accession", all.x = TRUE) -> total_col
total_col$percentCol <- total_col$sumarea / total_col$rootArea
# order accessions by percent colonization
 
total_col[order(total_col$percentCol, decreasing = TRUE),] -> total_col
pc_order <- total_col$accession
factor(percent_col$accession, levels = pc_order) -> percent_col$accession

# order accessions by count density
table_frq[!table_frq$annotations == "root" & ! table_frq$annotations == "others",] %>% 
  group_by(accession) %>% 
  summarise(sumcount = sum(freq)) -> totalfungalcount

merge(totalfungalcount, rootarea, by = "accession", all.x = TRUE) -> totalfungalcount
totalfungalcount$totalcd <- totalfungalcount$sumcount / totalfungalcount$rootArea

totalfungalcount[order(totalfungalcount$totalcd, decreasing = TRUE),] -> totalfungalcount
cd_order <- totalfungalcount$accession
factor(percent_col$accession, levels = cd_order) -> percent_col$accession

# visualization
library(ggplot2)
library(viridis)
ggplot(percent_col, aes(fill=annotations, y=percentCol, x=accession)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_viridis_d()

# visualization

ggplot(percent_col, aes(fill=annotations, y=percentCol, x=accession)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_viridis_d()

# visual 3
library(tidyr)
totalfungalcount$totalcd2 <- totalfungalcount$totalcd *10^4
merge(total_col, totalfungalcount, by = "accession") -> totals
totals[,c(1,4,8)] %>% gather(key = "phenotype", value = "value", -accession) -> total2

factor(total2$accession, levels = pc_order) -> total2$accession
ggplot(total2, aes(fill=phenotype, y=value, x=accession)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_viridis_d()





filenames <- strsplit(table_frq$filename, split = "/")
#filenames <- substr(table_frq$filename, nchar(table_frq$filename)-19,nchar(table_frq$filename)-4)
plantid = c()
for (i in 1:length(filenames)){
  tmp = filenames[[i]][7]
  id = substr(tmp, 1, nchar(tmp)-4)
  plantid = c(plantid, id)
}
table_frq$plantid <- plantid 
length(unique(plantid))

# fix plant ids
table_frq[table_frq$plantid == "e2.1_PI656065_C0", ]$plantid <- "PI656065"

# fix accession name
accession <- substr(table_frq$plantid, 1,8)
unique(accession)
accession[accession == "PI92270_"] <- "PI92270"
accession[accession == "PI19770_"] <- "PI19770"
table_frq$accession <- accession

# 

