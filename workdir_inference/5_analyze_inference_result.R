library(dplyr)
library(stringr)

getwd()
setwd("/Users/lovely_shufan/Dropbox (Edison_Lab@UGA)/AMF/AMF Imaging 2022/0_inference_using_MaskRCNN_2021/2_infer_result/GA_GWAS_2022/")
list.files()
data_blk2 = read.csv("block2_segmentation_image_complete.txt", header = TRUE, stringsAsFactors = FALSE) 
data_blk8 = read.csv("block8_segmentation_image_complete.txt", header = TRUE, stringsAsFactors = FALSE) 
data_blk10 = read.csv("block10_segmentation_image_complete.txt", header = TRUE, stringsAsFactors = FALSE) 
rbind(data_blk2, data_blk8, data_blk10) -> raw
##=============================================================================================================
##                                Data cleaning
##=============================================================================================================
# check for NAs
colSums(is.na(raw))
raw[!is.na(raw$area),] -> raw # removed segmentation with 0 area

# get accession names
accession <- str_extract(raw$filename, "(?<=/)[^/_]+(?=_[^/]+)")
raw$accession <- accession

# get field position
field <- str_extract(raw$filename, "C\\d+_R\\d+")
raw$fieldposition <- field 
table(nchar(field))

# clean up NAs in accession names
table(is.na(accession))
unique(raw[is.na(accession),]$filename)
accession[is.na(accession)] <- 'PI521280'

# examine abnormal accession names
abn_accession <- c("RecognizedCode-109-1", "GWAS",                
                   "RecognizedCode-77-2", "e2.1")
unique(raw[accession %in% abn_accession,]$filename)

# fix abnormal accession name 'e2.1' by replacing it with 'PI656065'
accession[accession == "e2.1"] <- "PI656065"

# fix abnormal accession name "GWAS" for 21 accessions using regex pattern
raw[!raw$accession == "GWAS", ] -> part1
raw[raw$accession == "GWAS", ] -> part2
part2$accession <- str_extract(part2$filename, "(?<=/)[^/]+(?=/[^/]+$)") 
rbind(part1, part2) -> raw

# temporarily remove "RecognizedCode-109-1" "RecognizedCode-77-2"
abn_accession <- c("RecognizedCode-109-1", "RecognizedCode-77-2")
raw[!raw$accession %in% abn_accession, ] -> seg
colSums(is.na(seg))

# get field position
field <- str_extract(seg$filename, "C\\d+_R\\d+")
table(nchar(field))
table(is.na(field))

# field position is na; Keep those as NAs; check field position latter
unique(seg[is.na(field),]$filename)
seg[is.na(seg$fieldposition),] -> tmp

# field positions with 6 characters are turned to NAs
field[nchar(field) == 6]<- NA
seg$fieldposition <- field

# remove replicates of the same plant
plant <- str_extract(seg$filename, "(?<=/)[^/]+(?=\\.czi)")
table(nchar(plant))

# rows with 62 character plant id has zero reps; keep
unique(plant[nchar(plant) == 62])

# rows with 20 character plant id are "PI641909_C17_R05-4-5" "PI641909_C17_R05-4-6", are duplicates, removed
unique(plant[nchar(plant) == 20])
seg[!nchar(plant) == 20, ] -> seg_nodup

plant <- str_extract(seg_nodup$filename, "(?<=/)[^/]+(?=\\.czi)")
table(nchar(plant))

# rows with 18 character plant id are duplicates; removed
unique(plant[nchar(plant) == 18])
seg_nodup[!nchar(plant) == 18, ] -> seg_nodup

plant <- str_extract(seg_nodup$filename, "(?<=/)[^/]+(?=\\.czi)")
table(nchar(plant))

# rows with 17 character plant id are duplicates; removed
unique(plant[nchar(plant) == 17])
seg_nodup[!nchar(plant) == 17, ] -> seg_nodup

plant <- str_extract(seg_nodup$filename, "(?<=/)[^/]+(?=\\.czi)")
table(nchar(plant))

# rows with 16 character plant id are not duplicates; KEEP
seg_nodup[seg_nodup$accession == "PI521280",] -> tmp # seg_nodup[plant == "_PI521280_C17_R1",]
unique(tmp$filename)
seg_nodup[seg_nodup$accession == "PI656065",] -> tmp # seg_nodup[plant == "e2.1_PI656065_C0",]
unique(tmp$filename)

# rows with 15 character plant id are not duplicates; KEEP
unique(plant[nchar(plant) == 15])

##=============================================================================================================
##              Initial exploration
##=============================================================================================================

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

