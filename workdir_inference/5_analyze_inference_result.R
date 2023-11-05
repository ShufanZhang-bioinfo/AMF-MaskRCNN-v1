library(dplyr)
library(stringr)
library(ggplot2)

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
raw$accession <- accession

# examine abnormal accession names
abn_accession <- c("RecognizedCode-109-1", "GWAS",                
                   "RecognizedCode-77-2", "e2.1")
unique(raw[accession %in% abn_accession,]$filename)

# fix abnormal accession name 'e2.1' by replacing it with 'PI656065'
accession[accession == "e2.1"] <- "PI656065"
raw$accession <- accession

# fix abnormal accession name "GWAS" for 21 accessions using regex pattern
raw[!raw$accession == "GWAS", ] -> part1
raw[raw$accession == "GWAS", ] -> part2
unique(part2$filename)
part2$accession <- str_extract(part2$filename, "(?<=/)[^/]+(?=/[^/]+$)") 
part2[part2$accession == "Block10",]$accession <- "PI569433"

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

# remove others
seg_nodup[!seg_nodup$annotations == "others",] -> seg_clean

# remove segmentation with 0 size
length(seg_clean$area[seg_clean$area == 0]) 
seg_clean <- seg_clean[!seg_clean$area == 0,]

# order fungal structures
factor(seg_clean$annotations, levels = c( "AMF vesicle", "AMF spore", "AMF internal hypha",
                                          "AMF external hypha","AMF arbuscule",
                                   "root")) -> seg_clean$annotations
# convert area unit from pixel to mm^2
seg_clean$area_mm2 <- seg_clean$area * 0.1202399 *10^-6

# some segmentations are very small; after examine tiles, decided to remove segmentations that are smaller than 250 pixels
seg_clean[seg_clean$area<1000,] -> small

ggplot(small, aes(x = area,fill=annotations)) + 
  geom_histogram(binwidth = 10) +
  labs(title = "Histogram of Area by Annotations")+
  facet_grid(rows = vars(annotations))+
  scale_fill_viridis_d()+
  theme_minimal()

seg_clean[!seg_clean$area < 250,] -> seg_clean_v2

# filter out tiles that have too much root and contains root only
seg_clean_v2 %>% group_by(filename, scene, tile, annotations) %>% summarise(freq = n(),rootArea = sum(area_mm2)) -> table_frq
table_frq_root <- table_frq[table_frq$annotations == "root",]

table_frq_root[table_frq_root$rootArea > 0.5910032,] -> toomuchroot
toomuchroot$key <- paste(toomuchroot$filename,toomuchroot$scene,toomuchroot$tile,
                         sep = "~")
seg_clean_v2[!seg_clean_v2$annotations == "root",] -> seg_amf
seg_amf$key <- paste(seg_amf$filename,seg_amf$scene,seg_amf$tile,sep = "~")
seg_clean_v2$key <- paste(seg_clean_v2$filename,seg_clean_v2$scene,seg_clean_v2$tile,sep = "~")

seg_clean_v2[seg_clean_v2$key %in% unique(toomuchroot$key),] -> tmp
merge(tmp, toomuchroot[,-c(1:4)], by = "key", all.x = TRUE) -> tmp1

tmp1[tmp1$freq %in% c(7,8,9,10,15),] -> tmp2 
seg_clean_v2[!seg_clean_v2$key %in% unique(tmp2$key),] -> seg_clean_v2

# remove duplicated segmentations by groups of file,scene,tile,annotations
seg_clean_v3 <- seg_clean_v2 %>%
  group_by(key,annotations) %>%
  filter(row_number() == 1) %>%
  ungroup()

# 10% tiles contain fungal structures but no roots; didn't remove any tiles
seg_clean_v3 %>% group_by(filename,scene,tile,annotations) %>% summarise(freq = n(), sumArea = sum(area_mm2)) -> freq_table 
freq_table$key <- paste(freq_table$filename,freq_table$scene,freq_table$tile,sep = "~")

keys_without_root <- freq_table %>%
  group_by(filename,scene,tile) %>%
  filter(all(annotations != "root")) %>%
  distinct(key) %>%
  ungroup()

length(keys_without_root$key) / length(unique(freq_table$key))
##=============================================================================================================
##              Initial exploration
##=============================================================================================================
ggplot(seg_clean_v3, aes(x = area,fill=annotations)) + 
  geom_histogram(binwidth = 10000) +
  labs(title = "Histogram of Area by Annotations")+
  facet_grid(rows = vars(annotations))+
  scale_fill_viridis_d()+
  theme_minimal()

ggplot(seg_clean_v3, aes(x = area_mm2,fill=annotations)) + 
  geom_histogram(binwidth = 0.001) +
  labs(title = "Histogram of Area by Annotations")+
  facet_grid(rows = vars(annotations))+
  scale_fill_viridis_d()+
  theme_minimal()

ggplot(seg_clean_v3, aes(x = confidenceScore,fill=annotations)) + 
  geom_histogram(binwidth = 0.005) +
  facet_grid(rows = vars(annotations), scales = "free_y")+
  labs(title = "Histogram of ConfidenceScore by Annotations")+
  scale_fill_viridis_d()+
  theme_minimal()

ggplot(seg_clean_v3, aes(x = area_mm2, y = confidenceScore, color=annotations)) + 
  geom_point() +
  facet_grid(rows = vars(annotations), scales = "free")+
  scale_color_viridis_d()+
  labs(title = "Scatterplot of area vs. confidenceScore")+
  theme_minimal()

##=============================================================================================================
##              Data wrangling
##=============================================================================================================
# calculate fungal structure frequency by accession and block
seg_clean_v3 %>% group_by(accession, block, annotations) %>% summarise(freq = n()) -> table_frq
table_frq$key <- paste(table_frq$accession, table_frq$block, sep = "_")

# calculate total root area per accession and block
seg_clean_v3[seg_clean_v3$annotations == "root",] %>% group_by(accession,block) %>% summarise(rootArea = sum(area)) -> rootarea
rootarea$key <- paste(rootarea$accession, rootarea$block, sep = "_")
rootarea$rootArea_mm2 <- rootarea$rootArea * 0.1202399 *10^-6

# calculate count density by block and annotations
merge(table_frq, rootarea[,-c(1,2)], by = "key", all.x = TRUE) -> count_density_blk_anno
count_density_blk_anno[!count_density_blk_anno$annotations == "root",] -> count_density_blk_anno
count_density_blk_anno$countDensity <- count_density_blk_anno$freq / count_density_blk_anno$rootArea_mm2

#===================================================================
# calculate fungal frequency by accession and block
seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% group_by(accession, block) %>% summarise(fungalFreq = n()) -> table_frq
table_frq$key <- paste(table_frq$accession, table_frq$block, sep = "_")

# calculate count density by block and annotations
merge(table_frq, rootarea[,-c(1,2)], by = "key", all.x = TRUE) -> count_density_blk
count_density_blk$countDensity <- count_density_blk$fungalFreq / count_density_blk$rootArea_mm2

#===================================================================

# calculate fungal structure percent colonization by accession and block
seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
  group_by(block, accession, annotations) %>% 
  summarise(FungalArea_mm2 = sum(area_mm2)) -> fungalarea
fungalarea$key <- paste(fungalarea$accession, fungalarea$block, sep = "_")

merge(fungalarea, rootarea[,-c(1:3)], by = "key", all.x = TRUE) -> percent_col_blk_anno
percent_col_blk_anno$percentCol <- percent_col_blk_anno$FungalArea_mm2 / percent_col_blk_anno$rootArea_mm2

# calculate fungal percent colonization by accession and block
seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
  group_by(block, accession) %>% 
  summarise(FungalArea_mm2 = sum(area_mm2)) -> fungalarea
fungalarea$key <- paste(fungalarea$accession, fungalarea$block, sep = "_")

merge(fungalarea, rootarea[,-c(1:3)], by = "key", all.x = TRUE) -> percent_col_blk
percent_col_blk$percentCol <- percent_col_blk$FungalArea_mm2 / percent_col_blk$rootArea_mm2

#===================================================================
# calculate fungal structure percent colonization by accession 
seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
  group_by(accession, annotations) %>% 
  summarise(FungalArea_mm2 = sum(area_mm2)) -> fungalarea

# calculate root area per accession
seg_clean_v3[seg_clean_v3$annotations == "root",] %>% group_by(accession) %>% summarise(rootArea = sum(area)) -> rootarea
rootarea$rootArea_mm2 <- rootarea$rootArea * 0.1202399 *10^-6

merge(fungalarea, rootarea, by = "accession", all.x = TRUE) -> percent_col_anno
percent_col_anno$percentCol <- percent_col_anno$FungalArea_mm2 / percent_col_anno$rootArea_mm2

# calculate fungal count density by accession 
seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
  group_by(accession, annotations) %>% 
  summarise(freq = n())-> fungalcount

merge(fungalcount, rootarea, by = "accession", all.x = TRUE) -> count_density_anno
count_density_anno$countDensity <- count_density_anno$freq / count_density_anno$rootArea_mm2

#===================================================================

seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
  group_by(accession) %>% 
  summarise(sumFungalArea = sum(area_mm2)) -> totalfungalarea

seg_clean_v3[seg_clean_v3$annotations == "root",] %>% 
  group_by(accession) %>% 
  summarise(sumRootArea = sum(area_mm2)) -> totalrootarea

merge(totalfungalarea, totalrootarea, by = "accession", all.x = TRUE) -> total_pc
total_pc$percentCol <- total_pc$sumFungalArea / total_pc$sumRootArea

# order accessions by percent colonization
total_pc[order(total_pc$percentCol, decreasing = TRUE),] -> total_pc
pc_order <- total_pc$accession

# order accessions by count density
seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% group_by(accession) %>% summarise(sumFungalCount = n()) -> totalfungalcount

merge(totalfungalcount, totalrootarea, by = "accession", all.x = TRUE) -> total_cd
total_cd$countDensity <- total_cd$sumFungalCount / total_cd$sumRootArea

total_cd[order(total_cd$countDensity, decreasing = TRUE),] -> total_cd
cd_order <- total_cd$accession

# visualization
library(ggplot2)
library(viridis)
factor(percent_col$accession, levels = pc_order) -> percent_col$accession
#factor(percent_col$accession, levels = cd_order) -> percent_col$accession

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


save.image(file = "GWAS_2022_inference_data_wrangling_EDA.RData")

# block effect
model1 <- lm(log(countDensity)~accession+block+annotations+annotations*block, data = count_density_blk_anno)
a <- aov(model1)
tukey_result <- TukeyHSD(a)
