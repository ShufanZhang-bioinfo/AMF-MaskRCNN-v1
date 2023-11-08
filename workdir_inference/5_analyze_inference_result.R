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
##              Data wrangling: calculate overall colonization
##=============================================================================================================
# calculate fungal structure frequency by accession and block
#seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
#  group_by(accession, block, annotations) %>% 
#  summarise(freq = n()) -> fungal_frq

seg_clean_v3[!seg_clean_v3$annotations %in% c("root", "AMF spore", "AMF vesicle"),] %>% 
  group_by(accession, block, annotations) %>% 
  summarise(freq = n()) -> fungal_frq
fungal_frq$key <- paste(fungal_frq$accession, fungal_frq$block, sep = "_")
spread(fungal_frq, key = "annotations",value = "freq") -> fungal_frq_wide
fungal_frq_wide[is.na(fungal_frq_wide)] <- 0

# calculate overall fungal frequency by accession and block
seg_clean_v3[!seg_clean_v3$annotations %in% c("root", "AMF spore", "AMF vesicle"),] %>% group_by(accession, block) %>% summarise(freq = n()) -> fungal_frq_overall
fungal_frq_overall$key <- paste(fungal_frq_overall$accession, fungal_frq_overall$block, sep = "_")

# calculate total root area per accession and block
seg_clean_v3[seg_clean_v3$annotations == "root",] %>% group_by(accession,block) %>% summarise(rootArea = sum(area)) -> rootarea
rootarea$key <- paste(rootarea$accession, rootarea$block, sep = "_")
rootarea$rootArea_mm2 <- rootarea$rootArea * 0.1202399 *10^-6

# calculate overall count density by accession and block
merge(fungal_frq_overall, rootarea[,-c(1:3)], by = "key", all.x = TRUE) -> overall_count_density_blk
overall_count_density_blk$countDensity <- overall_count_density_blk$freq / overall_count_density_blk$rootArea_mm2
merge(overall_count_density_blk, fungal_frq_wide[,-c(1,2)], by="key",all.x = TRUE) -> overall_count_density_blk 
#colnames(overall_count_density_blk) <- c("key","accession","block",             
#                                       "fungalCount","rootArea","rootArea_mm2",      
#                                       "countDensity","ves_count","sp_count",         
#                                       "inH_count", "exH_count", "arb_count")

colnames(overall_count_density_blk) <- c("key","accession","block",             
                                         "fungalCount","rootArea_mm2",      
                                         "countDensity",        
                                         "inH_count", "exH_count", "arb_count")
#===================================================================
# calculate fungal size by accession and block

#seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
#  group_by(accession, block, annotations) %>% 
#  summarise(fungalArea = sum(area_mm2)) -> fungal_area

seg_clean_v3[!seg_clean_v3$annotations %in% c("root", "AMF spore", "AMF vesicle"),] %>% 
  group_by(accession, block, annotations) %>% 
  summarise(fungalArea = sum(area_mm2)) -> fungal_area
fungal_area$key <- paste(fungal_area$accession, fungal_area$block, sep = "_")
spread(fungal_area, key = "annotations",value = "fungalArea") -> fungal_area_wide
fungal_area_wide[is.na(fungal_area_wide)] <- 0

# calculate overall fungal area by accession and block
#seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
#  group_by(accession, block) %>% 
#  summarise(fungalArea = sum(area_mm2)) -> fungal_area_overall

seg_clean_v3[!seg_clean_v3$annotations %in% c("root", "AMF spore", "AMF vesicle"),] %>% 
  group_by(accession, block) %>% 
  summarise(fungalArea = sum(area_mm2)) -> fungal_area_overall
fungal_area_overall$key <- paste(fungal_area_overall$accession, fungal_area_overall$block, sep = "_")

# calculate overall percent colonization by accession and block
merge(fungal_area_overall, rootarea[,-c(1:3)], by = "key", all.x = TRUE) -> overall_percent_col_blk
overall_percent_col_blk$percentCol <- overall_percent_col_blk$fungalArea / overall_percent_col_blk$rootArea_mm2
merge(overall_percent_col_blk, fungal_area_wide[,-c(1,2)], by="key",all.x = TRUE) -> overall_percent_col_blk 
#colnames(overall_percent_col_blk) <- c("key","accession","block",             
#                                       "fungalArea","rootArea","rootArea_mm2",      
#                                       "percentCol","ves_area_mm2","sp_area_mm2",         
#                                       "inH_area_mm2", "exH_area_mm2", "arb_area_mm2")

colnames(overall_percent_col_blk) <- c("key","accession","block",             
                                       "fungalArea","rootArea_mm2",      
                                       "percentCol",        
                                       "inH_area_mm2", "exH_area_mm2", "arb_area_mm2")

#===================================================================
# all phenotypes in one dataframe
merge(overall_percent_col_blk, overall_count_density_blk[,-c(2,3,5)], by = "key") -> pheno

##=============================================================================================================
##              Data wrangling: calculate class-level colonization
##=============================================================================================================

# calculate class level count density
merge(fungal_frq, rootarea[,-c(1:3)], by = "key", all.x = TRUE) -> cd_blk_anno
cd_blk_anno$countDensity <- cd_blk_anno$freq / cd_blk_anno$rootArea_mm2
cd_blk_anno$tf_countDensity <- cd_blk_anno$countDensity^0.2

# calculate class level percent colonization
merge(fungal_area, rootarea[,-c(1:3)], by = "key", all.x = TRUE) -> pc_blk_anno
pc_blk_anno$percentCol <- pc_blk_anno$fungalArea / pc_blk_anno$rootArea_mm2
pc_blk_anno$tf_percentCol <- pc_blk_anno$percentCol^0.2

##=============================================================================================================
##              Data wrangling: total colonization by accession
##=============================================================================================================
# total percent colonization

#seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
#  group_by(accession) %>% 
#  summarise(sumFungalArea = sum(area_mm2)) -> totalfungalarea

seg_clean_v3[!seg_clean_v3$annotations %in% c("AMF spore", "AMF vesicle", "root"),] %>% 
  group_by(accession) %>% 
  summarise(sumFungalArea = sum(area_mm2)) -> totalfungalarea

seg_clean_v3[seg_clean_v3$annotations == "root",] %>% 
  group_by(accession) %>% 
  summarise(sumRootArea = sum(area_mm2)) -> totalrootarea

merge(totalfungalarea, totalrootarea, by = "accession", all.x = TRUE) -> total_pc
total_pc$percentCol <- total_pc$sumFungalArea / total_pc$sumRootArea

# total count density
#seg_clean_v3[!seg_clean_v3$annotations == "root",] %>% 
#  group_by(accession) %>% 
#  summarise(sumFungalCount = n()) -> totalfungalcount

seg_clean_v3[!seg_clean_v3$annotations %in% c("AMF spore", "AMF vesicle", "root"),] %>% 
  group_by(accession) %>% 
  summarise(sumFungalCount = n()) -> totalfungalcount

merge(totalfungalcount, totalrootarea, by = "accession", all.x = TRUE) -> total_cd
total_cd$countDensity <- total_cd$sumFungalCount / total_cd$sumRootArea


##=============================================================================================================
##              Order accessions
##=============================================================================================================

# order accessions by percent colonization
total_pc[order(total_pc$percentCol, decreasing = TRUE),] -> total_pc
order_pc <- total_pc$accession

# order accessions by count density
total_cd[order(total_cd$countDensity, decreasing = TRUE),] -> total_cd
order_cd <- total_cd$accession

# for each block, order accession by percent colonization
pc_blk_anno %>% 
  group_by(accession, block) %>% 
  summarise(overallPercentCol = sum(percentCol)) -> pc_blk

pc_blk[pc_blk$block == "Block2",] -> tmp
tmp$trans_percentCol <- tmp$overallPercentCol^0.2
tmp[order(tmp$overallPercentCol,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_pc_blk2

overall_percent_col_blk[overall_percent_col_blk$block == "Block2",] -> tmp 
tmp[order(tmp$percentCol,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_pc_blk2

overall_percent_col_blk[overall_percent_col_blk$block == "Block8",] -> tmp 
tmp[order(tmp$percentCol,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_pc_blk8

overall_percent_col_blk[overall_percent_col_blk$block == "Block10",] -> tmp 
tmp[order(tmp$percentCol,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_pc_blk10

cd_blk_anno %>% 
  group_by(accession, block) %>% 
  summarise(overallCountDensity = sum(countDensity)) -> cd_blk

cd_blk[cd_blk$block == "Block2",] -> tmp 
tmp[order(tmp$overallCountDensity,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_cd_blk2

pheno[pheno$block == "Block2",] -> tmp 
tmp[order(tmp$countDensity,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_cd_blk2

pheno[pheno$block == "Block8",] -> tmp 
tmp[order(tmp$countDensity,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_cd_blk8

pheno[pheno$block == "Block10",] -> tmp 
tmp[order(tmp$countDensity,decreasing = TRUE), ] -> tmp  
tmp$accession -> order_cd_blk10
##=============================================================================================================
##              Order Blocks
##=============================================================================================================
factor(seg_clean_v3$block, levels = c("Block2", "Block8", "Block10")) -> seg_clean_v3$block
factor(pheno$block, levels = c("Block2", "Block8", "Block10")) -> pheno$block
factor(cd_blk_anno$block, levels = c("Block2", "Block8", "Block10")) -> cd_blk_anno$block
factor(pc_blk_anno$block, levels = c("Block2", "Block8", "Block10")) -> pc_blk_anno$block
factor(overall_count_density_blk$block, levels = c("Block2", "Block8", "Block10")) -> overall_count_density_blk$block
factor(overall_percent_col_blk$block, levels = c("Block2", "Block8", "Block10")) -> overall_percent_col_blk$block

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


