library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(R.utils)
library(lmerTest)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(emmeans)
library(ggheatmap)
library(reshape2)
library(magrittr)
library(multipanelfigure)
getwd()
setwd("/Users/lovely_shufan/Downloads/workdir 2/")
# Must load RDS file
load(file = "/Users/lovely_shufan/Dropbox (Edison_Lab@UGA)/AMF/AMF Imaging 2021/3_analysis/inference_on_14_RILS_accessions/inference using Georgia model/workspace_modeltable.rds")


my_colors <- RColorBrewer::brewer.pal(n = 7, name = "Dark2")
# First figure in the data visualization model.R is using the col.ll data 

colors <- c("Extraradical hypha" = "#1B9E77", "Intraradical hypha" = "#7570B3", "Arbuscule" = "#66A61E", "Vesicle" = "#E6AB02", "Spore" = "#A6761D", "Root" = "#D95F02") 

seg.test.ms <- seg_testset %>%
  dplyr::filter(annotations != "others") %>%
  #dplyr::filter(class != "root") %>%
  dplyr::mutate(annotations = str_replace_all(string = annotations, pattern = "AMF external hypha", replacement = "Extraradical hypha")) %>%
  dplyr::mutate(annotations = str_replace_all(string = annotations, pattern = "AMF internal hypha", replacement = "Intraradical hypha")) %>%
  dplyr::mutate(annotations = str_replace_all(string = annotations, pattern = "AMF arbuscule", replacement = "Arbuscule")) %>%
  dplyr::mutate(annotations = str_replace_all(string = annotations, pattern = "AMF vesicle", replacement = "Vesicle")) %>%
  dplyr::mutate(annotations = str_replace_all(string = annotations, pattern = "AMF spore", replacement = "Spore")) %>%
  dplyr::mutate(annotations = str_replace_all(string = annotations, pattern = "root", replacement = "Root"))

seg.ms <- seg2 %>%
  dplyr::filter(class != "others") %>%
  #dplyr::filter(class != "root") %>%
  dplyr::mutate(class = str_replace_all(string = class, pattern = "AMF external hypha", replacement = "Extraradical hypha")) %>%
  dplyr::mutate(class = str_replace_all(string = class, pattern = "AMF internal hypha", replacement = "Intraradical hypha")) %>%
  dplyr::mutate(class = str_replace_all(string = class, pattern = "AMF arbuscule", replacement = "Arbuscule")) %>%
  dplyr::mutate(class = str_replace_all(string = class, pattern = "AMF vesicle", replacement = "Vesicle")) %>%
  dplyr::mutate(class = str_replace_all(string = class, pattern = "AMF spore", replacement = "Spore")) %>%
  dplyr::mutate(class = str_replace_all(string = class, pattern = "root", replacement = "Root")) %>%
  dplyr::mutate(region = str_replace_all(string = region, pattern = "TOP",replacement = "Top")) %>%
  dplyr::mutate(region = str_replace_all(string = region, pattern = "MID",replacement = "Middle")) %>%
  dplyr::mutate(region = str_replace_all(string = region, pattern = "BOT",replacement = "Bottom")) %>%
  dplyr::mutate(region = factor(region, levels = c("Top", "Middle", "Bottom")))

seg.test.ms$annotations <- factor(seg.test.ms$annotations, levels = c("Arbuscule" ,"Extraradical hypha", "Intraradical hypha",        
                                           "Spore","Vesicle", "Root" ))

p1 <- ggplot(data = seg.test.ms, aes(x = confidenceScore, fill = annotations)) +
  geom_density() +
  facet_wrap(~annotations, ncol = 6) +
  scale_fill_manual(values = colors) +
  xlab("Confidence score") +
  ylab("Density") +
  theme_bw()+
  theme(legend.position = "none")

# Visualization
dat1 <- as.data.frame(p1$data) %>%
  dplyr::mutate(Panel = "P1")

# ouputs
if(all(colnames(dat1) == colnames(dat2))){
  dat <- rbind(dat1,dat2)
}
data.table::fwrite(file = "results/Figure_2_panel_1_2.txt", x= dat)
ggsave(filename = "~/Dropbox (Edison_Lab@UGA)/AMF/AMF Imaging 2021/4_writing/Nature Scientific Report/Figure_2.pdf", plot = p1, device = "pdf", width = 12, height = 6, units = 'in')


# R version 4.3.2 (2023-10-31)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.7.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] mbr_0.0.1              GGally_2.2.0           multipanelfigure_2.1.5 magrittr_2.0.3         reshape_0.8.9          ggheatmap_2.2          emmeans_1.9.0          ggpubr_0.6.0           RColorBrewer_1.1-3    
# [10] gridExtra_2.3          R.utils_2.12.3         R.oo_1.25.0            R.methodsS3_1.8.2      lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1          purrr_1.0.2            readr_2.1.4           
# [19] tidyr_1.3.0            tibble_3.2.1           ggplot2_3.4.4          tidyverse_2.0.0        dplyr_1.1.4            lmerTest_3.1-3         lme4_1.1-35.1          Matrix_1.6-1.1        
# 
# loaded via a namespace (and not attached):
#   [1] rlang_1.1.2             matrixStats_1.1.0       compiler_4.3.2          png_0.1-8               systemfonts_1.0.5       vctrs_0.6.4             RcppZiggurat_0.1.6      pkgconfig_2.0.3         crayon_1.5.2           
# [10] fastmap_1.1.1           backports_1.4.1         magick_2.8.2            labeling_0.4.3          utf8_1.2.4              rmarkdown_2.25          tzdb_0.4.0              nloptr_2.0.3            ragg_1.2.6             
# [19] xfun_0.41               Rfast_2.1.0             cachem_1.0.8            aplot_0.2.2             progress_1.2.2          broom_1.0.5             parallel_4.3.2          prettyunits_1.2.0       R6_2.5.1               
# [28] stringi_1.8.2           car_3.1-2               boot_1.3-28.1           numDeriv_2016.8-1.1     estimability_1.4.1      Rcpp_1.0.11             knitr_1.45              splines_4.3.2           timechange_0.2.0       
# [37] tidyselect_1.2.0        rstudioapi_0.15.0       abind_1.4-5             yaml_2.3.7              lattice_0.21-9          plyr_1.8.9              withr_2.5.2             signal_1.8-0            evaluate_0.23          
# [46] gridGraphics_0.5-1      RcppParallel_5.1.7      ggstats_0.5.1           pillar_1.9.0            carData_3.0-5           ggfun_0.1.3             generics_0.1.3          hms_1.1.3               assertive.numbers_0.0-2
# [55] munsell_0.5.0           scales_1.2.1            minqa_1.2.6             glue_1.6.2              tools_4.3.2             dplR_1.7.6              data.table_1.14.8       assertive.files_0.0-2   ggsignif_0.6.4         
# [64] XML_3.99-0.15           fs_1.6.3                mvtnorm_1.2-4           cowplot_1.1.2           grid_4.3.2              colorspace_2.1-0        nlme_3.1-163            patchwork_1.1.3         cli_3.6.1              
# [73] textshaping_0.3.7       fansi_1.0.5             gtable_0.3.4            rstatix_0.7.2           yulab.utils_0.1.1       digest_0.6.33           assertive.base_0.0-9    ggrepel_0.9.4           ggplotify_0.1.2        
# [82] farver_2.1.1            memoise_2.0.1           htmltools_0.5.7         factoextra_1.0.7        lifecycle_1.0.4         MASS_7.3-60 
# 
# 
# 
