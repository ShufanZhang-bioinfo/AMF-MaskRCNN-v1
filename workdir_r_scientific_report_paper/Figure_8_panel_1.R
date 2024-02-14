library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(R.utils)
library(lmerTest)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(emmeans)


# Must load RDS file
load(file = "data/workspace_modeltable.rds")

colors <- c("Arbuscule"= "#084594", "External Hyphae" = "#2B8CBE", "Spore" = "#4DAF4A", "Vesicle"= "#FFFF33")

col.class.count3.sub.ms <- col.class.count3.sub %>%
  dplyr::mutate(class2 = str_replace_all(string = class2, pattern = "arb", replacement = "Arbuscule")) %>%
  dplyr::mutate(class2 = str_replace_all(string = class2, pattern = "exH", replacement = "External Hyphae")) %>%
  dplyr::mutate(class2 = str_replace_all(string = class2, pattern = "sp", replacement = "Spore")) %>%
  dplyr::mutate(class2 = str_replace_all(string = class2, pattern = "ves", replacement = "Vesicle"))

p <- ggplot2::ggplot(col.class.count3.sub.ms, aes(x=factor(accession, levels = accession_order_pc),y=count_density^0.5,fill=class2)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  ylab(expression("Class count density"^(1/2))) +
  xlab("") +
  labs(fill = "AMF Structure") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
dat <- as.data.frame(p$data)
data.table::fwrite(file = "results/Figure_8_panel_1.txt", x= dat)
ggsave(filename = "results/Figure_8_panel_1.pdf", plot = p, device = "pdf", width = 6, height = 6, units = 'in')


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
#   [1] emmeans_1.9.0      ggpubr_0.6.0       RColorBrewer_1.1-3 R.utils_2.12.3     R.oo_1.25.0        R.methodsS3_1.8.2  lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1     
# [10] purrr_1.0.2        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1       tidyverse_2.0.0    dplyr_1.1.4        ggplot2_3.4.4      gridExtra_2.3      lmerTest_3.1-3    
# [19] lme4_1.1-35.1      Matrix_1.6-1.1    
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4        xfun_0.41           rstatix_0.7.2       lattice_0.21-9      tzdb_0.4.0          numDeriv_2016.8-1.1 vctrs_0.6.4         tools_4.3.2        
# [9] generics_0.1.3      fansi_1.0.5         pkgconfig_2.0.3     data.table_1.14.8   lifecycle_1.0.4     compiler_4.3.2      farver_2.1.1        textshaping_0.3.7  
# [17] munsell_0.5.0       carData_3.0-5       htmltools_0.5.7     yaml_2.3.7          pillar_1.9.0        car_3.1-2           nloptr_2.0.3        MASS_7.3-60        
# [25] boot_1.3-28.1       abind_1.4-5         nlme_3.1-163        tidyselect_1.2.0    digest_0.6.33       mvtnorm_1.2-4       stringi_1.8.2       labeling_0.4.3     
# [33] splines_4.3.2       fastmap_1.1.1       grid_4.3.2          colorspace_2.1-0    cli_3.6.1           magrittr_2.0.3      utf8_1.2.4          broom_1.0.5        
# [41] withr_2.5.2         scales_1.2.1        backports_1.4.1     estimability_1.4.1  timechange_0.2.0    rmarkdown_2.25      ggsignif_0.6.4      ragg_1.2.6         
# [49] hms_1.1.3           evaluate_0.23       knitr_1.45          viridisLite_0.4.2   mgcv_1.9-0          rlang_1.1.2         Rcpp_1.0.11         glue_1.6.2         
# [57] rstudioapi_0.15.0   minqa_1.2.6         R6_2.5.1            systemfonts_1.0.5