library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(R.utils)
library(lmerTest)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)

# Must load RDS file
load(file = "data/workspace_modeltable.rds")

my_colors <- RColorBrewer::brewer.pal(4, "Set2")[2:4]
# First figure in the data visualization model.R is using the col.ll data 
colors <- c("Top" = "#FC8D62", "Middle" = "#8DA0CB", "Bottom" = "#E78AC3")

col.class.count2.ms <- col.class.count2 %>%
  dplyr::mutate(region = str_replace_all(string = region, pattern = "TOP",replacement = "Top")) %>%
  dplyr::mutate(region = str_replace_all(string = region, pattern = "MID",replacement = "Middle")) %>%
  dplyr::mutate(region = str_replace_all(string = region, pattern = "BOT",replacement = "Bottom")) %>%
  dplyr::mutate(region = factor(region, levels = c("Top", "Middle", "Bottom")))

p <- ggplot(col.class.count2.ms, aes(x=region, y=count_density^(1/4), fill = region)) +
  geom_boxplot() +
  geom_jitter(size = 0.2, width =0.2)+
  labs(fill = "Region") +
  xlab("") +
  ylab(expression("Count Density"^(1/4))) +
  scale_fill_manual(values = colors) +
  facet_grid(~class) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dat <- as.data.frame(p$data)
data.table::fwrite(file = "results/Figure_7_panel_4.txt", x= dat)
ggsave(filename = "results/Figure_7_panel_4.pdf", plot = p, device = "pdf", width = 6, height = 6, units = 'in')


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
#   [1] ggpubr_0.6.0       RColorBrewer_1.1-3 R.utils_2.12.3     R.oo_1.25.0        R.methodsS3_1.8.2  lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1     
# [9] purrr_1.0.2        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1       tidyverse_2.0.0    dplyr_1.1.4        ggplot2_3.4.4      gridExtra_2.3     
# [17] lmerTest_3.1-3     lme4_1.1-35.1      Matrix_1.6-1.1    
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4        xfun_0.41           rstatix_0.7.2       lattice_0.21-9      tzdb_0.4.0          numDeriv_2016.8-1.1 vctrs_0.6.4         tools_4.3.2        
# [9] generics_0.1.3      fansi_1.0.5         pkgconfig_2.0.3     data.table_1.14.8   lifecycle_1.0.4     compiler_4.3.2      farver_2.1.1        textshaping_0.3.7  
# [17] munsell_0.5.0       carData_3.0-5       htmltools_0.5.7     yaml_2.3.7          pillar_1.9.0        car_3.1-2           nloptr_2.0.3        MASS_7.3-60        
# [25] boot_1.3-28.1       abind_1.4-5         nlme_3.1-163        tidyselect_1.2.0    digest_0.6.33       stringi_1.8.2       labeling_0.4.3      splines_4.3.2      
# [33] fastmap_1.1.1       grid_4.3.2          colorspace_2.1-0    cli_3.6.1           magrittr_2.0.3      utf8_1.2.4          broom_1.0.5         withr_2.5.2        
# [41] scales_1.2.1        backports_1.4.1     timechange_0.2.0    rmarkdown_2.25      ggsignif_0.6.4      ragg_1.2.6          hms_1.1.3           evaluate_0.23      
# [49] knitr_1.45          viridisLite_0.4.2   rlang_1.1.2         Rcpp_1.0.11         glue_1.6.2          rstudioapi_0.15.0   minqa_1.2.6         R6_2.5.1           
# [57] systemfonts_1.0.5

