if(!require("pacman")){
  install.packages("pacman")
  require("pacman")
}
pacman::p_load(raster, rgdal)

source("utils/doe_mineral_utils.R")

fusion_files <- list.files(path="d:/mm_final/Fusion", pattern = "fusion.gri$", 
           all.files = TRUE, full.names = TRUE)
fusion_stack <- stack(fusion_files)

fusion_stack <- setMinMax(fusion_stack)
fusion_stack <- fusion_stack/maxValue(fusion_stack)
spplot(fusion_stack)
sam_rule_fusion <- fusion_stack[["sam_rule_fusion"]] 
spplot(sam_rule_fusion)

###### Trying and failing to obtain better results using infeasibility 
mttcimf <- stack("d:/mm_final/HyMap/HyMapFull_mttcimf")
mttcimf_brady <- crop(mttcimf, extent_tall_brady)
mttcimf_brady_infeasibility <- mttcimf_brady[[6:10]]
names(mttcimf_brady_infeasibility) <- layer_names
mttcimf_brady_infeasibility <- mttcimf_brady_infeasibility[[1]]
infeasibility_curve <- (sam_rule_fusion*100 - mttcimf_brady_infeasibility)>0
infeasibility_curve <- mttcimf_brady_infeasibility < 25
sam_fusion <- sam_rule_fusion*infeasibility_curve
spplot(sam_fusion)
