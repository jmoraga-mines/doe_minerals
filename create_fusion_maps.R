if(!require("pacman")){
  install.packages("pacman")
  require("pacman")
}
pacman::p_load(raster, rgdal)

source("utils/doe_mineral_utils.R")

orange2red <- (rev(heat.colors(12)))[6:12]
layer_names <- c("Chalcedony","Kaolinite", "Gypsum", "Hematite", "Epsomite")


s1 <- stack("../../doe_test/HyMapTargetDetectionFull/HyMap_full_mtmf")
s1
s2 <- stack("../../doe_test/HyMapTargetDetectionFull/HyMap_full_mttcimf")
s2
names(s1)
names(s2)
i1_chalcedony <- s1[[14]]
i2_chalcedony <- s2[[14]]
m1_chalcedony <- s1[[07]]
m2_chalcedony <- s2[[07]]
i1_chalcedony <- raster::setMinMax(i1_chalcedony)
i2_chalcedony <- raster::setMinMax(i2_chalcedony)
m1_chalcedony <- raster::setMinMax(m1_chalcedony)
m2_chalcedony <- raster::setMinMax(m2_chalcedony)
chalcedony_1 <- m1_chalcedony - i1_chalcedony
chalcedony_1 <- raster::setMinMax(chalcedony_1)

m1 <- m1_chalcedony  # use m1_chalcedony or abs(m1_chalcedony)
m1 <- m1/raster::maxValue(m1)
i1 <- i1_chalcedony/raster::maxValue(i1_chalcedony)
x1 <- m1>=(i1*1)
m1 <- m1*x1
#m1[i1_chalcedony<18.12090]<- NA
m1[m1<0.1]<- NA

raster::plot(m1, col=orange2red)
r1 <- raster::crop(m1, extent_tall_brady)
raster::plot(r1, col=orange2red)
rx1 <- raster::crop(x1, extent_tall_brady)
raster::plot(r1*rx1, col=orange2red)
r1_i <- raster::crop(i1_chalcedony, extent_tall_brady)
raster::plot(r1_i)



hymap <- stack("D:/HyMap_USA_2003_B_NV_HyVista/Brady_ref_L1C_mos.bil")
brady <- raster::crop(hymap, extent_tall_brady)


################ New Data ################

ace_th <- c(0.09,	0.17,	0.1, 0.13, 0.09)
mf_th <- c(0.04, 0.04, 0.023, 0.03, 0.009)
mtmf_th <- c(0.025, 0.022, 0.01, 0.015, 0.005)
mtmf_th <- mf_th*0.8

cem_th <- c(0.04, 0.04, 0.023, 0.03, 0.009)
osp_th <- c(0.14, 0.06, 0.06, 0.08, 0.02)
sam_th <- c(1, 1, 1, 1, 1)
tcimf_th <- c(0.07, 0.055, 0.040, 0.055, 0.015)
mttcimf_th <- tcimf_th*0.8



base_directory <- "D:/mm_final/HyMap"
analysis_file <- "HyMapFull_tcimf"
detection_file <- file.path(base_directory, analysis_file)
minerals_stack <- stack(detection_file)
names(minerals_stack)<-gsub(".*\\.(.*)\\.", "\\1", names(minerals_stack))

brady<-crop(minerals_stack, extent_tall_brady)
brady[brady<=0]<-NA
raster::spplot(brady[["Chalcedony"]])
b<-brady
b[["Chalcedony"]] <- b[["Chalcedony"]]*2
fusion <- raster::calc(b, fun=sum, na.rm=TRUE)

ogrListLayers("d:/Jim+Erika/BuiltUp_Brady_Desert.shp")
ogrInfo("d:/Jim+Erika/BuiltUp_Brady_Desert.shp")
BuildUp <- readOGR("d:/Jim+Erika/BuiltUp_Brady_Desert.shp", integer64 = "allow.loss", layer="BuiltUp_Brady_Desert")
BuildUp_id <- BuildUp["id"]
BuildUp_id <- spTransform(BuildUp_id, crs(hymap))
spplot(BuildUp_id)

# brady<-stack("d:/mm_final/HyMap")
BuiltUpRaster <- rasterize(BuildUp_id, hymap[[1]], na.rm=TRUE)
spplot(BuiltUpRaster)
spplot(BuildUp)
BuildUp
b <- BuiltUpRaster
b[!is.na(b)] <- 1
b[is.na(b)] <- 0
c <- !b
c
spplot(c)
c_brady <- crop(c, extent_tall_brady)
spplot(c_brady)
f <- fusion*c
f <- fusion*c_brady
f[f<=0] <- NA
spplot(f)
f_zero <- f
f_zero[f<0.1] <- NA
spplot(f_zero)

k <- stack("d:/mm_final/HyMap/HyMapFull_tcimf_target")
k <- stack("d:/mm_final/HyMap/HyMapFull_tcimf")
# names(k)<-gsub("(.*[\\.]+)+(.*)\\.\\.", "\\2", names(k))
names(k) <- layer_names
chalcedony <- k[["Chalcedony"]]
chalcedony[chalcedony<0.07]<-NA
chalcedony_brady <- crop(chalcedony, extent_tall_brady)
k_brady<-crop(k, extent_tall_brady)
k_brady[is.na(k_brady)]<-0
k_brady[k_brady<c(0.07, 0.055, 0.040, 0.055, 0.015)] <- NA
spplot(k_brady)
k_brady <- calc(k_brady, function(x){return (x*(x>c(0.07, 0.055, 0.040, 0.055, 0.02)))})
spplot(k_brady)
k_brady_fusion <- calc(k_brady, fun = sum, rm.na=TRUE)
k_brady_fusion[k_brady_fusion<=0] <- NA

spplot(k_brady_fusion)

names(k_brady_fusion) <- "tcimf_fusion"
f1<-doe_write_raster(k_brady_fusion, "d:/mm_final/Fusion/tcimf_brady_fusion")


sam_rule <- stack("d:/mm_final/HyMap/HyMapFull_sam_rule")
spplot(sam_rule)
names(sam_rule) <- layer_names
sam_rule_brady <- crop(sam_rule, extent_tall_brady)
sam_rule_brady[sam_rule_brady==0] <- NA
sam_rule_brady <- (pi-abs(sam_rule_brady))/pi
sam_rule_brady[sam_rule_brady<c(0.62, 0.60, 0.61, 0.58,0.60)] <- NA
names(sam_rule_brady) <- layer_names
spplot(sam_rule_brady)

doe_write_raster(sam_rule_brady, "d:/mm_final/Fusion/sam_rule_brady")

sam_rule_brady_fusion <- calc(sam_rule_brady, fun=sum, na.rm=TRUE)
names(sam_rule_brady_fusion) <- "sam_rule_fusion"
spplot(sam_rule_brady_fusion)
f1<-doe_write_raster(sam_rule_brady_fusion, "d:/mm_final/Fusion/sam_brady_fusion")


mf <- stack("d:/mm_final/HyMap/HyMapFull_mf")
#spplot(mf)
mf_brady <- crop(mf, extent_tall_brady)
# mf_brady <- calc(mf_brady, function(x){return (x*(x>mf_th))})
mf_brady[mf_brady<mf_th] <- NA
names(mf_brady) <- layer_names
spplot(mf_brady)
mf_brady_fusion <- calc(mf_brady, fun=sum, na.rm=TRUE)
names(mf_brady_fusion) <- "mf_fusion"
spplot(mf_brady_fusion)
f1<-doe_write_raster(mf_brady_fusion, "d:/mm_final/Fusion/mf_brady_fusion")

cem <- stack("d:/mm_final/HyMap/HyMapFull_cem")
#spplot(cem)
cem_brady <- crop(cem, extent_tall_brady)
# cem_brady <- calc(cem_brady, function(x){return (x*(x>cem_th))})
cem_brady[cem_brady<cem_th] <- NA
names(cem_brady) <- layer_names
spplot(cem_brady)
cem_brady_fusion <- calc(cem_brady, fun=sum, na.rm=TRUE)
names(cem_brady_fusion) <- "cem_fusion"
spplot(cem_brady_fusion)
f1<-doe_write_raster(cem_brady_fusion, "d:/mm_final/Fusion/cem_brady_fusion")


osp <- stack("d:/mm_final/HyMap/HyMapFull_osp")
#spplot(osp)
osp_brady <- crop(osp, extent_tall_brady)
# osp_brady <- calc(osp_brady, function(x){return (x*(x>osp_th))})
osp_brady[osp_brady<osp_th] <- NA
names(osp_brady) <- layer_names
spplot(osp_brady)
osp_brady_fusion <- calc(osp_brady, fun=sum, na.rm=TRUE)
names(osp_brady_fusion) <- "osp_fusion"
spplot(osp_brady_fusion)
f1<-doe_write_raster(osp_brady_fusion, "d:/mm_final/Fusion/osp_brady_fusion")


ace <- stack("d:/mm_final/HyMap/HyMapFull_ace")
#spplot(ace)
ace_brady <- crop(ace, extent_tall_brady)
# ace_brady <- calc(ace_brady, function(x){return (x*(x>ace_th))})
ace_brady[ace_brady<ace_th] <- NA
names(ace_brady) <- layer_names
spplot(ace_brady)
ace_brady_fusion <- calc(ace_brady, fun=sum, na.rm=TRUE)
names(ace_brady_fusion) <- "ace_fusion"
spplot(ace_brady_fusion)
f1<-doe_write_raster(ace_brady_fusion, "d:/mm_final/Fusion/ace_brady_fusion")


tcimf <- stack("d:/mm_final/HyMap/HyMapFull_tcimf")
#spplot(tcimf)
tcimf_brady <- crop(tcimf, extent_tall_brady)
# tcimf_brady <- calc(tcimf_brady, function(x){return (x*(x>tcimf_th))})
tcimf_brady[tcimf_brady<tcimf_th] <- NA
names(tcimf_brady) <- layer_names
spplot(tcimf_brady)
tcimf_brady_fusion <- calc(tcimf_brady, fun=sum, na.rm=TRUE)
names(tcimf_brady_fusion) <- "tcimf_fusion"
spplot(tcimf_brady_fusion)
f1<-doe_write_raster(tcimf_brady_fusion, "d:/mm_final/Fusion/tcimf_brady_fusion")


################### With Infeasibility #######
mtmf <- stack("d:/mm_final/HyMap/HyMapFull_mtmf")
#spplot(mtmf)
mtmf_brady <- crop(mtmf, extent_tall_brady)
# mtmf_brady <- calc(mtmf_brady, function(x){return (x*(x>mtmf_th))})

mtmf_brady_infeasibility <- mtmf_brady[[6:10]]
mtmf_brady_mf <- mtmf_brady[[1:5]]
names(mtmf_brady_mf) <- layer_names
names(mtmf_brady_infeasibility) <- layer_names

infeasibility_curve <- (mtmf_brady_mf*300 + 15 - mtmf_brady_infeasibility)>0
infeasibility_curve
mtmf_brady <- mtmf_brady_mf*infeasibility_curve
mtmf_brady[mtmf_brady<(mtmf_th)] <- NA
spplot(stack(mf_brady[[1]], mtmf_brady[[1]]), col=orange2red)
spplot(stack(mf_brady, mtmf_brady), col=orange2red)
names(mtmf_brady) <- layer_names
spplot(mtmf_brady)
mtmf_brady_fusion <- calc(mtmf_brady, fun=sum, na.rm=TRUE)
names(mtmf_brady_fusion) <- "mtmf_fusion"
spplot(mtmf_brady_fusion)
f1<-doe_write_raster(mtmf_brady_fusion, "d:/mm_final/Fusion/mtmf_brady_fusion")


###### mttcimf 
mttcimf <- stack("d:/mm_final/HyMap/HyMapFull_mttcimf")
#spplot(mttcimf)
mttcimf_brady <- crop(mttcimf, extent_tall_brady)
# mttcimf_brady <- calc(mttcimf_brady, function(x){return (x*(x>mttcimf_th))})

mttcimf_brady_infeasibility <- mttcimf_brady[[6:10]]
mttcimf_brady_mf <- mttcimf_brady[[1:5]]
names(mttcimf_brady_mf) <- layer_names
names(mttcimf_brady_infeasibility) <- layer_names

infeasibility_curve <- (mttcimf_brady_mf*300 + 15 - mttcimf_brady_infeasibility)>0
infeasibility_curve
mttcimf_brady <- mttcimf_brady_mf*infeasibility_curve
mttcimf_brady[mttcimf_brady<(mttcimf_th)] <- NA
spplot(stack(mf_brady[[1]], mttcimf_brady[[1]]), col=orange2red)
spplot(stack(mf_brady, mttcimf_brady), col=orange2red)
names(mttcimf_brady) <- layer_names
spplot(mttcimf_brady)
mttcimf_brady_fusion <- calc(mttcimf_brady, fun=sum, na.rm=TRUE)
names(mttcimf_brady_fusion) <- "mttcimf_fusion"
spplot(mttcimf_brady_fusion)
f1<-doe_write_raster(mttcimf_brady_fusion, "d:/mm_final/Fusion/mttcimf_brady_fusion")

###############
break_cuts <- 10
break_points <- quantile(fusion_mtmf_brady,probs=seq(0,1,1/break_cuts),na.rm=TRUE)
spplot(mtmf_brady_fusion, cuts=break_cuts, at=break_points, 
       col.regions=orange2red, colorkey=list(width=0.3, space="right", 
                                             tick.number=break_cuts, 
                                             labels=list(
                                               at=break_points,
                                               labels=break_points)))
