if(!require("pacman")){
  install.packages("pacman")
  require("pacman")
}
pacman::p_load(raster, rgdal)

source("utils/doe_mineral_utils.R")

orange2red <- (rev(heat.colors(15)))[6:15]
layer_names <- c("Chalcedony","Kaolinite", "Gypsum", "Hematite", "Epsomite")

################ New Data ################

ace_th <- c(0.09,	0.17,	0.1, 0.13, 0.09)
mf_th <- c(0.04, 0.04, 0.023, 0.03, 0.009)
mtmf_th <- c(0.025, 0.022, 0.01, 0.015, 0.005)
mtmf_th <- mf_th*0.8

cem_th <- c(0.04, 0.04, 0.023, 0.03, 0.009)
osp_th <- c(0.14, 0.06, 0.06, 0.08, 0.02)
sam_th <- c(0.65, 0.65, 0.65, 0.65, 0.65)
sam_th2 <- c(pi/2, pi/2, pi/2, pi/2, pi/2)
tcimf_th <- c(0.07, 0.055, 0.040, 0.055, 0.015)
mttcimf_th <- tcimf_th*0.8

base_directory <- "D:/mm_final/HyMap"
analysis_file <- "HyMapFull_tcimf"

mm_plot <- function(x, cuts=10, cut_colors=(rev(heat.colors(15)))[6:15], 
                    plot=TRUE)
{
  break_points <- quantile(x,probs=seq(0,1,1/cuts), na.rm=TRUE)
  p <- spplot(x, cuts=cuts, at=break_points, col.regions=cut_colors,
              colorkey=list(width=0.3, space="right", tick.number=cuts,
                       labels=list(at=break_points,labels=break_points)))
  if(plot==TRUE) plot(p)
  return(p)
}

create_fusion <- function(file_name, new_layer_name, thresholds,
                          layer_names=c("Chalcedony","Kaolinite", "Gypsum", 
                                        "Hematite", "Epsomite"),
                          plot_colors=(rev(heat.colors(12)))[6:12],
                          is_sam=FALSE, has_infeasibility=FALSE,
                          cropping_extent=NULL, plot_crop=FALSE, 
                          output_file_name=NULL, plot_fusion=FALSE,
                          normalize_fusion=FALSE)
{
  s <- stack(file_name)
  if(!exists("s")) return(NULL)
  n<-length(layer_names)
  if(n!=length(thresholds)){
    cat("Length of layer names and thresholds do not match")
    return(NULL)
  }
  if(!is.null(cropping_extent))
    s_crop <- crop(s, cropping_extent)
  else
    s_crop<-s
  l<-nlayers(s_crop)
  if(has_infeasibility){
    if((2*n)!=nlayers(s)){
      cat("Length of layer names and feasibility layers in file do not match\n")
      return(NULL)
    }
    s_feasible <- s_crop[[1:n]]
    s_infeasible <- s_crop[[(n+1):(2*n)]]
    infeasibility_curve <- ((s_feasible*300 + 15 - s_infeasible)>0)
    s_crop <- s_feasible*infeasibility_curve
  } else if(n!=nlayers(s)){
      cat("Length of layer names and layers in file do not match\n")
      return(NULL)
    }
  if(is_sam) {
    cat("SAM analysis:\n")
    s_crop[s_crop==0] <- NA
    print(s_crop)
    s_crop <- (pi-abs(s_crop))/pi
    cat("After applying threshold and pi normalization:\n")
    print(s_crop)
  }
  s_crop[s_crop<thresholds] <- NA
  names(s_crop) <- layer_names
  if(plot_crop) plot(spplot(s_crop, col=plot_colors))
  s_crop_fusion <- calc(s_crop, fun=sum, na.rm=TRUE)
  s_crop_fusion[s_crop_fusion==0] <- NA
  if(normalize_fusion){
    cat("Calculating maximum value to normalize to [0-1]\n")
    cat(paste("MaxValue:", raster::maxValue(s_crop_fusion)),"\n")
    s_crop_max<-raster::maxValue(s_crop_fusion)
    s_crop_fusion <- s_crop_fusion/(s_crop_max)
  }
  names(s_crop_fusion) <- new_layer_name
  if (plot_fusion) plot(spplot(s_crop_fusion, col=plot_colors))
  if (!is.null(output_file_name))
    doe_write_raster(s_crop_fusion, output_file_name)
  return(s_crop_fusion)
}


#################### Desert Peak Maps
fusion_ace_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_ace", "ace_fusion", ace_th,
                                   cropping_extent = extent_desert, 
                                   output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_ace")
fusion_cem_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_cem", "cem_fusion", cem_th,
                                   cropping_extent = extent_desert, 
                                   output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_cem")
fusion_mf_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_mf", "mf_fusion", mf_th,
                           cropping_extent = extent_desert, 
                           output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_mf")
fusion_mtmf_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_mtmf", "mtmf_fusion",
                             mf_th, cropping_extent = extent_desert, 
                             has_infeasibility = T, 
                             output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_mtmf")
fusion_osp_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_osp", "osp_fusion", osp_th,
                                   cropping_extent = extent_desert, 
                                   output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_osp")
fusion_sam_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_sam_rule", "sam_fusion", sam_th,
                                   cropping_extent = extent_desert, is_sam=TRUE, 
                                   plot_crop = TRUE, plot_fusion = TRUE, 
                                   output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_sam")
fusion_tcimf_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_tcimf", "tcimf_fusion", tcimf_th,
                                     cropping_extent = extent_desert, 
                                     output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_tcimf")
fusion_mttcimf_desert <- create_fusion("d:/mm_final/HyMap/HyMapFull_mttcimf", "mttcimf_fusion", mttcimf_th,
                                       cropping_extent = extent_desert, 
                                       has_infeasibility = T,  
                                       output_file_name="d:/mineral_markers_manuscript/fusion/desert_peak/desert_fusion_mttcimf")
############ Brady maps
fusion_ace_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_ace", "ace_fusion", ace_th,
                                  cropping_extent = extent_tall_brady, 
                                  output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_ace")
fusion_cem_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_cem", "cem_fusion", cem_th,
                                  cropping_extent = extent_tall_brady, 
                                  output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_cem")
fusion_mf_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_mf", "mf_fusion", mf_th,
                                 cropping_extent = extent_tall_brady, 
                                 output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_mf")
fusion_mtmf_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_mtmf", "mtmf_fusion",
                                   mf_th, cropping_extent = extent_tall_brady, 
                                   has_infeasibility = T, 
                                   output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_mtmf")
fusion_osp_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_osp", "osp_fusion", osp_th,
                                  cropping_extent = extent_tall_brady, 
                                  output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_osp")
fusion_sam_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_sam_rule", "sam_fusion", sam_th,
                                  cropping_extent = extent_tall_brady, is_sam=TRUE,   
                                  plot_crop = TRUE, plot_fusion = TRUE, 
                                  output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_sam")
fusion_tcimf_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_tcimf", "tcimf_fusion", tcimf_th,
                                    cropping_extent = extent_tall_brady, 
                                    output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_tcimf")
fusion_mttcimf_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_mttcimf", "mttcimf_fusion", mttcimf_th,
                                      cropping_extent = extent_tall_brady, 
                                      has_infeasibility = T,  
                                      output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_mttcimf")
##### Normalization test





fusion_mf_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_mf", "mf_fusion", mf_th,
                           cropping_extent = extent_tall_brady, plot_crop = T, 
                           plot_fusion = T, normalize_fusion = F, 
                           output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_mf")

fusion_mtmf_brady <- create_fusion("d:/mm_final/HyMap/HyMapFull_mtmf", "mtmf_fusion", 
                             mtmf_th, cropping_extent = extent_tall_brady, 
                             plot_crop = T, has_infeasibility = T, 
                             plot_fusion = T, normalize_fusion = F, 
                             output_file_name="d:/mineral_markers_manuscript/fusion/brady/brady_fusion_mtmf")


########### Full image analysis and maps (HyMap for Brady and Desert Peak)
fusion_ace_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_ace", "ace_fusion", ace_th,
                                  output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_ace")
fusion_cem_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_cem", "cem_fusion", cem_th,
                                  output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_cem")
fusion_mf_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_mf", "mf_fusion", mf_th,
                                 output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_mf")
fusion_mtmf_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_mtmf", "mtmf_fusion", mtmf_th,
                                   has_infeasibility = T, 
                                   output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_mtmf")
#fusion_mtmf_hymap2 <- create_fusion("d:/mm_final/HyMap/HyMapFull_mtmf", "mtmf_fusion", mf_th,
#                                   has_infeasibility = T, 
#                                   output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_mtmf2")
fusion_osp_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_osp", "osp_fusion", osp_th,
                                  output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_osp")
fusion_sam_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_sam_rule", "sam_fusion", sam_th, is_sam=TRUE, 
                                  plot_crop = TRUE, plot_fusion = TRUE, 
                                  output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_sam")
fusion_tcimf_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_tcimf", "tcimf_fusion", tcimf_th,
                                    output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_tcimf")
fusion_mttcimf_hymap <- create_fusion("d:/mm_final/HyMap/HyMapFull_mttcimf", "mttcimf_fusion", mttcimf_th,
                                      has_infeasibility = T,  
                                      output_file_name="d:/mineral_markers_manuscript/fusion/hymap/hymap_fusion_mttcimf")




########### Plotting test (do not run)
break_cuts <- 10
break_points <- quantile(fusion_mtmf_brady,probs=seq(0,1,1/break_cuts),na.rm=TRUE)
spplot(fusion_mtmf_brady, cuts=break_cuts, at=break_points, 
       col.regions=orange2red, colorkey=list(width=0.3, space="right", 
                                             tick.number=break_cuts, 
                                             labels=list(
                                               at=break_points,
                                               labels=break_points)))
p<-mm_plot(fusion_mf_brady, cuts = 5, (rev(heat.colors(10)))[6:10])
plot(p)
spplot(fusion_mf_brady)
quantile(fusion_mf_brady, probs=seq(0,1,0.2), na.rm=TRUE)




##### Gather all fused results

brady_files <- list.files(path = "d:/mineral_markers_manuscript/fusion/brady", 
                          pattern = ".gri$", full.names = TRUE)
brady_stack <- stack(brady_files)
doe_write_raster(brady_stack, "d:/mineral_markers_manuscript/fusion/brady_fusion_all")
min_brady <-raster::minValue(brady_stack)
max_brady <-raster::maxValue(brady_stack)
brady_stack <- (brady_stack-min_brady)/(max_brady - min_brady)
min_brady <- max(raster::minValue(brady_stack))
brady_stack[brady_stack<=min_brady]<-NA
doe_write_raster(brady_stack, "d:/mineral_markers_manuscript/fusion/brady_fusion_all_normal")
#b <- brady_stack[[c(1,3:8)]]
b <- brady_stack[[c("mtmf_fusion", "mttcimf_fusion", "osp_fusion", "sam_fusion")]]
b <- raster::calc(b, fun = sum, na.rm=TRUE)
b[b<=min_brady] <- NA
b_windsor <- raster::quantile(b, na.rm=T, probs=c(0.005, 0.995))
b_max <- max(b_windsor) # Winsorizing to percentile 99
b_min <- min(b_windsor) # Winsorizing to percentile 99
b[b>=b_max]<-b_max
b[b>=b_max]<-b_max

b[b>1]<-1
names(b)<-"brady_fusion"
doe_write_raster(b, "d:/mineral_markers_manuscript/fusion/brady_fusion_final")
spplot(b)

#### ####
desert_files <- list.files(path = "d:/mineral_markers_manuscript/fusion/desert_peak", 
                          pattern = ".gri$", full.names = TRUE)
desert_stack <- stack(desert_files)
doe_write_raster(desert_stack, "d:/mineral_markers_manuscript/fusion/desert_fusion_all")

min_desert <-raster::minValue(desert_stack)
max_desert <-raster::maxValue(desert_stack)
desert_stack <- (desert_stack-min_desert)/(max_desert - min_desert)
min_desert <- max(raster::minValue(desert_stack))
desert_stack[desert_stack<=min_desert]<-NA
doe_write_raster(desert_stack, "d:/mineral_markers_manuscript/fusion/desert_fusion_all_normal")
#d <- desert_stack[[1,3:8]]
d <- desert_stack[[c("mtmf_fusion", "mttcimf_fusion", "osp_fusion", "sam_fusion")]]
d <- raster::calc(d, fun = sum, na.rm=TRUE)
d_max <- max(raster::quantile(d, na.rm=T, probs=c(0.99))) # Winsorizing to percentile 99
d <- d/d_max
d[d>1]<-1
names(d)<-"desert_peak_fusion"
doe_write_raster(d, "d:/mineral_markers_manuscript/fusion/desert_fusion_final")
spplot(d)


#### ####
hymap_files <- list.files(path = "d:/mineral_markers_manuscript/fusion/hymap", 
                          pattern = ".gri$", full.names = TRUE)
hymap_stack <- stack(hymap_files)
h_names <- names(hymap_stack)
doe_write_raster(hymap_stack, "d:/mineral_markers_manuscript/fusion/hymap_fusion_all")

min_hymap <-raster::minValue(hymap_stack)
max_hymap <-raster::maxValue(hymap_stack)
hymap_stack <- (hymap_stack-min_hymap)/(max_hymap - min_hymap)
min_hymap <- max(raster::minValue(hymap_stack))
hymap_stack[hymap_stack<=min_hymap]<-NA
names(hymap_stack) <- h_names
doe_write_raster(hymap_stack, "d:/mineral_markers_manuscript/fusion/hymap_fusion_all_normal")

h <- hymap_stack[[c("mtmf_fusion", "mttcimf_fusion", "osp_fusion", "sam_fusion")]]
h <- raster::calc(h, fun = sum, na.rm=TRUE)
h_max <- max(raster::quantile(h, na.rm=T, probs=c(0.99))) # Winsorizing to percentile 99
h <- h/h_max
h[h>1]<-1
names(h)<-"hymap_fusion"
doe_write_raster(h, "d:/mineral_markers_manuscript/fusion/hymap_fusion_final")
spplot(h)
h


################### Gather all fused results
#### ####

fuse_layers <- function(input_directory, output_directory, file_basename, 
                        file_extension = ".gri", 
                        layers_to_fuse = c("mtmf_fusion", "mttcimf_fusion", 
                                           "osp_fusion", "sam_fusion"), 
                        windsor_percentile = 0.99)
{
  input_files <- list.files(path = input_directory, 
                            pattern = paste0(file_extension,"$"), 
                            full.names = TRUE)
  input_stack <- stack(input_files)
  input_stack_names <- names(input_stack)
  cat("Input file/layer names:", fill = TRUE)
  cat(input_stack_names, sep = ", ", fill = TRUE)
  if(!prod(unlist(lapply(layers_to_fuse, 
                         function(x) x %in% input_stack_names)))){
    cat("Members of subset layers don't exist: ", layers_to_fuse, fill=TRUE)
    return(-1)
  }

  files_basename <- file.path(output_directory, file_basename)
  stack_filename <- paste0(files_basename, "_all")
  normal_filename <- paste0(files_basename, "_all_normal")
  final_filename <- paste0(files_basename, "_final")
  doe_write_raster(input_stack, stack_filename) # Save as single file
  input_quant <- c(quantile(input_stack, prob=c(windsor_percentile)))
  input_stack <- (input_stack)/(input_quant) # Winsorizing on the positive side only
  input_stack[input_stack>1] <- 1
  names(input_stack) <- input_stack_names
  doe_write_raster(input_stack, normal_filename)
  b <- input_stack[[layers_to_fuse]]
  b <- raster::calc(b, fun = sum, na.rm=TRUE)
  b[b<=0] <- NA # Eliminate zeros (result of adding NA's)
  b_windsor <- max(raster::quantile(b, na.rm=T, probs=c(windsor_percentile)))
  b <- b/b_windsor
  b[b<=0] <- NA # Eliminate zeros (result of adding NA's)
  b[b>1] <- 1
  names(b)<-paste0(file_basename, "_fusion")
  doe_write_raster(b, final_filename)
  return(b)
}

b <- fuse_layers("d:/mineral_markers_manuscript/fusion/brady", 
                 "d:/mineral_markers_manuscript/fusion", "brady_fusion")
spplot(b)
d <- fuse_layers("d:/mineral_markers_manuscript/fusion/desert_peak", 
                 "d:/mineral_markers_manuscript/fusion", "desert_fusion")
spplot(d)
h <- fuse_layers("d:/mineral_markers_manuscript/fusion/hymap", 
                 "d:/mineral_markers_manuscript/fusion", "hymap_fusion")
spplot(h)

