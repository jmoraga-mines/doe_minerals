
layer_names <- c("Chalcedony","Kaolinite", "Gypsum", "Hematite", "Epsomite")
layer_names <- c("Calcite", "Chalcedony",  "Chlorite", "Epidote", "Hematite",
                 "Illite", "Kaolinite", "Opal", "Quartz")
usable_layers <- c("Opal", "Kaolinite", "Epidote", "Hematite")
ace_th <- c(0.586, 0.585, 0.583, 0.603)
cem_th <- c(0.744, 0.7, 0.503, 0.527)
mf_th <- c(0.744, 0.7, 0.503, 0.527)
sam_th <- c(0.821, 0.827, 0.814, 0.812)
osp_th <- c(0.662, 0.516, 0.728, 0.468)
mtmf_th <- c(0.041, 0.053, 0.0425, 0.059)
tcimf_th <- c(0.634, 0.544, 0.667, 0.46)
mttcimf_th <- c(0.681, 0.316, 0.291, 0.337)



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
                          usable_layers=NULL,
                          plot_colors=(rev(heat.colors(12)))[6:12],
                          is_sam=FALSE, has_infeasibility=FALSE,
                          cropping_extent=NULL, plot_crop=FALSE, 
                          output_file_name=NULL, plot_fusion=FALSE,
                          normalize_fusion=FALSE, verbose=FALSE)
{
  ### Validate parameters
  # File is valid
  if(!file.exists(file_name)){
    cat(paste("File does not exist:", file_name))
    return(NULL)
  }
  s <- stack(file_name)
  # File data is of type raster or stack
  if(!exists("s")) {
    cat("Error loading file into raster object\n")
    return(NULL)
  }
  # No usable_layers parameter, use all layers as defaults
  if (is.null(usable_layers)){
    usable_layers <- layer_names
  }
  n<-length(layer_names)
  # Validate correct number of layer names and layers in file
  if((has_infeasibility & ((2*n)!=nlayers(s))) # has infeasibility
     | (!has_infeasibility &(n!=nlayers(s)))){   # Number of bands names check
    cat("Length of layer names and bands in file do not match")
    return(NULL)
  }
  # Verify all band names that are usable belong to layer names
  if(sum(usable_layers %in% layer_names) != length(usable_layers)) {
    cat("Usable layers don't match members of known layers\n")
    cat("Offending layers:\n")
    cat(usable_layers[!(usable_layers %in% layer_names)])
    return(NULL)
  }
  # Check each usable layer has a threshold
  if (length(thresholds)==1){
    cat("Using threshold: percentile ", sprintf("%0.2f%%", thresholds*100), "\n")
  } else if(length(usable_layers)!=length(thresholds)){
    cat("Length of usable layer names and thresholds do not match")
    return(NULL)
  }
  if (has_infeasibility){
    cat("Splitting between feasible and infeasible\n")
    s_infeasible <- s[[(n+1):(2*n)]]
    # selects only usable layers from infeasibility layers
    names(s_infeasible)<-layer_names
    s <- s[[1:n]]
    names(s) <- layer_names
    cat("Subsetting usable layers\n")
    s_infeasible <- s_infeasible[[usable_layers]]
    # Crops infeasibility layers if needed
    cat("Cropping infeasibility layers\n")
    if(!is.null(cropping_extent))
      s_infeasible <- crop(s_infeasible, cropping_extent)
  }
  # Uses only usable layers from now on
  names(s) <- layer_names
  s <- s[[usable_layers]]
  layer_names <- usable_layers
  n<-length(layer_names)
  # Crops usable layers if needed
  cat("Cropping layers\n")
  if(!is.null(cropping_extent))
    s_crop <- crop(s, cropping_extent)
  else
    s_crop<-s
  if(has_infeasibility){
    cat("Applying infeasibility curve\n")
    infeasibility_curve <- ((s_crop*300 + 15 - s_infeasible)>0)
    s_crop <- s_crop*infeasibility_curve
    rm(s_infeasible)
  }
  l<-nlayers(s_crop)
  if(is_sam) {
    cat("SAM analysis:\n")
    s_crop[s_crop==0] <- NA
    print(s_crop)
    s_crop <- (pi-abs(s_crop))/pi
    cat("After applying threshold and pi normalization:\n")
    print(s_crop)
  }
  names(s_crop) <- usable_layers
  if ((length(thresholds)==1) & (length(usable_layers)>1)){
    q <- quantile(s_crop, rm.na=T, probs=c(thresholds))
    cat("Calculated thresholds:\n")
    print(q)
    thresholds <- as.numeric(q)
  }

  s_crop[s_crop<thresholds] <- NA
  names(s_crop) <- usable_layers
  if(plot_crop) plot(spplot(s_crop, col=plot_colors))
  cat("Fuse all bands\n")
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



### Coso
############ Coso maps
fusion_ace_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_ace", 
                                 "ace_fusion", ace_th, 
                                 layer_names =layer_names,
                                 usable_layers = usable_layers,
                                 output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_ace")
fusion_cem_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_cem", 
                                  "cem_fusion", cem_th,
                                  layer_names =layer_names,
                                  usable_layers = usable_layers,
                                  output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_cem")
fusion_mf_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_mf", 
                                 "mf_fusion", mf_th,
                                 layer_names =layer_names,
                                 usable_layers = usable_layers,
                                 output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_mf")
fusion_mtmf_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_mtmf", 
                                   "mtmf_fusion", mtmf_th, 
                                   layer_names =layer_names,
                                   usable_layers = usable_layers,
                                   has_infeasibility = T, 
                                   output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_mtmf")
fusion_osp_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_osp", 
                                 "osp_fusion", osp_th,
                                 layer_names =layer_names,
                                 usable_layers = usable_layers,
                                  output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_osp")
fusion_sam_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_sam_rule", 
                                  "sam_fusion", sam_th,
                                  layer_names =layer_names,
                                  usable_layers = usable_layers, is_sam=TRUE,   
                                  plot_crop = TRUE, plot_fusion = TRUE, 
                                  output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_sam")
fusion_tcimf_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_tcimf", 
                                    "tcimf_fusion", tcimf_th,
                                    layer_names =layer_names,
                                    usable_layers = usable_layers, 
                                    output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_tcimf")
fusion_mttcimf_coso <- create_fusion("e:/Coso/CosoMinerals/Trial_07_mttcimf", 
                                      "mttcimf_fusion", mttcimf_th,
                                      layer_names =layer_names,
                                      usable_layers = usable_layers,
                                      has_infeasibility = T,  
                                      output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_mttcimf")
###############   Coso maps, using 99% threshold
fusion_ace_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_ace", 
                                 "ace_fusion", 0.99, 
                                 layer_names =layer_names,
                                 usable_layers = usable_layers,
                                 output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_ace_99")
fusion_cem_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_cem", 
                                 "cem_fusion", 0.99,
                                 layer_names =layer_names,
                                 usable_layers = usable_layers,
                                 output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_cem_99")
fusion_mf_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_mf", 
                                "mf_fusion", 0.99,
                                layer_names =layer_names,
                                usable_layers = usable_layers,
                                output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_mf_99")
fusion_osp_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_osp", 
                                 "osp_fusion", 0.99,
                                 layer_names =layer_names,
                                 usable_layers = usable_layers,
                                 output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_osp_99")
fusion_sam_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_sam_rule", 
                                 "sam_fusion", 0.99,
                                 layer_names =layer_names,
                                 usable_layers = usable_layers, is_sam=TRUE,   
                                 plot_crop = TRUE, plot_fusion = TRUE, 
                                 output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_sam_99")
fusion_tcimf_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_tcimf", 
                                   "tcimf_fusion", 0.99,
                                   layer_names =layer_names,
                                   usable_layers = usable_layers, 
                                   output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_tcimf_99")
fusion_mtmf_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_mtmf", 
                                     "mtmf_fusion", 0.99, 
                                     layer_names =layer_names,
                                     usable_layers = usable_layers,
                                     has_infeasibility = T, 
                                     output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_mtmf_99")
fusion_mttcimf_coso_99 <- create_fusion("e:/Coso/CosoMinerals/Trial_07_mttcimf", 
                                     "mttcimf_fusion", 0.99,
                                     layer_names =layer_names,
                                     usable_layers = usable_layers,
                                     has_infeasibility = T,  
                                     output_file_name="e:/Coso/CosoFusion/Trial_07_fusion_mttcimf_99")
#
coso_files <- list.files(path = "e:/Coso/CosoFusion", 
                          pattern = "99.gri$", full.names = TRUE)
coso_stack <- stack(coso_files)
coso_names <- names(coso_stack)
doe_write_raster(coso_stack, "e:/Coso/coso_fusion_all")
min_coso <-raster::minValue(coso_stack)
max_coso <-raster::maxValue(coso_stack)
coso_stack <- (coso_stack-min_coso)/(max_coso - min_coso)
min_coso <- max(raster::minValue(coso_stack))
coso_stack[coso_stack<=min_coso]<-NA
names(coso_stack) <- coso_names
doe_write_raster(coso_stack, "e:/Coso/coso_fusion_all_normal")
coso_stack  <- stack("e:/Coso/coso_fusion_all_normal")

b <- coso_stack[[c("mf_fusion", "tcimf_fusion", "osp_fusion", "sam_fusion")]]
b <- raster::calc(b, fun = sum, na.rm=TRUE)
b[b<=min_coso] <- NA
b_windsor <- raster::quantile(b, na.rm=T, probs=c(0.005, 0.995))
b_max <- max(b_windsor) # Winsorizing to percentile 99
b_min <- min(b_windsor) # Winsorizing to percentile 99
b[b>=b_max]<-b_max
b[b<=b_min]<-b_min

b[b>1]<-1
names(b)<-"coso_fusion"
doe_write_raster(b, "e:/Coso/coso_fusion_final_nomt")
spplot(b, main="Fusion of Fusions")
