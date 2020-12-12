if(!require(pacman)){
  installed.packages("pacman")
  library("pacman")
}
pacman::p_load(raster, rgdal, sp)


#extent_brady <- raster(xmn=327499.1, xmx=329062.1, 
#                       ymn = 4405906, ymx=4409320, res=c(3,3), 
#                       crs="+init=epsg:32611")

extent_brady <- raster(xmn=326440.1, xmx=329062.1, 
                       ymn = 4405094, ymx=4409320, res=c(3,3), 
                       crs="+init=epsg:32611")


extent_desert <- raster(xmn=330955.3, xmx=335989.3, 
                        ymn = 4401681, ymx=4406211, res=c(3,3), 
                        crs="+init=epsg:32611")

extent_hymap <- raster(xmn=325209, xmx=336531, ymn = 4395438, ymx=4412103, 
                       res=c(3,3), crs=crs('+init=epsg:32611'))

doe_writeRaster <- function(x, filename, format="raster", overwrite=TRUE, bandorder="BSQ"){
  if(tools::file_ext("filename.grd") != "grd") {
    filename <- tools::file_path_sans_ext(filename)
    filename <- paste(filename, ".grd", sep="")
  }
  f1<-writeRaster(x=x, filename=filename, bandorder=bandorder,
                  format=format, overwrite=overwrite)
  hdr(f1, "ENVI")
  return(f1)
}


test_all_thresholds <- function(x) {
  x = as.matrix(x)*1000
  mode(x) <- "integer"
  
  the_good_methods <-c("IJDefault", "Intermodes", #"IsoData", 
                       "Minimum", "Moments", "Otsu", "RenyiEntropy")
  r <- vector(mode="list", length = 0)
  
  for (th_method in the_good_methods){ 
    b_t <- auto_thresh(x, method = th_method, 
                       ignore_black = TRUE, ignore_na = TRUE)
    r[th_method] = (as.numeric(b_t)/1000.0)
    print(paste("Method:", th_method, "; threshold = ", 
                as.numeric(b_t)/1000.0) )
  }
  return (r)
}

plot_threshold <- function(x, threshold, e=NULL){
  if(!is.null(e)){
    x <- crop(x, e)
  }
  top_x <- brick(x)
  top_x[top_x<(as.numeric(threshold))] <- NA
  p <- spplot(top_x, main=paste0("Threshold Applied (", threshold, ")", sep=""))
  print(p)
  return(p)
}


minerals_directory = "../../HyMapTargetDetectionFull"
sam_file = paste(minerals_directory, "HyMap_full_sam_rule", sep="/")
sam <- stack(sam_file)
# names(sam) <- c("OpalizedTuff_CU00_15E", "KrattOpal", "Epsomite_GDS149", "Gypsum_SU2202", "Hematite_FE2602", "Kaolinite_CM3", "Chalcedony_CU91_6a")
band_names <- c("OpalizedTuff", "KrattOpal", "Epsomite", "Gypsum", "Hematite", "Kaolinite", "Chalcedony")
names(sam) <- band_names
strict = 0.1
moderate = 0.2
relaxed = 0.3
center_point = 1.0

min_limit = center_point-strict
max_limit = center_point+strict

filtered_sam <- stack(sam_file)

sam_filter <- function(x, center_point, threshold, e=NULL){
  threshold = abs(threshold)
  min_limit = center_point-threshold
  max_limit = center_point+threshold
  if(!is.null(e)){
    x_crop <- crop(x, e)
  } else {
      x_crop <- x
  }
  x_crop[(x_crop<min_limit)|(x_crop>max_limit)] <- NA
  return (x_crop)
  
}

sam_normalize <- function(x){
  x <- abs(x-1)
  x <- 1-x
  return (x)
}

raster_sum <- function(x){
  idx = c(1:nlayers(x))*0+1
  r <- stackApply(x, indices = idx, fun = sum, na.rm = TRUE)
  return(r)
}

stack2matrix <- function(x, layer_name){
  b = as.matrix(x[[layer_name]])
  mode(b) <- "integer"
  return(b)
}


low_normalize <- function(x){
  x[x==0]<-NA
  x <- abs(x)
  x <- 1-x
  return (x)
}


# brady_sam <- crop(filtered_sam, extent_brady)
brady_sam <- filtered_sam
brady_sam <- sam_normalize(brady_sam)
brady_relaxed <- sam_filter(brady_sam, center_point, relaxed)  ### , extent_brady)
brady_moderate <- sam_filter(brady_sam, center_point, moderate)  ### , extent_brady)
brady_strict <- sam_filter(brady_sam, center_point, strict)  ### , extent_brady)
# spplot(brady_sam, main = "SM Brady (strict)")
names(brady_relaxed) <- band_names
names(brady_moderate) <- band_names
names(brady_strict) <- band_names
spplot(brady_relaxed, main = "SAM Brady (relaxed)")
spplot(brady_moderate, main = "SAM Brady (moderate)")
spplot(brady_strict, main = "SAM Brady (strict)")

f1 <- doe_writeRaster(brady_relaxed, paste(minerals_directory, "SAM_relaxed", sep="/"))
rm(f1)
dir.create("results/SAM", recursive = TRUE, mode = "0775")

brady_minerals_s <- raster_sum(brady_strict)
brady_minerals_m <- raster_sum(brady_moderate)
brady_minerals_r <- raster_sum(brady_relaxed)
# brady_minerals_sum <- stackApply(brady_minerals, indices = (c(1:nlayers(brady_sam))*0+1),fun = sum, na.rm = TRUE)
spplot(brady_minerals_s, main = "Brady SAM Summary (strict)")
spplot(brady_minerals_m, main = "Brady SAM Summary (moderate)")
spplot(brady_minerals_r, main = "Brady SAM Summary (relaxed)")
# spplot(brady_minerals_sum, main = "Brady SAM Summary (relaxed)")

limit = cellStats(brady_minerals_r, max)
brady_minerals_sum = brady_minerals_r / limit
spplot(brady_minerals_sum, main = "Brady SAM Relaxed (normalized)")

f1 <- doe_writeRaster(brady_minerals_sum, "results/SAM/SAM_Mineral_relaxed")
rm(f1)

brady_sam <- crop(brady_minerals_sum, extent_brady)
desert_sam <- crop(brady_minerals_sum, extent_desert)

top_sam <- brick(brady_sam)
top_sam[top_sam<0.4] <- NA
spplot(top_sam, main="SAM Brady Only")


pacman::p_load(autothresholdr)

b <- as.matrix(brady_sam)*1000
b = as.matrix(brady_minerals_sum)*1000
mode(b) <- "integer"

# The available methods are "IJDefault", "Huang", "Huang2", "Intermodes", 
# "IsoData", "Li", "MaxEntropy", "Mean", "MinErrorI", "Minimum", "Moments",
# "Otsu", "Percentile","RenyiEntropy", "Shanbhag", "Triangle" and "Yen". 
th_methods <-c("IJDefault", "Huang", "Huang2", "Intermodes",
               "IsoData", "Li", "MaxEntropy", "Mean", "MinErrorI", 
               "Minimum", "Moments", "Otsu", "Percentile",
               "RenyiEntropy", "Shanbhag", "Triangle", "Yen")
th_good_methods <-c("IJDefault", "Intermodes","IsoData", "Minimum", 
                    "Moments", "Otsu", "RenyiEntropy")

for (th_method in th_good_methods){ 
  b_t <- auto_thresh(b, method = th_method, ignore_black = TRUE, ignore_na = TRUE)
  print(paste("Method:", th_method, "; threshold = ", 
              as.numeric(b_t)/1000.0) )
}
b_t <- auto_thresh(b, method = "Otsu", ignore_black = TRUE, ignore_na = TRUE)
top_sam <- brick(brady_sam)
top_sam[top_sam<(as.numeric(b_t)/1000.0)] <- NA
spplot(top_sam, main=paste0("SAM Brady Threshold (", th_method, ")", sep=""))

rm(brady_relaxed, brady_moderate, brady_strict, brady_minerals)
rm(brady_minerals_m, brady_minerals_r, brady_minerals_s)
rm(sam, filtered_sam, b, brady_minerals_sum)
rm(brady_sam, top_sam, desert_sam)

################## Done with SAM #####################

ace_file = paste(minerals_directory, "HyMap_full_ace", sep="/")
ace_stack <- low_normalize(stack(ace_file))
# ace_stack <- crop(ace_stack, extent_brady)
names(ace_stack) <- band_names
spplot(ace_stack)


# b <- stack2matrix(ace_stack, "Chalcedony")
r<-test_all_thresholds(ace_stack[["Chalcedony"]])
summary(ace_stack[["Chalcedony"]])

plot_threshold(ace_stack[["Chalcedony"]], 0.843, extent_brady)
