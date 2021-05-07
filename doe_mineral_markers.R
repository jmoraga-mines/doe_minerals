if (!require("pacman")) {
  installed.packages("pacman")
  library("pacman")
}
pacman::p_load(raster, rgdal, sp, rgeos)
pacman::p_load(autothresholdr)
pacman::p_load(latticeExtra)

source("LST/doe_utilities.R")

### START # Defining geographical extent from Areas of Interest ---------------
###

# Original Brady extent: (xmn=327499.1, xmx=329062.1, ymn = 4405906,
#                         ymx=4409320, res=c(3,3), crs="+init=epsg:32611")

extent_tall_brady <- raster::raster(
  xmn = 326440.1,
  xmx = 329062.1,
  ymn = 4405094,
  ymx = 4409320,
  res = c(3, 3),
  crs = "+init=epsg:32611"
)

poly_extent_brady <- as(extent(extent_tall_brady), 'SpatialPolygons')
crs(poly_extent_brady) <- crs(extent_tall_brady)

extent_desert <- raster::raster(
  xmn = 330955.3,
  xmx = 335989.3,
  ymn = 4401681,
  ymx = 4406211,
  res = c(3, 3),
  crs = "+init=epsg:32611"
)

poly_extent_desert <- as(extent(extent_desert), 'SpatialPolygons')
crs(poly_extent_desert) <- crs(extent_desert)

extent_hymap <- raster::raster(
  xmn = 325209,
  xmx = 336531,
  ymn = 4395438,
  ymx = 4412103,
  res = c(3, 3),
  crs = "+init=epsg:32611"
)
### END  # Defining geographical extent from Areas of Interest ----------------



doe_write_raster <- function(x,
                             filename,
                             format = "raster",
                             overwrite = TRUE,
                             bandorder = "BSQ") {
  if (tools::file_ext(filename) != "grd") {
    filename <- tools::file_path_sans_ext(filename)
    filename <- paste0(filename, ".grd")
  }
  f1 <- raster::writeRaster(
    x = x,
    filename = filename,
    bandorder = bandorder,
    format = format,
    overwrite = overwrite
  )
  raster::hdr(f1, "ENVI")
  return(f1)
}


test_all_thresholds <- function(x, ignore_black=TRUE) {
  x <- as.matrix(x) * 1000
  mode(x) <- "integer"
  the_good_methods <- c(
    "IJDefault",
    "Intermodes",
    # "IsoData",
    "Minimum",
    "Moments",
    "Otsu",
    "RenyiEntropy"
  )
  r <- vector(mode = "list", length = 0)
  for (th_method in the_good_methods) {
    b_t <- tryCatch(
      {
        autothresholdr::auto_thresh(b, method = th_method,
                                    ignore_black = ignore_black,
                                    ignore_na = TRUE)
      }, error=function(cond){return(-1)}, 
      warning=function(cond){return(-1)}
    )
    r[th_method] <- (as.numeric(b_t) / 1000.0)
    print(paste(
      "Method:",
      th_method,
      "; threshold = ",
      as.numeric(b_t) / 1000.0
    ))
  }
  return(r)
}


select_largest_threshold <- function(x, ignore_black=TRUE) {
  x <- as.matrix(x) * 1000
  mode(x) <- "integer"
  the_good_methods <- c(
    "IJDefault",
    "Intermodes",
    "IsoData",
    "Minimum",
    "Moments",
    "Otsu",
    "RenyiEntropy"
  )
  r <- list("method"="None", "th"=0)
  for (th_method in the_good_methods) {
    b_t <- tryCatch(
      {
        autothresholdr::auto_thresh(b, method = th_method,
                                    ignore_black = ignore_black,
                                    ignore_na = TRUE)
      }, error=function(cond){return(-1)}, 
      warning=function(cond){return(-1)}
    )
    test_th <- (as.numeric(b_t)/ 1000.0)
    if(test_th > r$th) {
      r$th <- test_th
      r$method <- th_method
    }
  }
  return(r)
}

plot_threshold <- function(x,
                           threshold,
                           e = NULL,
                           print_plot = FALSE) {
  if (!is.null(e)) {
    x <- raster::crop(x, e)
  }
  top_x <- raster::brick(x)
  top_x[top_x < (as.numeric(threshold))] <- NA
  if (print_plot) {
    p <- raster::spplot(top_x,
      main = paste0("Threshold Applied (", threshold, ")")
    )
    print(p)
  }
  return(top_x)
}

minerals_directory <- "../../doe_test/HyMapTargetDetectionFull"
sam_file <- file.path(minerals_directory, "HyMap_full_sam_rule")
sam <- raster::stack(sam_file)
# All thresholding methods: "OpalizedTuff_CU00_15E", "KrattOpal",
#    "Epsomite_GDS149", "Gypsum_SU2202", "Hematite_FE2602", "Kaolinite_CM3",
#    "Chalcedony_CU91_6a"
band_names <- c(
  "OpalizedTuff",
  "KrattOpal",
  "Epsomite",
  "Gypsum",
  "Hematite",
  "Kaolinite",
  "Chalcedony"
)
names(sam) <- band_names
strict <- 0.1
moderate <- 0.2
relaxed <- 0.3
center_point <- 1.0

min_limit <- center_point - strict
max_limit <- center_point + strict

brady_sam <- raster::stack(sam_file)

sam_filter <- function(x, center_point, threshold, e = NULL) {
  threshold <- abs(threshold)
  min_limit <- center_point - threshold
  max_limit <- center_point + threshold
  if (!is.null(e)) {
    x_crop <- raster::crop(x, e)
  } else {
    x_crop <- x
  }
  x_crop[(x_crop < min_limit) | (x_crop > max_limit)] <- NA
  return(x_crop)
}

#' Normalizes a SAM raster to a number between 0-1, where 1 is better
#' The results of the Spectral Angle Mapper (SAM) correspond to
#' angles in radians, where a low number indicates closeness between
#' the spectrum tested and the pixel analyzed.
#' To compare against other methods, this function converts these 
#' low numbers (positive or negative) into values in the range 0-1
#' where a 1 is perfect match and 0 is no match.
#' 
#' @param x Raster. A raster with the results of the SAM analysis
#'
#' @return Raster. A normalized raster
#'
#' @examples
#' n_sam <- sam_normalize(sam_raster)
#' 
sam_normalize <- function(x) {
  x <- abs(x - 1)
  x <- 1 - x
  return(x)
}


raster_sum <- function(x) {
  idx <- c(1:nlayers(x)) * 0 + 1
  r <- raster::stackApply(x, indices = idx, fun = sum,na.rm = TRUE)
  return(r)
}

raster_euclidean_collapse <- function(x) {
  x <- x*x
  idx <- rep(1, nlayers(x))
  r <- raster::stackApply(x, indices = idx, fun = sum,na.rm = TRUE)
  r <- sqrt(r)
  return(r)
}


stack2matrix <- function(x, layer_name) {
  b <- as.matrix(x[[layer_name]])
  mode(b) <- "integer"
  return(b)
}


low_normalize <- function(x) {
  # x[x == 0] <- NA
  x <- abs(x)
  return(x)
}

names(brady_sam) <- band_names
brady_sam <- sam_normalize(brady_sam)
brady_relaxed <- sam_filter(brady_sam, center_point, relaxed)
brady_moderate <- sam_filter(brady_sam, center_point, moderate)
brady_strict <- sam_filter(brady_sam, center_point, strict)

names(brady_relaxed) <- band_names
names(brady_moderate) <- band_names
names(brady_strict) <- band_names
spplot(brady_relaxed, main = "SAM Brady (relaxed)")
spplot(brady_moderate, main = "SAM Brady (moderate)")
spplot(brady_strict, main = "SAM Brady (strict)")

f1 <- doe_write_raster(
  brady_relaxed,
  file.path(minerals_directory, "SAM_relaxed")
)
f1 <- doe_write_raster(
  brady_moderate,
  file.path(minerals_directory, "SAM_moderate")
)
f1 <- doe_write_raster(
  brady_strict,
  file.path(minerals_directory, "SAM_strict")
)
rm(f1)

# Create results directory if it does not exist
if (!dir.exists("results/SAM"))
  dir.create("results/SAM", recursive = TRUE, mode = "0775")

# Save shapefiles with boundaries for both sites
raster::shapefile(poly_extent_brady, 
                  filename="results/extent_brady_frame.shp", 
                  overwrite=TRUE)
raster::shapefile(poly_extent_desert, 
                  filename="results/extent_desert_frame.shp", 
                  overwrite=TRUE)

# Collapse multiple layers into one by summing
brady_minerals_s <- raster_sum(brady_strict)
brady_minerals_m <- raster_sum(brady_moderate)
brady_minerals_r <- raster_sum(brady_relaxed)
brady_minerals_s2 <- raster_euclidean_collapse(brady_strict)
brady_minerals_m2 <- raster_euclidean_collapse(brady_moderate)
brady_minerals_r2 <- raster_euclidean_collapse(brady_relaxed)
raster::spplot(brady_minerals_r2, main = "Brady SAM Euclidean Collapse (relaxed)", 
               col.regions=rev(heat.colors(100))) + 
  latticeExtra::layer(sp.polygons(poly_extent_brady, lwd=1, col = "blue")) + 
  latticeExtra::layer(sp.polygons(poly_extent_desert, lwd=1, col = "blue"))


raster::spplot(brady_minerals_s, main = "Brady SAM Summary (strict)")
raster::spplot(brady_minerals_m, main = "Brady SAM Summary (moderate)")
raster::spplot(brady_minerals_r, main = "Brady SAM Summary (relaxed)", 
               col.regions=rev(heat.colors(100))) + 
  latticeExtra::layer(sp.polygons(poly_extent_brady, lwd=1, col = "blue")) + 
  latticeExtra::layer(sp.polygons(poly_extent_desert, lwd=1, col = "blue"))

limit <- raster::cellStats(brady_minerals_r, max)
brady_minerals_sum <- brady_minerals_r / limit
raster::spplot(brady_minerals_sum, 
               main = "Brady SAM Relaxed (normalized)", 
               col.regions=rev(heat.colors(100))) + 
  latticeExtra::layer(sp.polygons(poly_extent_brady, lwd=1, col = "blue")) + 
  latticeExtra::layer(sp.polygons(poly_extent_desert, lwd=1, col = "blue"))

limit <- raster::cellStats(brady_minerals_r2, max)
brady_minerals_euclidean <- brady_minerals_r2 / limit
raster::spplot(brady_minerals_euclidean, 
               main = "Brady SAM Euclidean (relaxed, normalized)", 
               col.regions=rev(heat.colors(100))) + 
  latticeExtra::layer(sp.polygons(poly_extent_brady, lwd=1, col = "blue")) + 
  latticeExtra::layer(sp.polygons(poly_extent_desert, lwd=1, col = "blue"))

f1 <- doe_write_raster(brady_minerals_sum, "results/SAM/SAM_Mineral_relaxed")
f1 <- doe_write_raster(brady_minerals_euclidean, "results/SAM/SAM_Mineral_euclidean_relaxed")
rm(f1)

brady_sam <- raster::crop(brady_minerals_euclidean, extent_tall_brady)
spplot(brady_sam, main="Brady SAM complete (relaxed)", 
       col.regions=rev(heat.colors(100)))
desert_sam <- raster::crop(brady_minerals_euclidean, extent_desert)
spplot(desert_sam, main="Desert Peak SAM complete (relaxed)", 
       col.regions=rev(heat.colors(100)))

top_sam <- raster::brick(brady_sam)
top_sam[top_sam < 0.4] <- NA
plot(top_sam, main = "SAM Brady Only", axes=TRUE)


b <- as.matrix(brady_sam) * 1000
b <- as.matrix(brady_minerals_euclidean) * 1000
mode(b) <- "integer"

# The available methods are "IJDefault", "Huang", "Huang2", "Intermodes",
# "IsoData", "Li", "MaxEntropy", "Mean", "MinErrorI", "Minimum", "Moments",
# "Otsu", "Percentile","RenyiEntropy", "Shanbhag", "Triangle" and "Yen".
th_methods <- c(
  "IJDefault",
  "Huang",
  "Huang2",
  "Intermodes",
  "IsoData",
  "Li",
  "MaxEntropy",
  "Mean",
  "MinErrorI",
  "Minimum",
  "Moments",
  "Otsu",
  "Percentile",
  "RenyiEntropy",
  "Shanbhag",
  "Triangle", 
  "Yen"
)

threshold <- select_largest_threshold(b, ignore_black = TRUE)
b_t <- threshold$th
th_method <- threshold$method
top_sam <- raster::crop(brady_minerals_euclidean, extent_tall_brady)
top_sam[top_sam < b_t] <- 0
raster::spplot(top_sam, main = paste0("SAM Brady Threshold (", th_method, ")"), 
               col.regions=rev(heat.colors(100)))
plot(top_sam, main = paste0("SAM Brady Threshold (", th_method, ")"), 
     axes=TRUE,
     col=rev(heat.colors(100)))
rm(brady_relaxed, brady_moderate, brady_strict, brady_minerals)
rm(brady_minerals_m, brady_minerals_r, brady_minerals_s)
rm(sam, filtered_sam, b, brady_minerals_sum)
rm(brady_sam, top_sam, desert_sam)

################## Done with SAM #####################

ace_file <- file.path(minerals_directory, "HyMap_full_ace")
ace_stack <- low_normalize(stack(ace_file))
names(ace_stack) <- band_names
spplot(ace_stack)


r <- test_all_thresholds(ace_stack[["Chalcedony"]])
summary(ace_stack[["Chalcedony"]])

plot_threshold(ace_stack[["Chalcedony"]], 0.843, extent_tall_brady)

#### Generic Function ####

full_detection <- function(base_dir, file_name, band_names, e = NULL, is_sam = FALSE) {
  detection_file <- file.path(base_dir, file_name)
  print(paste0("Reading and nomalizing: ", detection_file))
  detection_stack <- low_normalize(stack(detection_file))
  if (!is.null(e)) {
    print("Cropping...")
    detection_stack <- raster::crop(detection_stack, e)
  }
  names(detection_stack) <- band_names
  print("Checking thresholds for each layer")
  for (n in band_names) {
    print(paste0("Layer: ", n))
    t <- select_largest_threshold(detection_stack[[n]])
    plot_threshold(detection_stack[[n]], t$th)
  }
}

otsu_detection <- function(base_dir, file_name, band_names, e = NULL) {
  detection_file <- file.path(base_dir, file_name)
  print(paste("Reading:", detection_file))
  detection_stack <- raster::stack(detection_file)
  if (!is.null(e)) {
    print("Cropping...")
    detection_stack <- raster::crop(detection_stack, e)
  }
  print("Normalizing...")
  detection_stack <- low_normalize(detection_stack)
  names(detection_stack) <- band_names
  s <- raster_sum(detection_stack)
  names(s) <- file_name
  t <- test_all_thresholds(s)
  p <-
    plot_threshold(detection_stack[[n]], t$Otsu, print_plot = FALSE)
  raster::spplot(p,
    main = paste("Detection on", file_name, "Threshold:", t$Otsu)
  )
  return(s)
}


b <- quick_detection(minerals_directory, "HyMap_full_cem", band_names)

#################
b_cem <- full_detection(minerals_directory, "HyMap_full_cem", band_names)
doe_write_raster(b_cem, "results/Sum_cem")
b_mf <- full_detection(minerals_directory, "HyMap_full_mf", band_names)
doe_write_raster(b_mf, "results/Sum_mf")
b_osp <- full_detection(minerals_directory, "HyMap_full_osp", band_names)
doe_write_raster(b_osp, "results/Sum_osp")

b_tcimf <- full_detection(minerals_directory, "HyMap_full_tcimf", band_names)
doe_write_raster(b_tcimf, "results/Sum_tcimf")

###################
b_tcimf <- raster::stack(file.path(minerals_directory, "HyMap_full_tcimf"))
b_mttcimf <- raster::stack(file.path(minerals_directory, "HyMap_full_mttcimf"))
b_mtmf <- raster::stack(file.path(minerals_directory, "HyMap_full_mtmf"))
