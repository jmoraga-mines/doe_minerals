if (!require("pacman")) {
  installed.packages("pacman")
  library("pacman")
}
pacman::p_load(raster, rgdal, sp, rgeos, sf)
pacman::p_load(autothresholdr, igraph)
pacman::p_load(latticeExtra)

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


test_all_thresholds <- function(x, ignore_black=TRUE, ignore_white=FALSE) {
  x <- as.integer(as.matrix(x) * 1000)
  the_good_methods <- c(
    "IJDefault",
    "Intermodes",
    # "IsoData",
    "Minimum",
    "Moments",
    "Otsu",
    "RenyiEntropy"
  )
  all_methods <- c(
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
  r <- vector(mode = "list", length = 0)
  for (th_method in all_methods) {
    b_t <- tryCatch(
      {
        autothresholdr::auto_thresh(x, method = th_method,
                                    ignore_black = ignore_black,
                                    ignore_white = ignore_white,
                                    ignore_na = TRUE)
      }, error=function(cond){cat("\n"); return(NA)}, 
      warning=function(cond){cat("\n"); return(NA)}
    )
    r[th_method] <- (as.numeric(b_t) / 1000.0)
    cat(paste(
      "Method:",
      th_method,
      "; threshold = ",
      as.numeric(b_t) / 1000.0,
      "\n"
    ))
  }
  return(unlist(r))
}


select_largest_threshold <- function(x, ignore_black=TRUE, test_all=FALSE, ignore_white=FALSE) {
  x <- round(as.matrix(x) * 1000, digits = 0)
  the_good_methods <- c(
    "IJDefault",
    "Intermodes",
    "IsoData",
    "Moments",
    "RenyiEntropy",
    "Otsu"
  )
  all_methods <- c(
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
  r <- list("method"="None", "th"=0)
  if (test_all)  {
    test_methods = all_methods
  } else
    test_methods = the_good_methods
  for (th_method in test_methods) {
    b_t <- tryCatch(
      {
        autothresholdr::auto_thresh(x, method = th_method,
                                    ignore_black = ignore_black,
                                    ignore_white = ignore_white,
                                    ignore_na = TRUE)
      }, error=function(cond){cat("\n"); return(NA)}, 
      warning=function(cond){cat("\n"); return(NA)}
    )
    test_th <- (as.numeric(b_t)/ 1000.0)
    if(test_th > r$th) {
      r$th <- test_th
      r$method <- th_method
    }
  }
  return(r)
}


plot_threshold <- function(x, threshold, e = NULL, print_plot = FALSE) {
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


raster_filter <- function(x, center_point, threshold, e = NULL) {
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
  m <- as.matrix(x[[layer_name]])
  mode(m) <- "integer"
  return(m)
}

#' Unit normalization, it converts the values of a raster to the
#' range [0,1]. Original zero values are turned to NA.
#'
#' @param x Raster. An input raster
#'
#' @return Raster. a normalized raster with values from 0 to 1
#'
#' @examples
#' > library(raster)
#' > r <- raster::raster(nrows = 50, ncols = 50)
#' > set.seed(42)
#' > set.seed(42)
#' > r <- setValues(r,runif(ncell(r)))*5
#' > r
#' class      : RasterLayer 
#' dimensions : 50, 50, 2500  (nrow, ncol, ncell)
#' resolution : 7.2, 3.6  (x, y)
#' extent     : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#' crs        : +proj=longlat +datum=WGS84 +no_defs 
#' source     : memory
#' names      : layer 
#' values     : 0.001194483, 4.998816  (min, max)
#' > unit_normalization(r)
#' class      : RasterLayer 
#' dimensions : 50, 50, 2500  (nrow, ncol, ncell)
#' resolution : 7.2, 3.6  (x, y)
#' extent     : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#' crs        : +proj=longlat +datum=WGS84 +no_defs 
#' source     : memory
#' names      : layer 
#' values     : 0, 1  (min, max)
unit_normalization <- function(x) {
  # x[x==0] <- NA
  the_names <- names(x)
  x_min <- raster::minValue(x) # vector with min for all layers
  x_max <- raster::maxValue(x) # vector with max for all layers
  x_span <- x_max-x_min
  x <- (x-x_min)/x_span
  names(x) <- the_names
  return(x)
}

minmax_normalization <- function(x) {
  # x[x==0] <- NA
  the_names <- names(x)
  # Obtain global minimum and maximum
  x_min <- min(raster::cellStats(x, stat = 'min', rm.na = TRUE, asSample= FALSE))
  x_max <- max(raster::cellStats(x, stat = 'max', rm.na = TRUE, asSample= FALSE))
  x_span <- x_max-x_min
  x <- (x-x_min)/x_span
  names(x) <- the_names
  return(x)
}


sigmoid_normalizer <- function(x, sigmas = 2){
  the_names <- names(x)
  x_mean <- raster::cellStats(x, stat=mean, na.rm=TRUE, asSample=FALSE)
  x_std <-  raster::cellStats(x, stat=sd, na.rm=TRUE, asSample=FALSE)
  x <- 1/(1+exp((x_mean-x)/(sigmas*x_std)))
  names(x) <- the_names
  return(x)
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
sam_normalize <- function(x, zero2na = FALSE) {
  if (zero2na)
    x[x==0] <- NA
  x <- abs(x - pi/2)/pi
  x <- raster::setMinMax(x)
  x <- unit_normalization(x)
  return(x)
}


#' Normalizes a non-SAM raster to a number between 0-1, where 1 is better
#' To compare against other methods, this function converts values to 
#' absolute low numbers and normalizes values to the range 0-1
#' where a 1 is perfect match and 0 is no match.
#' 
#' @param x Raster. A raster with the results of the SAM analysis
#'
#' @return Raster. A normalized raster
#'
#' @examples
#' n_sam <- nonsam_normalize(sam_raster)
#' 
nonsam_normalize <- function(x, zero2na = FALSE) {
  if (zero2na)
    x[x==0] <- NA
  x <- abs(x)
  x <- raster::setMinMax(x)
  x <- unit_normalization(x)
  return(x)
}


extract_best_matches <- function(x, top_percentage = 0.025) {
  top_percentage = abs(top_percentage)
  for(name in names(x)){
    aLayer <- x[[name]]
    q <- quantile(aLayer, 
                  probs = seq((1-top_percentage), 1, top_percentage))
    c <- q[2]
    t <- c-q[1]
    aLayer <- raster_filter(aLayer, center_point = c, threshold = t)
    x[[name]] <- aLayer
  }
  return(x)
}


load_minerals <- function(minerals_file, 
                          band_names = c("OpalizedTuff", "KrattOpal",
                                         "Epsomite", "Gypsum", 
                                         "Hematite", "Kaolinite",
                                         "Chalcedony"),
                          target_list = c("Chalcedony", "Kaolinite",
                                           "Gypsum", "Epsomite"),
                          is_sam = FALSE, e=NULL){
  minerals_stack <- stack(minerals_file)
  if (!is.null(e)) {
    print("Cropping...")
    minerals_stack <- raster::crop(minerals_stack, e)
  }
  names(minerals_stack) <- band_names
  minerals_stack <- minerals_stack[[target_list]]
  raster::setMinMax(minerals_stack)
  if (is_sam)
    minerals_stack <- sam_normalize(minerals_stack)
  else
    minerals_stack <- unit_normalization(minerals_stack)
  return(minerals_stack)
}


full_detection <- function(base_dir, file_name, band_names, 
                           target_list, e = NULL, is_sam = FALSE,
                           ignore_black = FALSE, 
                           top_threshold = 0.05) {
  detection_file <- file.path(base_dir, file_name)
  print(paste0("Reading and normalizing: ", detection_file))
  detection_stack <- load_minerals(detection_file, 
                                   band_names = band_names,
                                   target_list = target_list,
                                   is_sam = is_sam, e= e)
  if (is_sam){
    detection_stack <- sam_normalize(detection_stack)
  } else
    detection_stack <- unit_normalization(detection_stack)
  print(detection_stack)
  print("Checking thresholds for each layer")
  names(detection_stack) <- target_list
  b_t <- 0
  for (n in names(detection_stack)) {
    print(paste0("Layer: ", n))
    aLayer <- detection_stack[[n]]
    print("Select threshold")
    q <- quantile(aLayer, probs=seq((1-top_threshold), 1, top_threshold))
    c <- q[2]
    t <- c-q[1]
    aLayer <- raster_filter(aLayer, center_point = c, threshold = t)
    aLayer <- unit_normalization(aLayer)
    t <- select_largest_threshold(aLayer,
                                  ignore_black = ignore_black)
    print(paste0("Method: ", t$method, "; Threshold: ", t$th))
    aLayer[aLayer < b_t] <- 0
    aLayer <- unit_normalization(aLayer)
    b_t <- t$th
    th_method <- t$method
    top_sam <- detection_stack[[n]]
    top_sam[top_sam < b_t] <- 0
    top_sam <- unit_normalization(top_sam)
    detection_stack[[n]] <- top_sam
  }
  return(detection_stack)
}

