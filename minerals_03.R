source("utils/doe_mineral_utils.R")

###############################################################################
## Instructions and setup
###############################################################################
####
####
#### Target detection methods
detection_methods <- c("ace", "cem", "mf", "osp","sam_rule","tcimf")
####
###############################################################################
#### Start - User defined Parameters  #########################################
#### Replace by your system's relative directory containing the
#### results of ENVI's target detection
minerals_directory <- "../../doe_test/HyMapTargetDetectionFull"
#### Define crop area or NULL
crop_area <- NULL
# crop_area <- extent_tall_brady


####
#### Replace by the base name you used when creating the target
#### detection files
base_name <- "HyMap_full"
####
#### Replace by the names of target detection spectra used
band_names <- c("OpalizedTuff", "KrattOpal", "Epsomite", "Gypsum", 
                "Hematite", "Kaolinite", "Chalcedony")
####
#### Replace by the spectra that will be used in this analysis
target_list <- c("Chalcedony", "Kaolinite", "Gypsum", 
                 "Hematite", "Epsomite")
####
#### Target detection methods to be used in the analysis
selected_methods <- c("ace", "mf", "osp","sam_rule","tcimf")
####
#### End   - User defined Parameters  #########################################

###############################################################################
#### Base name to use 
full_base_name <- paste0(minerals_directory, "/", base_name, "_")

### Example of filenames for HyMap project
# input_files <- c("HyMap_full_ace", "HyMap_full_mf", "HyMap_full_osp",
#                 "HyMap_full_sam_rule","HyMap_full_tcimf")
####
####
###############################################################################

raw_TargetDetection <- function(base_name_full, selected_methods, band_names, 
                                 target_list = NULL, crop_area = NULL)
{
  input_files <- unlist(lapply(selected_methods,
                               function(x){paste0(base_name_full,  x)}))
  lapply(input_files, function(x) if(!file.exists(x)){stop(paste0(x, " does not exist\n"), call. = F)})
  if (is.null(target_list)) 
    target_list <- band_names
  all_stacks <- vector(mode="list", length = 0)
  for (file_name in input_files){
    s <- stack(file_name)
    
    if(!is.null(crop_area))
      s <- raster::crop(s, crop_area)
    names(s) <- band_names
    s <- s[[target_list]]
    s[s==0] <- NA
    s <- raster::setMinMax(s)
    # After reading from files, no min-max are defined
    
    #### if it's the output of the SAM algorithm, adjust results
    if(grepl("_sam", file_name, fixed = TRUE)){
      cat(paste0(file_name, " is SAM\n"))
      ### Apply SAM regularization
      s <- abs(s - pi/2)/pi
    } else {
      cat(paste0(file_name, " is not SAM\n"))
      s <- abs(s)
    }
    s <- minmax_normalization(s)
    names(s) <- target_list
    cat("Appending results...\n")
    print(s)
    all_stacks <- append(all_stacks, s)
    rm(s)
  }
  names(all_stacks) <- selected_methods
  return(all_stacks)
}

###############################################################################
### Step 01 - Load all files and normalize                                   ##
###############################################################################
### Read all files and apply min-max normalization

#' Read the results of ENVI's target detection
#'
#' @param base_name_full String. Full base name used for target detection
#' @param selected_methods String. List of target detection methods to load
#' @param band_names String. Names of all bands in the detection file
#' @param target_list String. List of bands to load from each file
#' @param crop_area raster. NULL or raster with extent to crop each image
#'
#' @return raster
#' @export
#'
#' @examples
read_TargetDetection <- function(base_name_full, selected_methods, band_names, 
                                 target_list = NULL, crop_area = NULL)
{
  input_files <- unlist(lapply(selected_methods,
                               function(x){paste0(base_name_full,  x)}))
  lapply(input_files, function(x) if(!file.exists(x)){stop(paste0(x, " does not exist\n"), call. = F)})
  if (is.null(target_list)) 
    target_list <- band_names
  all_stacks <- vector(mode="list", length = 0)
  for (file_name in input_files){
    s <- stack(file_name)
    
    if(!is.null(crop_area))
      s <- raster::crop(s, crop_area)
    names(s) <- band_names
    s <- s[[target_list]]
    s[s==0] <- NA
    # After reading from files, no min-max are defined
    cat("Applying min-max normalization... \n")
    
    #### if it's the output of the SAM algorithm, adjust results
    if(grepl("_sam", file_name, fixed = TRUE)){
      cat(paste0(file_name, " is SAM\n"))
      ### Apply SAM regularization and min-max normalization
      s <- sam_normalize(s)
    } else {
      cat(paste0(file_name, " is not SAM\n"))
      ### Apply min-max normalization
      s <- nonsam_normalize(s)
    }
    s[s<=0] <- NA
    names(s) <- target_list
    cat("Appending results...\n")
    print(s)
    all_stacks <- append(all_stacks, s)
    rm(s)
  }
  names(all_stacks) <- selected_methods
  return(all_stacks)
}

all_stacks <- read_TargetDetection(full_base_name, selected_methods, 
                                   band_names, target_list, crop_area)

###############################################################################
### Step 2 - 
###############################################################################
#### START  ##### Create Baseline for Analysis and Validation #################
# Create results directory if it does not exist
if (!dir.exists("results/baseline"))
  dir.create("results/baseline", showWarnings = FALSE, 
             recursive = TRUE, mode = "0775")

top_percentage <- 5
minerals_vector <- vector(mode = "list", length = 0)
for (mineral_name in target_list){
  m <- stack()
  cat("Analizing ", mineral_name, " layers\n", sep="")
  for (i in 1:length(selected_methods)){
    m_i <- all_stacks[[i]][[mineral_name]]
    m_i_size <- as.integer((m_i@ncols*m_i@nrows)*
                             (top_percentage/100)) # Top items
    top_pixels <- sort(as.vector(m_i), 
                       decreasing = TRUE, na.last = TRUE)[m_i_size]
    
    cat(paste0("Top ", top_percentage, "% > ", 
               format(top_pixels, digits = 1, nsmall = 4), " ", 
               selected_methods[i], " \n"))
    m_i[m_i<top_pixels] <- NA
    detection_method <- selected_methods[i]
    names(m_i) <- paste0(mineral_name,"_", detection_method)
    m <- stack(m, m_i)
  }
  cat("Max-min unit normalization\n")
  m_min   <- raster::minValue(m)
  m_max   <- raster::maxValue(m)
  m_range <- m_max-m_min
  m <- (m-m_min)/m_range
  # names(m) <- mineral_name
  minerals_vector <- append(minerals_vector, m)
  cat("Saving: ", mineral_name, "... Full ... ")
  f1 <- doe_write_raster(m, paste0("results/baseline/All_normalized_", mineral_name))
#  cat("Brady ... ")
#  b <- crop(m, extent_tall_brady)
#  f1 <- doe_write_raster(b, paste0("results/intersection/Brady_", mineral_name))
#  cat("Desert Peak ...")
#  d <- crop(m, extent_desert)
#  f1 <- doe_write_raster(b, paste0("results/intersection/Desert_", mineral_name))
#  rm(b, d, f1)
  cat("\n")
  rm(m)
}

all_data <- stack(minerals_vector)
f1 <- doe_write_raster(all_data, 
                       paste0("results/baseline/All_normalized_rasters"))
rm(f1)
# all_data <- stack("results/baseline/All_normalized_rasters")

minerals_baseline <- stack() # vector( mode="list", length = 0 )
##### Baseline estimates
for(mineral_brick in minerals_vector){
  i    <- calc(mineral_brick, fun = sum, na.rm = TRUE)
  minerals_baseline<-stack(minerals_baseline, i)
  i_th <- (i@ncols*i@nrows)*(top_percentage/100) # Threshold
  itop <- sort(as.vector(i), decreasing = TRUE, na.last = TRUE)[i_th]
  i[i<itop] <- NA
  cat(">>> ", format(i_th, big.mark = ","), "items @threshold:", itop, "\n")
  for(l in seq(nlayers(mineral_brick))){
    layer_name <- names(mineral_brick[[l]])
    l <- mineral_brick[[l]]*i
    cat(format(100*sum(as.vector(!is.na(l)))/i_th, digits = 2, nsmall = 2),
        "% of pixels match from", layer_name, "\n")
    rm(l)
  }
  names(minerals_baseline) <- target_list
  rm(i)
}
minerals_baseline <- unit_normalization(minerals_baseline)
names(minerals_baseline) <- target_list
f1 <- doe_write_raster(minerals_baseline, 
                       paste0("results/baseline/All_minerals_baseline"))
rm(f1)
# rm(all_stacks)

# chalcedony <- stack()
# for (i in 1:length(all_stacks))
#   chalcedony <- stack(chalcedony, all_stacks[[i]][["Chalcedony"]])
# chalcedony <- calc(chalcedony, sum, na.rm = TRUE)
# chalcedony[is.na(all_stacks[[5]][["Chalcedony"]])] <- NA
# plot(chalcedony)

#all_stacks <- lapply(input_files, stack) # Load all files from scratch
#all_stacks <- lapply(all_stacks, raster::setMinMax)
#for(i in seq(length(selected_methods))) {
#  a <- all_stacks[[i]]
#  names(a) <- band_names
#  all_stacks[[i]] <- a[[target_list]]
#}
#### END    ##### Create Baseline for Analysis and Validation #################

cat("Collating minerals\n")
minerals_stack <- stack()
for (mineral_name in target_list){
  m <- stack()
  for (i in 1:length(all_stacks))
    m <- stack(m, all_stacks[[i]][[mineral_name]])
  m <- calc(m, sum, na.rm = TRUE)
  names(m) <- mineral_name
  minerals_stack <- stack(minerals_stack, m)
}
cat("Plotting minerals\n")

plot(minerals_stack, main="Mineral Maps before normalization")

hydrothermal <- calc(minerals_stack, sum, na.rm = TRUE)
names(hydrothermal) <- "Hydrothermal"

plot(hydrothermal, main="Fusion of all available minerals")

target_list2 <- c("Chalcedony", "Kaolinite", "Epsomite")
minerals_stack2 <- minerals_stack[[target_list2]]
hydrothermal2 <- calc(minerals_stack2, sum, na.rm = TRUE)
names(hydrothermal2) <- "Geothermal"
plot(crop(hydrothermal2, extent_tall_brady), main="Fusion of top hydrothermal minerals (Brady)")
plot(crop(hydrothermal2, extent_desert), main="Fusion of top hydrothermal minerals (Desert Peak)")

cat("Finding thresholds\n")
geo_m <- as.matrix(hydrothermal2)*1000
mode(geo_m) <- "integer"
geo<-autothresholdr::auto_thresh(geo_m, method = "Otsu", ignore_black = TRUE, ignore_na = TRUE)/1000.0
print(as.numeric(geo))


geo2 <- sigmoid_normalizer(hydrothermal2)
plot(geo2)

geo2_m <- as.matrix(geo2)*1000
mode(geo2_m) <- "integer"
geo2_th<-autothresholdr::auto_thresh(geo2_m, method = "Otsu",
                                 ignore_black = TRUE, ignore_na = TRUE)/1000.0
cat(paste0("New threshold: ", geo2_th, "\n"))
geo2[geo2<as.numeric(geo2_th)] <- NA
geo2 <- unit_normalization(geo2)
plot(geo2, main="Geothermal alterations after normalization and thresholding", col=rev(heat.colors(20)))

# f1 <- doe_write_raster(geo2, "results/filtered/BradyDesert_Markers")
f1 <- doe_write_raster(minerals_stack, "results/filtered/Minerals_Stack")
rm(f1)
####################### END ######################

##### Load and run basic hydrothermal thresholds ####
##### Using top 2.5% and Otsu for thresholding
source("utils/doe_mineral_utils.R")
minerals_baseline <- stack("results/baseline/All_minerals_baseline")
hydrothermal_baseline <- minerals_baseline
brady_hydro_baseline <- minerals_baseline
hydrothermal_baseline <- minerals_baseline[[c("Chalcedony", "Kaolinite")]] # Only keeps Chalcedony and Kaolinite
brady_hydro_baseline <- crop(hydrothermal_baseline, extent_tall_brady)
brady_hydro_baseline <- calc(brady_hydro_baseline, fun=sum)
plot(brady_hydro_baseline, main="Collapsed mineral detection")
total_points <- brady_hydro_baseline@nrows*brady_hydro_baseline@ncols
top_points <- (2.5/100)*total_points    # top 2.5% of the points
threshold_01 <- sort(brady_hydro_baseline[], decreasing = TRUE)[top_points]
brady_th01 <- brady_hydro_baseline
brady_th01[brady_hydro_baseline<threshold_01] <- 0
plot(brady_th01, main="Collapsed mineral detection (th01)")
brady_th01_despeckle <- focal(brady_th01, w=matrix(c(1,1,1,1,1,1,1,1,1), nrow=3), median, na.rm=T)
brady_th01_despeckle[brady_th01_despeckle<=0] <- NA
brady_th01_despeckle <- unit_normalization(brady_th01_despeckle)
plot(brady_th01_despeckle, main="Collapsed mineral detection (th01 despeckled 3x3 Queen)")
brady_th01_count <- sum((brady_th01_despeckle>=0)[], na.rm=TRUE)
threshold_02 <- as.numeric(autothresholdr::auto_thresh(as.integer(as.matrix(brady_hydro_baseline)*1000),
                                                       method = 'Otsu')/1000.0)
brady_th02 <- brady_hydro_baseline
brady_th02[brady_hydro_baseline<threshold_02] <- 0
plot(brady_th02, main="Collapsed mineral detection (th02:Otsu)")
brady_th02_despeckle <- focal(brady_th02, w=matrix(c(1,1,1,1,1,1,1,1,1), nrow=3), median, na.rm=T)
brady_th02_despeckle[brady_th02_despeckle<=0] <- NA
brady_th02_despeckle <- unit_normalization(brady_th02_despeckle)
plot(brady_th02_despeckle, main="Collapsed mineral detection (th02 despeckled 3x3 Queen)")
brady_th02_count <- sum((brady_th02_despeckle>=0)[], na.rm=TRUE)

bth01 <- brady_th01_despeckle*(brady_hydro_baseline>0)
bth02 <- brady_th02_despeckle*(brady_hydro_baseline>0)
plot(crop(bth01, extent_tall_brady), main="Clean mineral detection (th01 - brady only)")
plot(crop(bth02, extent_tall_brady), main="Clean mineral detection (th02 - brady only)")
plot(crop(bth01, extent_desert), main="Clean mineral detection (th01 - Desert Peak only)")
plot(crop(bth02, extent_desert), main="Clean mineral detection (th02 - Desert Peak only)")



### Don't run!!!
# ace_stack <- stack("results/filtered/ace_stack")
# cem_stack <- stack("results/filtered/cem_stack")
# mf_stack <- stack("results/filtered/mf_stack")
# osp_stack <- stack("results/filtered/osp_stack")
# sam_stack <- stack("results/filtered/sam_stack")
# tcimf_stack <- stack("results/filtered/tcimf_stack")
# 
# chalcedony <- stack(sam_stack[["Chalcedony"]], ace_stack[["Chalcedony"]], cem_stack[["Chalcedony"]], mf_stack[["Chalcedony"]], osp_stack[["Chalcedony"]], tcimf_stack[["Chalcedony"]])
# chalcedony <- calc(chalcedony, sum, na.rm = TRUE)
# 
# f1 <- doe_write_raster(chalcedony, "results/filtered/Brady_Chalcedony")

all_stacks2 <- raw_TargetDetection(full_base_name, selected_methods, 
                                   band_names, target_list, crop_area)
minerals_vector <- vector(mode = "list", length = 0)
for (mineral_name in target_list){
  m <- stack()
  for (i in 1:length(all_stacks2)){
    m_i <- all_stacks2[[i]][[mineral_name]]
    detection_method <- selected_methods[i]
    detection_method <- detection_method[[min(3, length(detection_method))]]
    names(m_i) <- paste0(mineral_name,"_", detection_method)
    m <- stack(m, m_i)
  }
  # names(m) <- mineral_name
  minerals_vector <- append(minerals_vector, m)
  cat("Saving: ", mineral_name, "... Full ... ")
  f1 <- doe_write_raster(m, paste0("results/filtered/All_", mineral_name))
  cat("Brady ... ")
  b <- crop(m, extent_tall_brady)
  f1 <- doe_write_raster(b, paste0("results/filtered/Brady_", mineral_name))
  cat("Desert Peak ...\n")
  d <- crop(m, extent_desert)
  f1 <- doe_write_raster(b, paste0("results/filtered/Desert_", mineral_name))
  rm(m, b, d, f1)
}

