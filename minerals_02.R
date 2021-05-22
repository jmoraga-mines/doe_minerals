source("utils/doe_mineral_utils.R")

#### Replace by your system's relative directory containing the
#### results of ENVI's target detection
minerals_directory <- "../../doe_test/HyMapTargetDetectionFull"

#### Define crop area or NULL
crop_area <- NULL
crop_area <- extent_tall_brady

#### Replace by the base name you used when creating the target
#### detection files
base_name <- "HyMap_full"

#### Replace by the names of target detection spectra used
band_names <- c("OpalizedTuff", "KrattOpal", "Epsomite", "Gypsum", 
                "Hematite", "Kaolinite", "Chalcedony")

#### Replace by the spectra that will be used in this analysis
target_list <- c("Chalcedony", "Kaolinite", "Gypsum", 
                 "Hematite", "Epsomite")

#### Target detection methods
detection_methods <- c("ace", "cem", "mf", "osp","sam_rule","tcimf")

#### Target detection methods to be used in the analysis
selected_methods <- c("ace", "mf", "osp","sam_rule","tcimf")

#### List with file names to use
input_files <- unlist(lapply(selected_methods, 
                             function(x){paste0(base_name, "_", x)}))

### Example of filenames for HyMap project
# input_files <- c("HyMap_full_ace", "HyMap_full_mf", "HyMap_full_osp",
#                 "HyMap_full_sam_rule","HyMap_full_tcimf")

### Read all files and apply min-max normalization
all_stacks <- lapply(input_files, stack)
for (s in all_stacks){
  if(!is.null(crop_area))
    s <- raster::crop(s, crop_area)
  # s <- raster::crop(s, extent_tall_brady)
  # s <- raster::crop(s, extent_desert)
  names(s) <- band_names
  s <- s[[target_list]]
  #### if it's the output of the SAM algorithm, adjust results
  cat("Setting new min-max... ")
  s <- raster::setMinMax(s)
  cat("Applying min-max normalization... ")
  if(grepl("_sam", file_name, fixed = TRUE)){
    cat(paste0(file_name, " is SAM\n"))
    s <- sam_normalize(s)
  } else {
    cat(paste0(file_name, " is not SAM\n"))
    s[s<=0] <- NA
    ### Apply min-max normalization
    s <- unit_normalization(s)
  }
  names(s) <- target_list
  cat("Appending results...\n")
  print(s)
  all_stacks <- append(all_stacks, s)
  rm(s)
}

# chalcedony <- stack()
# for (i in 1:length(all_stacks))
#   chalcedony <- stack(chalcedony, all_stacks[[i]][["Chalcedony"]])
# chalcedony <- calc(chalcedony, sum, na.rm = TRUE)
# chalcedony[is.na(all_stacks[[5]][["Chalcedony"]])] <- NA
# plot(chalcedony)
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

minerals_stack2 <- minerals_stack[[target_list2]]
hydrothermal2 <- calc(minerals_stack2, sum, na.rm = TRUE)
names(hydrothermal2) <- "Geothermal"
plot(hydrothermal2, main="Fusion of top hydrothermal minerals")

cat("Finding thresholds\n")
geo_m <- as.matrix(hydrothermal2)*1000
mode(geo_m) <- "integer"
geo<-autothresholdr::auto_thresh(geo_m, method = "Otsu",
                                 ignore_black = TRUE, ignore_na = TRUE)
print(as.numeric(geo/1000.0))

sigmoid_normalizer <- function(x){
  x_mean <- raster::cellStats(x, stat=mean)
  x_std <-  raster::cellStats(x, stat=sd)
  x <- 1/(1+exp((x_mean-x)/(2*x_std)))
  return(x)
}
geo2 <- sigmoid_normalizer(hydrothermal2)
plot(geo2)

geo2_m <- as.matrix(geo2)*1000
mode(geo2_m) <- "integer"
geo2_th<-autothresholdr::auto_thresh(geo2_m, method = "Otsu",
                                 ignore_black = TRUE, ignore_na = TRUE)/1000.0
cat(paste0("New threshold: ", geo2_th, "\n"))
geo2[geo2<as.numeric(geo2_th)] <- NA
geo2 <- unit_normalization(geo2)
plot(geo2, main="Geothermal alterations after normalization and thresholding", col=rev(heat.colors(100)))

# f1 <- doe_write_raster(geo2, "results/filtered/BradyDesert_Markers")
f1 <- doe_write_raster(minerals_stack, "results/filtered/Minerals_Stack")

####################### END ######################

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

minerals_vector <- vector(mode = "list", length = 0)
for (mineral_name in target_list){
  m <- stack()
  for (i in 1:length(all_stacks)){
    m_i <- all_stacks[[i]][[mineral_name]]
    detection_method <- strsplit(input_files[i], "_")[[1]]
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



############ Create Baseline for Analysis and Validation ############
# Create results directory if it does not exist
if (!dir.exists("results/intersection"))
  dir.create("results/intersection", showWarnings = FALSE, 
             recursive = TRUE, mode = "0775")

top_percentage <- 5
minerals_vector <- vector(mode = "list", length = 0)
for (mineral_name in target_list){
  m <- stack()
  cat("Analizing ", mineral_name, " layers\n", sep="")
  for (i in 1:length(all_stacks)){
    m_i <- all_stacks[[i]][[mineral_name]]
    m_i_size <- (m_i@ncols*m_i@nrows)*(top_percentage/100) # Top items
    top_pixels <- sort(as.vector(m_i), 
                       decreasing = TRUE, na.last = TRUE)[m_i_size]

    cat(paste0("Top ", top_percentage, "% > ", top_pixels, "\n"))
    m_i[m_i<top_pixels] <- NA
    detection_method <- strsplit(input_files[i], "_")[[1]]
    detection_method <- detection_method[[min(3, length(detection_method))]]
    names(m_i) <- paste0(mineral_name,"_", detection_method)
    m <- stack(m, m_i)
  }
  cat("Max-min unit normalization")
  m_min   <- raster::minValue(m)
  m_max   <- raster::maxValue(m)
  m_range <- m_max-m_min
  m <- (m-m_min)/m_range
  # names(m) <- mineral_name
  minerals_vector <- append(minerals_vector, m)
  cat("Saving: ", mineral_name, "... Full ... ")
  f1 <- doe_write_raster(m, paste0("results/intersection/All_", mineral_name))
  cat("Brady ... ")
  b <- crop(m, extent_tall_brady)
  f1 <- doe_write_raster(b, paste0("results/intersection/Brady_", mineral_name))
  cat("Desert Peak ...\n")
  d <- crop(m, extent_desert)
  f1 <- doe_write_raster(b, paste0("results/intersection/Desert_", mineral_name))
  rm(m, b, d, f1)
}

##### Baseline estimates
for(mineral_brick in minerals_vector){
  i    <- calc(mineral_brick, fun = sum, na.rm = TRUE)
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
  rm(i)
}
