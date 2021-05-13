source("utils/doe_mineral_utils.R")

minerals_directory <- "../../doe_test/HyMapTargetDetectionFull"
band_names <- c("OpalizedTuff", "KrattOpal", "Epsomite", "Gypsum", 
                "Hematite", "Kaolinite", "Chalcedony")
target_list <- c("Chalcedony", "OpalizedTuff", "Kaolinite", "Gypsum", 
                 "Hematite", "Epsomite", "KrattOpal")
target_list2 <- c("Chalcedony", "Kaolinite", "Gypsum", 
                 "Hematite", "Epsomite")

input_files <- c("HyMap_full_ace", "HyMap_full_cem", "HyMap_full_mf",
                 "HyMap_full_osp", "HyMap_full_sam_rule","HyMap_full_tcimf")

### Read all files and apply min-max normalization
all_stacks <- vector(mode="list", length = 0)
for (file_name in input_files){
  s <- stack(file.path(minerals_directory, file_name))
  # s <- raster::crop(s, extent_desert)
  names(s) <- band_names
  s <- s[[target_list]]
  n <- is.na(s)
  s[is.na(s)] <- 0
  s <- raster::setMinMax(s)
  #### if it's the output of the SAM algorithm, adjust results
  if(grepl("_sam", file_name, fixed = TRUE)){
    cat(paste0(file_name, " is SAM\n"))
    s <- abs(s-1)
    s[s>1]  <- NA
    s <- 1-s
  } else {
    cat(paste0(file_name, " is not SAM\n"))
  }
  ### Apply min-max normalization
  s_min <- raster::cellStats(s, min)
  s_max <- raster::cellStats(s, max)
  s <- (s-s_min)/s_max
  s[n] <- NA
  names(s) <- target_list
  all_stacks <- append(all_stacks, s)
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
geo/1000.0

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
