source("utils/doe_mineral_utils.R")

minerals_directory <- "../../doe_test/HyMapTargetDetectionFull"
band_names <- c("OpalizedTuff", "KrattOpal", "Epsomite", "Gypsum", 
                "Hematite", "Kaolinite", "Chalcedony")
target_list <- c("Chalcedony", "Kaolinite", "Gypsum", "Epsomite")

input_files <- c("HyMap_full_ace", "HyMap_full_cem", "HyMap_full_mf",
                 "HyMap_full_osp", "HyMap_full_sam_rule","HyMap_full_tcimf")

all_stacks <- vector(mode="list", length = 0)
for (file_name in input_files){
  s <- stack(file.path(minerals_directory, file_name))
  names(s) <- band_names
  s <- s[[target_list]]
  s <- raster::setMinMax(s)
  n <- is.na(s)
  if(grepl("sam_rule", file_name, fixed = TRUE, ignore.case = TRUE)){
    print(paste0(file_name, " is SAM"))
    s <- abs(s)
    s[s>1]  <- NA
    s <- 1-s
  }
  s_min <- raster::cellStats(s, min)
  s_max <- raster::cellStats(s, max)
  s <- (s-s_min)/s_max
  s[n] <- NA
  all_stacks <- append(all_stacks, s)
}

chalcedony <- stack()
for (i in 1:length(all_stacks))
  chalcedony <- stack(chalcedony, all_stacks[[i]][["Chalcedony"]])
chalcedony <- calc(chalcedony, sum, na.rm = TRUE)
chalcedony[is.na(all_stacks[[3]][["Chalcedony"]])] <- NA
plot(chalcedony)
ace_stack <- stack("results/filtered/ace_stack")
cem_stack <- stack("results/filtered/cem_stack")
mf_stack <- stack("results/filtered/mf_stack")
osp_stack <- stack("results/filtered/osp_stack")
sam_stack <- stack("results/filtered/sam_stack")
tcimf_stack <- stack("results/filtered/tcimf_stack")

chalcedony <- stack(sam_stack[["Chalcedony"]], ace_stack[["Chalcedony"]], cem_stack[["Chalcedony"]], mf_stack[["Chalcedony"]], osp_stack[["Chalcedony"]], tcimf_stack[["Chalcedony"]])
chalcedony <- calc(chalcedony, sum, na.rm = TRUE)

f1 <- doe_write_raster(chalcedony, "results/filtered/Brady_Chalcedony")
