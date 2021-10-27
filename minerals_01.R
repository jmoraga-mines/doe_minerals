source("utils/doe_mineral_utils.R")

minerals_directory <- "../../doe_test/HyMapTargetDetectionFull"
# All thresholding methods: "OpalizedTuff_CU00_15E", "KrattOpal",
#    "Epsomite_GDS149", "Gypsum_SU2202", "Hematite_FE2602", "Kaolinite_CM3",
#    "Chalcedony_CU91_6a"
band_names <- c("OpalizedTuff", "KrattOpal", "Epsomite", "Gypsum", 
                "Hematite", "Kaolinite", "Chalcedony")
target_list <- c("Chalcedony", "Kaolinite", "Gypsum", "Epsomite")


# Create results directory if it does not exist
if (!dir.exists("results"))
  dir.create("results", showWarnings = FALSE, recursive = TRUE, mode = "0775")

# Save shapefiles with boundaries for both sites
raster::shapefile(poly_extent_brady, 
                  filename="results/extent_brady_frame.shp", 
                  overwrite=TRUE)
raster::shapefile(poly_extent_desert, 
                  filename="results/extent_desert_frame.shp", 
                  overwrite=TRUE)

# Collapse multiple layers into one by summing
brady_minerals_s <- raster_euclidean_collapse(sam)

raster::spplot(brady_minerals_s, main = "Brady SAM Euclidean Collapse (relaxed)", 
               col.regions=rev(heat.colors(100))) + 
  latticeExtra::layer(sp.polygons(poly_extent_brady, lwd=1, col = "blue")) + 
  latticeExtra::layer(sp.polygons(poly_extent_desert, lwd=1, col = "blue"))

# Normalize
limit <- raster::cellStats(brady_minerals_s, max)
brady_minerals_euclidean <- brady_minerals_s / limit

brady_sam <- raster::crop(brady_minerals_euclidean, extent_tall_brady)
plot(brady_sam, main="Brady SAM complete (normalized)", 
     axes=TRUE,
     col=rev(heat.colors(100)))
desert_sam <- raster::crop(brady_minerals_euclidean, extent_desert)
plot(desert_sam, main="Desert Peak SAM complete (normalized)", 
       axes=TRUE,
       col=rev(heat.colors(100)))

#b <- as.matrix(brady_minerals_euclidean) * 1000
#mode(b) <- "integer"
threshold <- select_largest_threshold(brady_minerals_euclidean, ignore_black = FALSE)
b_t <- threshold$th
th_method <- threshold$method
top_sam <- raster::crop(brady_minerals_euclidean, extent_tall_brady)
top_sam <- brady_minerals_euclidean
top_sam[top_sam < b_t] <- 0
top_sam <- unit_normalization(top_sam)
plot(top_sam, main = paste0("SAM Brady Threshold (", th_method, ")"), 
     axes=TRUE,
     col=rev(heat.colors(100)))

clean_sam <- top_sam
rc <- raster::clump(clean_sam, directions=8)
f <- as.data.frame(freq(rc))
invalid_clumps <- f[which(f$count<3),]$value
rc[invalid_clumps]<-NA

#####################################################################

#####################################################################
# Create results directory if it does not exist
if (!dir.exists("results/SAM"))
  dir.create("results/SAM", recursive = TRUE, mode = "0775")
dir.create("results/filtered", showWarnings = FALSE, recursive = TRUE, mode = "0775")

cem_stack <- full_detection(base_dir = minerals_directory, 
                            file_name = "HyMap_full_cem",
                            band_names = band_names,
                            target_list = target_list,
                            e = extent_tall_brady,
                            is_sam = FALSE, ignore_black = FALSE)
ace_stack <- full_detection(base_dir = minerals_directory, 
                            file_name = "HyMap_full_ace",
                            band_names = band_names,
                            target_list = target_list,
                            e = extent_tall_brady,
                            is_sam = FALSE, ignore_black = TRUE)

osp_stack <- full_detection(base_dir = minerals_directory, 
                            file_name = "HyMap_full_osp",
                            band_names = band_names,
                            target_list = target_list,
                            e = extent_tall_brady,
                            is_sam = FALSE, ignore_black = FALSE)

sam_stack <- full_detection(base_dir = minerals_directory, 
                            file_name = "HyMap_full_sam_rule",
                            band_names = band_names,
                            target_list = target_list,
                            e = extent_tall_brady,
                            is_sam = TRUE, ignore_black = FALSE)

mf_stack <- full_detection(base_dir = minerals_directory, 
                            file_name = "HyMap_full_mf",
                            band_names = band_names,
                            target_list = target_list,
                            e = extent_tall_brady,
                            is_sam = FALSE, ignore_black = FALSE)

tcimf_stack <- full_detection(base_dir = minerals_directory, 
                            file_name = "HyMap_full_tcimf",
                            band_names = band_names,
                            target_list = target_list,
                            e = extent_tall_brady,
                            is_sam = FALSE, ignore_black = FALSE)

f1 <- doe_write_raster(cem_stack, "results/filtered/cem_stack")
f1 <- doe_write_raster(ace_stack, "results/filtered/ace_stack")
f1 <- doe_write_raster(osp_stack, "results/filtered/osp_stack")
f1 <- doe_write_raster(mf_stack, "results/filtered/mf_stack")
f1 <- doe_write_raster(sam_stack, "results/filtered/sam_stack")
f1 <- doe_write_raster(tcimf_stack, "results/filtered/tcimf_stack")
rm(f1)

names(cem_stack) <- band_names
cem_stack <- cem_stack[[c("Chalcedony", "Kaolinite", "Gypsum",
                          "Epsomite")]]
raster::setMinMax(cem_stack)
cem_stack <- unit_normalization(get_stack)

spplot(ace_stack)
