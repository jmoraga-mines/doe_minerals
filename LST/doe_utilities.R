#### Utility functions, used for DOE project #
####

if (!require("pacman")) {
  install.packages("pacman")
  require("pacman")
}
p_load(raster, rgdal)

# Landsat-8 LST files use Albers projection in meters
# Resolution is 30 by 30 mts
extent_hymap_lst <-
  raster(
    xmn = 325209,
    xmx = 336531,
    ymn = 4395438,
    ymx = 4412103,
    res = c(30, 30),
    crs = "+init=epsg:32611"
  )
extent_hymap <-
  raster(
    xmn = 325209,
    xmx = 336531,
    ymn = 4395438,
    ymx = 4412103,
    res = c(3, 3),
    crs = "+init=epsg:32611"
  )

doe_write_raster <-
  function(x,
           filename,
           format = "raster",
           overwrite = TRUE,
           bandorder = "BSQ") {
    if (tools::file_ext("filename.grd") != "grd") {
      filename <- tools::file_path_sans_ext(filename)
      filename <- paste(filename, ".grd", sep = "")
    }
    f1 <- writeRaster(
      x = x,
      filename = filename,
      bandorder = bandorder,
      format = format,
      overwrite = overwrite
    )
    hdr(f1, "ENVI")
    return(f1)
  }


coords2spatial <-
  function(coord,
           sp_data,
           projection = crs("+init=epsg:32611")) {
    print("Creating Points")
    d0 <- SpatialPoints(coord, proj4string = projection)
    print("converting")
    sp_data <- as.data.frame(sp_data)
    print("Assigning values")
    d0$sp_data <- sp_data
    return(d0)
  }

coords2spatialdf <-
  function(coord,
           sp_data,
           projection = crs("+init=epsg:32611")) {
    d0 <-
      SpatialPointsDataFrame(coord, sp_data, proj4string = projection)
    return(d0)
  }

cluster_fun <- function(df_data) {
  clusters <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(clusters) <-
    c("image", "avg", "min", "max", "stddev", "range")
  for (i in names(df_data)) {
    sub_cluster <- df_data[i]
    s_mean <- sapply(sub_cluster, mean, na.rm = TRUE)
    s_min <- sapply(sub_cluster, min, na.rm = TRUE)
    s_max <- sapply(sub_cluster, max, na.rm = TRUE)
    s_sd  <- sapply(sub_cluster, sd, na.rm = TRUE)
    clusters[i, 1:6] <-
      list(i, s_mean, s_min, s_max, s_sd, s_max - s_min)
  }
  return(clusters)
}

df_plot <- function(df_img, save_file = FALSE) {
  plt_img <- df_img
  coordinates(plt_img) <- ~ x + y
  proj4string(plt_img)  <- "+init=epsg:32611"
  gridded(plt_img) <- TRUE
  plot(plt_img)
  if (save_file != FALSE) {
    writeRaster(
      raster(plt_img),
      filename = save_file,
      format = "GTiff",
      overwrite = TRUE
    )
  }
  return(plt_img)
}

get_subset <- function(a_df_lst, a_df_mask, img_name, a_class) {
  return(a_df_lst[(a_df_mask[img_name] == a_class), c("x", "y", img_name)])
}

save_plot <- function(a_raster, a_name) {
  png(paste0(a_name, ".png"), width = 600, height = 800)
  print(spplot(
    a_raster[[a_name]],
    col.regions = c("Yellow", "Cyan", "Magenta",
                    "Green", "Red", "Black"),
    main = paste0("K-means clustering (", a_name, ")")
  ))
  if (FALSE) {
    plot(a_raster[[a_name]],
         col = c("Black", "Yellow", "Cyan",
                 "Magenta", "Green", "Red"))
  }
  dev.off()
}

