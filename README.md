#### DOE Project: DE-EE0008760

```
# This code makes use of the CSV output files from SARPROZ
#
# CSV file fields:
#
# ID - point reference id (sequential)
# LAT, LON - latitude and longitude
# SVET, LVET – sample and line, pixel SAR coordinates
# HEIGHT – relative to sea level (m)
# HEIGHT WRT DEM – relative to the external DEM (m)
# SIGMA HEIGHT – height estimation standard deviation
# VEL – deformation velocity in mm/year
# SIGMA VEL – velocity standard deviation
# SEASONAL - seasonal trend coefficient
# CUMUL.DISP. - point cumulative displacement (mm)
# COHER - tepomral coherence
# IN, FIN - target point on-off model
# STDEV - phase standard deviation from estimated model
#
# psi_std_headers: ("ID", "X", "Y", "LAT", "LON", "HEIGHT", "HEIGHT.WRT.DEM",
#                    "SIGMA.HEIGHT", "VEL", "SIGMA.VEL", "SEASONAL",
#                    "CUMUL.DISP.", "COHER", "SVET", "LVET", "IN",
#                    "FIN", "STDEV")
#
# Input data is generated in SARPROZ by using SARPROZ Sparse Data Export
# [Excerpt from https://sarproz.com/manual/dbf_sito.html]
#
# With this function, the tool exports data on a point-wise basis.
# Data can be exported in 3 main formats. Choose:
# - csv files (text tables in csv format)
#
# Data are color-coded according to the selected parameter ("Color Code").
# Look at "Load Mask" for the full list of parameters.
# The color can be saturated passing minimum and maximum values.
# A spatial linear trend can be removed by selecting the corresponding check-box
# A reference point can be selected as well (it should have been saved by the
# "APS estimation" module).
# Default reference point is automatically selected.
# By clicking on "No Orthorect." you can geocode points (convert from SAR to
# geographic coordinates) without using the estimate height. This is useful 
# in cases where the accuracy of height estimation is low
# (e.g. in case of Sentinel data) and if you want to have your points aligned
# with the geocoded SAR images/data.
#
# Finally, time series ***must*** be included in the exported data by selecting
# the corresponding check-box.
# The time series export operation is available only if MISP processing options
# are already loaded. It means, if after MISP you access this module, the 
# software will use the same options to generate and export time series.
# If no processing options are loaded, use the "time series module" to generate
# and export time series.
```