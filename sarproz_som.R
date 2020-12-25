source("LST/doe_utilities.R")





#### Main code -----------------
brady_csv_file <- ("../data/SARPROZ/Brady_72img_TS_7.csv")
brady_csv_data <- read.csv(brady_csv_file)
df_csv <- as.data.frame(brady_csv_data)

### Using extent_hymap to filter data


#### End Main code --------------