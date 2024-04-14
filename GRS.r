library(devtools)
source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/install.R")

# Import the file noltegrstxt2.txt into the current directory
    system2("dx", args = c("download", "/Alexander/noltegrstxt2.txt"))
    system2("dx", args = c("download", "/Alexander/hrv_11_04_24.csv"))

# Read the GRS and HRV files as dataframes
    noltegrs_df <- read.table("noltegrstxt2.txt", header = TRUE)
    hrv_df <- read.csv("hrv_11_04_24.csv")

# Generate GRS using eid column of hrv_df as input
    a <- generate_grs(hrv_df)
    #grs_list <- list(grs1 = a$grs, grs2 = a$grs)
    #hrv_df$grs1 <- a$grs
    #hrv_df$grs2 <- a$grs

# Print the updated dataframe
    hrv_df
