library(devtools)
source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/install.R")

# Import the GRS and HRV files into the current directory
    system2("dx", args = c("download", "/Alexander/nolte_mfi.txt"))
    system2("dx", args = c("download", "/Alexander/hrv_11_04_24.csv"))

# Read the HRV file as a dataframe

    hrv_df <- read.csv("hrv_11_04_24.csv")

# Generate GRS 
    a <- generate_grs('nolte_mfi.txt')

saveRDS(a, file = "GRS.rds")
system('dx upload GRS.rds --path /Alexander/')

