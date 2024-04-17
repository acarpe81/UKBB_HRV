## Generates a GRS from the MFI file and uploads it to the DNAnexus project
## Author: Alexander Carpenter

library(devtools)
source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/install.R")

# Import the GRS file into the current directory
    system2("dx", args = c("download", "/Alexander/nolte_mfi.txt"))

# Generate GRS 
    a <- generate_grs('nolte_mfi.txt')

saveRDS(a, file = "GRS.rds")
system('dx upload GRS.rds --path /Alexander/')
