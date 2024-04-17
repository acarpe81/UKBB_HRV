## Appends the HRV dataframe with the GRS data to allow regression analysis
## Author: Alexander Carpenter

library("dplyr")

# Import the GRS file as a dataframe
    system2("dx", args = c("download", "/Alexander/GRS.rds"))
    GRS <- readRDS("GRS.rds")

# Import then read the HRV file as a dataframe
    system2("dx", args = c("download", "/Alexander/hrv_11_04_24.csv"))
    hrv_df <- read.csv("hrv_11_04_24.csv")

# Perform inner join to merge the GRS and hrv_df dataframes based on matching eid
    merged_df <- inner_join(GRS$grs, hrv_df, by = "eid")

# Save the merged dataframe to a CSV file and upload to DNAnexus
    write.csv(merged_df, "hrv_grs.csv", row.names = TRUE)
    system('dx upload hrv_grs.csv --path /Alexander/')

# Generate model to perform regression analysis between RMSSD and GRS
    model <- glm(RMSSD ~ grs, data = merged_df)
    summary(model)

    F <- anova(model, test = "F")
    summary(F)