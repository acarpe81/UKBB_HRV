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

# Create a new dataframe with rows that have mean_rr values between 600 and 1200 ms (HR 50-100) 
# to exclude non-physiological HRs and outliers caused by poor signal quality
    HR_clip_df <- merged_df %>% filter(mean_rr >= 0.6 & mean_rr <= 1.2)

    library(ggplot2)

# Create a frequency histogram for snr and snr_int values

    ggplot(HR_clip_df, aes(x = snr_db)) +
        geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
        labs(title = "Frequency Histogram of SNR_int Values",
            x = "SNR_int",
            y = "Frequency") +
        scale_x_continuous(breaks = seq(0, 20, 1))

    ggplot(HR_clip_df, aes(x = snr_db_int)) +
        geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
        labs(title = "Frequency Histogram of SNR_int Values",
                 x = "SNR_int",
                 y = "Frequency")

# Filter to exclude SNR_int values 
    HR_clip_df <- merged_df %>% filter(mean_rr >= 0.6 & mean_rr <= 1.2)
    
# Filter to include only rows with SNR values above 1 (exlcuding negative values)
    HRV_SNR <- HR_clip_df %>% filter(snr_db > 1)
    
# Save the merged dataframe to a CSV file and upload to DNAnexus
    write.csv(HR_clip_df, "HR_clip_df.csv", row.names = TRUE)
    system('dx upload HR_clip_df.csv --path /Alexander/')

# Generate model to perform regression analysis between RMSSD and GRS
    model <- glm(RMSSD ~ grs, data = HR_clip_df)
    summary(model)

    F <- anova(model, test = "F")
    summary(F)

# Calculate R^2
    # Calculate the total sum of squares with error handling
        sst <- sum((HR_clip_df$RMSSD[!is.na(HR_clip_df$RMSSD)] - mean(HR_clip_df$RMSSD[!is.na(HR_clip_df$RMSSD)]))^2)

    # Calculate the residual sum of squares
        sse <- sum(residuals(model)^2)

    # Calculate R-squared
        r_squared <- 1 - (sse / sst)

    print(r_squared)