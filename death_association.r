## Correlation of death (all-cause) with HRV (RMSSD)
## Author: Alexander Carpenter

library("dplyr")

# Import the GRS file as a dataframe
    system2("dx", args = c("download", "/Alexander/GRS.rds"))
    GRS <- readRDS("GRS.rds")

# Import then read the HRV file as a dataframe
    system2("dx", args = c("download", "/Alexander/hrv_11_04_24.csv"))
    hrv_df <- read.csv("hrv_11_04_24.csv")

# Perform inner join to merge the GRS and hrv_df dataframes based on matching eid
    L_merged_df <- left_join(hrv_df, GRS$grs, by = "eid")

# Create a new dataframe with rows that have mean_rr values between 600 and 1200 ms (HR 50-100) 
# to exclude non-physiological HRs and outliers caused by poor signal quality
    HR_clip_df <- L_merged_df %>% filter(mean_rr >= 0.6 & mean_rr <= 1.2)

# Load the UKBB health records library
    library(devtools) 
    source_url("https://raw.githubusercontent.com/hdg204/UKBB/main/UKBB_Health_Records_New_Project.R") 

# Load death data
    system('dx download file-GZKXVx8J9jFp1qpBQZ8z5PbJ')
    death <- read.csv("death_death.csv")
    df_hrv_death <- left_join(HR_clip_df, death, by = "eid")
    df_hrv_death$has_death <- ifelse(!is.na(df_hrv_death$date_of_death), 1, 0)

#Save and upload the merged dataframe to DNAnexus
    write.csv(df_hrv_death, "df_hrv_death.csv", row.names = TRUE)
    system('dx upload df_hrv_death.csv --path /Alexander/')

# GLM model to perform regression analysis between RMSSD and death
    model <- glm('has_death ~ RMSSD',data=df_hrv_death,family='binomial')
    summary(model)
    # Calculate R-squared value
    r_squared <- summary(model)$r.squared
    r_squared

        F <- anova(model, test = "F")
    summary(F)

# RMSSD correction for mean HR (RR interval) possible with the following code:
    df_hrv_death$RMSSD_normalised <- df_hrv_death$RMSSD / df_hrv_death$mean_rr


    # Quantify the correlation between RMSSD and has_death
    correlation <- cor(df_hrv_death$RMSSD, df_hrv_death$has_death, use = "complete.obs")
    correlation
    
    correlation_test <- cor.test(df_hrv_death$RMSSD, df_hrv_death$has_death, use = "complete.obs")
    correlation_test$p.value

# Regression analysis between GRS and death
    model <- glm('has_death ~ grs',data=df_hrv_death,family='binomial')
    summary(model)


