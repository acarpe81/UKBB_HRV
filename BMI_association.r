## Correlation of BMI (or any other primary care variable) with HRV (RMSSD)
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
    HR_clip_df <- merged_df %>% filter(mean_rr >= 0.6 & mean_rr <= 1.2)

# Load the UKBB health records library
    library(devtools) 
    source_url("https://raw.githubusercontent.com/hdg204/UKBB/main/UKBB_Health_Records_New_Project.R") 

# Reads GP data on BMI, DOB and age and joins with HRV data by eid
    a <- read_GP('22K..')
    df_hrv_bmi <- inner_join(a, HR_clip_df, by = "eid")
    df_hrv_bmi$event_age <- abs(df_hrv_bmi$event_age)


# Generate a frequency histogram for age using ggplot2
    library(ggplot2)
    ggplot(df_hrv_bmi, aes(x = event_age)) +
        geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
        labs(title = "Age Frequency Histogram", x = "Age", y = "Frequency") +
        scale_x_continuous(limits = c(0, max(df_hrv_bmi$event_age)))

# Save the merged dataframe to a CSV file and upload to DNAnexus
    write.csv(df_hrv_bmi, "df_hrv_bmi.csv", row.names = TRUE)
    system('dx upload df_hrv_bmi.csv --path /Alexander/')

# Generate model to perform regression analysis between RMSSD and age
    model <- glm(RMSSD ~ event_age, data = df_hrv_bmi)
    summary(model)

    F <- anova(model, test = "F")
    summary(F)

    # Calculate R^2
        # Calculate the total sum of squares with error handling
            sst <- sum((df_hrv_bmi$RMSSD[!is.na(df_hrv_bmi$RMSSD)] - mean(df_hrv_bmi$RMSSD[!is.na(df_hrv_bmi$RMSSD)]))^2)

        # Calculate the residual sum of squares
            sse <- sum(residuals(model)^2)

        # Calculate R-squared
            r_squared <- 1 - (sse / sst)

        print(r_squared)

# Generate model to perform regression analysis between RMSSD and BMI

# Convert the df_hrv_bmi columns value1, value2, and value3 to numeric
    df_hrv_bmi[, c("value1", "value2", "value3")] <- lapply(df_hrv_bmi[, c("value1", "value2", "value3")], as.numeric)

# Create a new column 'bmi' in df_hrv_bmi dataframe which takes the first non-NA numeric value from value1, value2, value3
    df_hrv_bmi$bmi <- ifelse(!is.na(df_hrv_bmi$value1) & is.numeric(df_hrv_bmi$value1), df_hrv_bmi$value1,
                            ifelse(!is.na(df_hrv_bmi$value2) & is.numeric(df_hrv_bmi$value2), df_hrv_bmi$value2,
                                ifelse(!is.na(df_hrv_bmi$value3) & is.numeric(df_hrv_bmi$value3), df_hrv_bmi$value3, NA)))

# Generate a frequency histogram for BMI using ggplot2 --Doesn't work for some unknown reason
    library(ggplot2)
    ggplot(df_hrv_bmi, aes(x = bmi)) +
        geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
        labs(title = "BMI Frequency Histogram", x = "BMI", y = "Frequency")

# Generate model to perform regression analysis between RMSSD and BMI
    model <- glm(RMSSD ~ bmi, data = df_hrv_bmi)
    summary(model)

    F <- anova(model, test = "F")
    summary(F)

# Calculate R^2
    sst <- sum((df_hrv_bmi$RMSSD[!is.na(df_hrv_bmi$RMSSD)] - mean(df_hrv_bmi$RMSSD[!is.na(df_hrv_bmi$RMSSD)]))^2)
    sse <- sum(residuals(model)^2)
    r_squared <- 1 - (sse / sst)
    print(r_squared)

        model <- glm(RMSSD ~ bmi, data = df_hrv_bmi)
        summary(model)

        F <- anova(model, test = "F")
        summary(F)

        # Calculate R^2
            # Calculate the total sum of squares with error handling
                sst <- sum((df_hrv_bmi$RMSSD[!is.na(df_hrv_bmi$RMSSD)] - mean(df_hrv_bmi$RMSSD[!is.na(df_hrv_bmi$RMSSD)]))^2)

            # Calculate the residual sum of squares
                sse <- sum(residuals(model)^2)

            # Calculate R-squared
                r_squared <- 1 - (sse / sst)

            print(r_squared)

# Association with sudden (cardiac) or unexplained death
    ICD10_codes=c('I46', 'I46.0', 'I46.1', 'I46.9', 'R96', 'R96.0', 'R96.1', 'R98', 'R99')
    ICD10_records=read_ICD10(ICD10_codes)

    df_ICD10 <- left_join(df_hrv_bmi, ICD10_records, by = "eid")

    library(dplyr)

    df_ICD10 <- df_ICD10 %>%
    mutate(binary = ifelse(!is.na(diag_icd10), 1, 0))

    model <- glm('binary~RMSSD',data=df_ICD10,family='binomial')
    summary(model)

    F <- anova(model, test = "F")
    summary(F)

