## The final MR analysis - two-stage least squares (2SLS) regression
## Author: Alexander Carpenter supervised by Harry Green

library(dplyr)

# Load the UKBB health records library
    library(devtools) 
    source_url("https://raw.githubusercontent.com/hdg204/UKBB/main/UKBB_Health_Records_New_Project.R") 

    system2("dx", args = c("download", "/Alexander/df_hrv_death.csv"))
    df_hrv_death <- read.csv("df_hrv_death.csv")

# Add age

    a <- read_GP('22K..')
    a <- distinct(a, eid, .keep_all = TRUE)
    df_hrv_death <- left_join(df_hrv_death, a, by = "eid")
    df_hrv_death$event_age <- abs(df_hrv_death$event_age)

# Add BMI
    # Convert the df_hrv_death columns value1, value2, and value3 to numeric
    df_hrv_death[, c("value1", "value2", "value3")] <- lapply(df_hrv_death[, c("value1", "value2", "value3")], as.numeric)

    # Create a new column 'bmi' in df_hrv_death dataframe which takes the first non-NA numeric value from value1, value2, value3
    df_hrv_death$bmi <- ifelse(!is.na(df_hrv_death$value1) & is.numeric(df_hrv_death$value1), df_hrv_death$value1,
                                ifelse(!is.na(df_hrv_death$value2) & is.numeric(df_hrv_death$value2), df_hrv_death$value2,
                                    ifelse(!is.na(df_hrv_death$value3) & is.numeric(df_hrv_death$value3), df_hrv_death$value3, NA)))

# MR analysis with all cause mortality- two-stage least squares (2SLS) regression

    # Stage 1: GLM - exposure variable regressed on the genetic risk score (instrumental variable) and other covariates
        stage_1 = glm(df_hrv_death$RMSSD ~ df_hrv_death$grs + df_hrv_death$event_age)

    # Creates new data frame including the predicted values from the stage 1 model for each individual (eid)
        stage_1_frame=data.frame(eid=df_hrv_death$eid,pred=predict(stage_1,df_hrv_death))

    # Merges the new df with the original data, matching by eid
        stage_1_frame=stage_1_frame%>%inner_join(df_hrv_death)

    # Second stage of 2SlS regression - GLM where outcome variable regressed on the predicted values from the first stage 
    # and other covariates
        stage_2 = glm(has_death ~ pred + event_age, data = stage_1_frame, family=binomial)%>%summary()

# MR analysis with sudden or unexplained death

    # Association with sudden (cardiac) or unexplained death
        ICD10_codes=c('I46', 'I46.0', 'I46.1', 'I46.9', 'R96', 'R96.0', 'R96.1', 'R98', 'R99')
        ICD10_records=read_ICD10(ICD10_codes)

        df_hrv_icd <- left_join(df_hrv_death, ICD10_records, by = "eid")
       
        df_hrv_icd$has_ICD_death <- ifelse(!is.na(df_hrv_icd$date_of_death), 1, 0)
        
        #df_ICD10 <- df_ICD10 %>%
        #mutate(binary = ifelse(!is.na(diag_icd10), 1, 0))

        model <- glm('has_ICD_death ~ RMSSD', data = df_hrv_icd, family = 'binomial')
        summary(model)

        F <- anova(model, test = "F")
        summary(F)


    # Stage 1: GLM - exposure variable regressed on the genetic risk score (instrumental variable) and other covariates
        stage_1 = glm(df_hrv_icd$RMSSD ~ df_hrv_icd$grs + df_hrv_icd$event_age)

    # Creates new data frame including the predicted values from the stage 1 model for each individual (eid)
        stage_1_frame = data.frame(eid = df_hrv_icd$eid, pred = predict(stage_1, df_hrv_icd))

    # Merges the new df with the original data, matching by eid
        stage_1_frame = stage_1_frame%>%inner_join(df_hrv_icd)

    # Second stage of 2SlS regression - GLM where outcome variable regressed on the predicted values from the first stage 
    # and other covariates
        stage_2 = glm(has_ICD_death ~ pred + event_age, data = stage_1_frame, family=binomial)%>%summary()