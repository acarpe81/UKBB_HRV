# This is sample code for running analysis on a cardiac ECG in the UK Biobank in R You need to run the first two lines separately
install.packages("XML")
install.packages("signal")
install.packages("pracma")
library('dplyr')
library('XML')
library('ggplot2')
library('signal')

#This paragraph is an early attempt at creating some sort of persistent list to output to, although I haven't thought it through yet!
#Label_Directory <- c("1", "2")
#Label_File <- c("1", "2")
#Label_RMSSD <- c("1", "2")
#Summary <- data.frame(Label_Directory,Label_File,Label_RMSSD)

directories=system('ls ../../mnt/project/Bulk/Electrocardiogram/Resting', intern = TRUE)%>%as.numeric()
for (i in directories[1]){ #it's a bit slow, the 1 here and on line 14 are for testing
	path=paste0('../../mnt/project/Bulk/Electrocardiogram/Resting/',i)
	files=system(paste('ls ', path), intern = TRUE)
	for (j in files){
		filepath=paste0(path,'/',j)
		xmldoc <- xmlTreeParse(filepath, useInternalNodes = TRUE)
		stripdata_node <- getNodeSet(xmldoc, "//StripData")[[1]]
		waveform_node <- getNodeSet(stripdata_node, "./WaveformData[@lead='II']")[[1]]
		waveform_text <- xmlValue(waveform_node)
		signal <- as.numeric(strsplit(waveform_text, ",")[[1]])
		signal=signal-mean(signal)
		low_pass_cut=25 #reduce this number to strengthen the denoising
		high_pass_cut=1 # change this number to alter the detrending
		filter_order=4 # I forgot what this does, I haven't done signal processing since forever
		nyquist_frequency = 0.5 * 500  # Assuming a sampling rate of 500 Hz
		butterworth_filter <- butter(filter_order, low_pass_cut/nyquist_frequency, type = "low")
		butterworth_filter_2 <- butter(filter_order, high_pass_cut/nyquist_frequency, type = "high")
		v_2 <- filter(butterworth_filter_2, signal)
		v_3 <- filter(butterworth_filter, v_2)
		ecg_frame=data.frame(t=1:5000/500,v=signal,v2=v_2,v3=v_3)
		# the nonlinear energy operator is a useful way to characterise the 'energy' of a signal, i.e. the time the most stuff is happening
		nleo=(abs((ecg_frame$v3[2:(nrow(ecg_frame)-1)])^2)-(((ecg_frame$v3[1:(nrow(ecg_frame)-2)])^2)*((ecg_frame$v3[3:(nrow(ecg_frame))])^2)))
		ecg_frame$energy=c(NA,nleo,NA)

            # Peak finding

            x <- ecg_frame$t
            sample_rate <- 500

            peaks<-
            data.frame(
                pracma::findpeaks(ecg_frame$v3[1:nrow(ecg_frame)], 
                                minpeakdistance = .4*sample_rate,
                                minpeakheight=mean(ecg_frame$v3[1:nrow(ecg_frame)])+1.85*sd(ecg_frame$v3[1:nrow(ecg_frame)])
                )
            )
            peaks_s <- peaks * 0.002
            
            # Creates an ordered list of R wave timings and calculates RR intervals and converts them to ms

            RR_times <- peaks_s$X2
            RR_order <- RR_times[order(RR_times)]
            RR_times_plus1 <- RR_order[2:length(RR_order)]
            RR_order_trim <- head(RR_order,-1)
            RR_int <- RR_times_plus1 - RR_order_trim
            RR_ms <- RR_int * 1000

            # Starting to calculate the RMSSD - calculates differences between RR intervals, squares and averages them

            RR_ms_offset1 <- RR_ms[2:length(RR_ms)]
            RR_ms_trim <- head(RR_ms,-1)
            RR_diff <- RR_ms_offset1 - RR_ms_trim

            RR_diff_squ <- RR_diff^2
            RR_diff_squ_avg <- mean(RR_diff_squ)
            RMSSD <- sqrt(RR_diff_squ_avg)
            RMSSD

            #Summary$Label_Directory[i] <-i
            #Summary$Label_File <- j
            #Summary$Label_RMSSD <- RMSSD
        }
    }

