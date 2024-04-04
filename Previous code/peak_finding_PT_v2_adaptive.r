# This is sample code for running analysis on a cardiac ECG in the UK Biobank in R You need to run the first two lines separately
install.packages("XML")
install.packages("signal")
install.packages("pracma")
library('dplyr')
library('XML')
library('ggplot2')
library('signal')

directories=system('ls ../../mnt/project/Bulk/Electrocardiogram/Resting', intern = TRUE)%>%as.numeric()
for (i in directories[1]){ #it's a bit slow, the 1 here and on line 14 are for testing
	path=paste0('../../mnt/project/Bulk/Electrocardiogram/Resting/',i)
	files=system(paste('ls ', path), intern = TRUE) # nolint
	for (j in files[1]){
		filepath=paste0(path,'/',j)
		xmldoc <- xmlTreeParse(filepath, useInternalNodes = TRUE)
		stripdata_node <- getNodeSet(xmldoc, "//StripData")[[1]]
		waveform_node <- getNodeSet(stripdata_node, "./WaveformData[@lead='II']")[[1]]
		waveform_text <- xmlValue(waveform_node)
		sample_rate <- 500	# Sampling frequency (Hz)
		signal <- as.numeric(strsplit(waveform_text, ",")[[1]])
		signal=signal-mean(signal)
		low_pass_cut=25 # Low pass filter
		high_pass_cut=0.5 # High pass filter
		filter_order=4 
		nyquist_frequency = 0.5 * sample_rate  # Assuming a sampling rate of 500 Hz
		butterworth_filter <- butter(filter_order, low_pass_cut/nyquist_frequency, type = "low")
		butterworth_filter_2 <- butter(filter_order, high_pass_cut/nyquist_frequency, type = "high")
		v_2 <- filter(butterworth_filter_2, signal)
		v_3 <- filter(butterworth_filter, v_2)
		ecg_frame=data.frame(t=1:5000/500,v=signal,v2=v_2,v3=v_3)
		# the nonlinear energy operator is a useful way to characterise the 'energy' of a signal, i.e. the time the most stuff is happening
		nleo=(abs((ecg_frame$v3[2:(nrow(ecg_frame)-1)])^2)-(((ecg_frame$v3[1:(nrow(ecg_frame)-2)])^2)*((ecg_frame$v3[3:(nrow(ecg_frame))])^2)))
		ecg_frame$energy=c(NA,nleo,NA)

		#	Start of Pan-Tompkins algorithm
		diff_signal=diff(ecg_frame$v3)	#Differentiating the signal to highlight the R wave upslope
	
		squared_signal=diff_signal^2	#Squaring the signal to emphasise larger differences
		
		#	Creating a moving window integrator, effectively operating as an averaging filter
		window_width <- round(0.08 * sample_rate)  # 80 ms window width							## Original value: 150 ms
		integrated_signal <- filter(rep(1/window_width, window_width), 1, squared_signal)

		#	Adaptive thresholding: # initialize thresholds based on the integrated signal
		thr_sig <- 0.2 * max(integrated_signal)
		thr_noise <- 0.5 * thr_sig

		# Initialize lists to keep track of detected R-peaks and noise peaks
		r_peaks <- c()
		noise_peaks <- c()

		# Initialize signal and noise levels
		signal_level <- thr_sig
		noise_level <- thr_noise

		# Loop through the integrated signal to detect peaks. The thresholds are adaptively updated based on the signal and noise levels.
		for (i in 2:(length(integrated_signal) - 1)) {
		if (integrated_signal[i] > thr_sig && integrated_signal[i] > integrated_signal[i-1] && integrated_signal[i] > integrated_signal[i+1]) {
			r_peaks <- c(r_peaks, i)

		# Update signal level and thr_sig based on the newly detected R-peak
		signal_level <- 0.875 * signal_level + 0.125 * integrated_signal[i]
		thr_sig <- noise_level + 0.25 * (signal_level - noise_level)
	} else if (integrated_signal[i] > thr_noise && integrated_signal[i] > integrated_signal[i-1] && integrated_signal[i] > integrated_signal[i+1]) {
		noise_peaks <- c(noise_peaks, i)
		
		# Update noise level and thr_noise based on the newly detected noise peak
		noise_level <- 0.875 * noise_level + 0.125 * integrated_signal[i]
		thr_noise <- noise_level + 0.25 * (signal_level - noise_level)
		}
		}
		
		# Post-processing: Remove peaks that are too close together
		min_peak_distance <- round(0.5 * sample_rate)  # Minimum RR 500 ms (120 bpm: to account for instantaneous HRV) 
													   
		# Initialize filtered peaks list with the first peak
		filtered_peak_indices <- r_peaks[1]
		for(i in 2:length(r_peaks)) {
		if(r_peaks[i] - filtered_peak_indices[length(filtered_peak_indices)] >= min_peak_distance) {
			filtered_peak_indices <- c(filtered_peak_indices, r_peaks[i])
				}
			}

		peaks_s <- filtered_peak_indices * 0.002	# Converts to seconds	
	}
}

# Plots peaks against ECG signal

ggplot(ecg_frame, aes(x = t))+
	geom_line(aes(y = v3, color = "Denoised Signal")) +
    geom_vline(xintercept = peaks_s,linetype = 2) +
	scale_color_manual(values = c("Denoised Signal" = "blue", "R intercept" = "black")) +
	theme_minimal()

# Creates an ordered list of R wave timings and calculates RR intervals and converts them to ms

RR_times <- peaks_s
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

## Further goals: 
	# 1. Generate a table of all the RMSSD values for each ECG file which iterates through all the files in the directory
	# 2. Include RMSSD, mean RR interval, heart rate in the table and exclude HR >100 or <50
	# 3. Generate a histogram of RMSSD values
	# 4. Generate a histogram of mean RR intervals
	# 5. Generate a histogram of heart rates
	# 6. Generate a histogram of the number of R peaks detected
	# 7. Include info on the file and directory in the table
	# 8. Include some patient identifiable label to allow UKBB to link back to the patient
	# 9. Include a marker for the presence of atrial fibrillation
	# 10. Include some measure of signal quality