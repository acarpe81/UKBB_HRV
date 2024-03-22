# Load necessary libraries
install.packages(c("XML", "signal", "pracma", "dplyr", "ggplot2"))
library(XML)
library(signal)
library(pracma)
library(dplyr)
library(ggplot2)

# Define ECG signal processing and R peak detection function
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
		signal <- as.numeric(strsplit(waveform_text, ",")[[1]])
		signal=signal-mean(signal)

		
process_ecg_signal <- function(filepath) {
  # Read and preprocess ECG signal from XML
  xmldoc <- xmlTreeParse(filepath, useInternalNodes = TRUE)
  waveform_node <- getNodeSet(xmldoc, "//WaveformData[@lead='II']")[[1]]
  signal <- as.numeric(strsplit(xmlValue(waveform_node), ",")[[1]])
  signal <- signal - mean(signal) # Detrend
  
  # Filter design
  fs <- 500 # Sampling frequency
  nyquist <- fs / 2
  low_pass_cut <- 25 / nyquist
  high_pass_cut <- 0.5 / nyquist # Slightly adjust to avoid very low-frequency drift
  filter_order <- 4
  
  # Apply filters
  bp_filter <- butter(filter_order, c(high_pass_cut, low_pass_cut), type = "band")
  signal_filtered <- filtfilt(bp_filter, signal)
  
  # R Peak detection (simplified for clarity)
  # Here, we can use a more sophisticated algorithm such as Pan-Tompkins, but let's keep it straightforward
  peak_indices <- findpeaks(signal_filtered, minpeakheight = quantile(signal_filtered, 0.98), 
                            minpeakdistance = 0.2*fs)$loc
  
  # Compute RR intervals
  RR_intervals <- diff(peak_indices) / fs # Convert to seconds
  HRV_RMSSD <- sqrt(mean(diff(RR_intervals)^2)) # Calculate RMSSD for HRV
  
  list(peaks = peak_indices, RR_intervals = RR_intervals, HRV_RMSSD = HRV_RMSSD)
}

# Example usage
# Assuming 'filepath' is the path to your ECG XML file
# results <- process_ecg_signal(filepath)
# plot(signal, type = "l")
# points(results$peaks, signal[results$peaks], col = "red", pch = 19)

# Assuming signal_filtered and peak_indices are obtained from the previous processing steps

calculate_snr <- function(signal_filtered, peak_indices) {
  # Signal Power: Calculate the power of the R-peaks
  # Assuming R-peaks are significantly higher than the noise, thus taking them as the signal
  r_peak_values <- signal_filtered[peak_indices]
  signal_power <- mean(r_peak_values^2)
  
  # Noise Power: Calculate the power of the signal excluding R-peaks
  # This is a simplistic approach - consider enhancing for more accurate noise estimation
  noise_indices <- setdiff(1:length(signal_filtered), peak_indices)
  noise_power <- mean(signal_filtered[noise_indices]^2)
  
  # Calculate SNR in decibels (dB)
  snr_db <- 10 * log10(signal_power / noise_power)
  
  return(snr_db)
}

# Example usage:
# Assuming 'signal_filtered' is your filtered ECG signal and 'peak_indices' are the detected R-peak indices
# snr <- calculate_snr(signal_filtered, peak_indices)
# print(snr)