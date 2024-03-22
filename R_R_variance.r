# Calculate RR intervals
filtered_peak_indices_s <- filtered_peak_indices / 500 # Convert indices to seconds
rr_intervals <- diff(filtered_peak_indices_s)

# Calculate the mean RR interval (could also set a predefined expected minimum RR interval)
mean_rr <- mean(rr_intervals)

# Instead of using standard deviation, define an explicit threshold for "short" RR intervals
# For instance, any RR interval less than 80% of the mean RR might be considered short
short_rr_threshold <- 0.8 * mean_rr

# Identify short RR intervals that are likely due to premature ectopic beats
short_rr_indices <- which(rr_intervals < short_rr_threshold)

# Flag the R-peaks leading to short RR intervals for exclusion
# Note: Adding 1 to the indices because diff() reduces the length by 1
ectopic_r_peaks <- filtered_peak_indices_s[short_rr_indices + 1]

# Remove the ectopic R-peaks from the original list of R-peaks
r_peaks_filtered <- filtered_peak_indices_s[!filtered_peak_indices_s %in% ectopic_r_peaks]

# r_peaks_filtered now contains the indices of R-peaks excluding those leading to short RR intervals

# Optional: Plotting for verification (assuming a mock ECG signal for demonstration)


plot(ecg_frame$t, ecg_frame$V3, type = "l", xlab = "Time", ylab = "ECG Signal", main = "Filtered R-Peaks Excluding Premature Ectopic Beats")
points(filtered_peak_indices_s, ecg_frame$V3[filtered_peak_indices_s], col = "red", pch = 19, cex = 1.5, lwd = 2)
points(r_peaks_filtered, ecg_frame$V3[r_peaks_filtered], col = "blue", pch = 17, cex = 1, lwd = 2)
legend("topright", legend = c("Original R-Peaks", "Filtered R-Peaks"), col = c("red", "blue"), pch = c(19, 17))



#OR


# Assuming the sampling frequency (fs) is 500 Hz
fs <- 500  # Hz, samples per second

# Correct the time scale for plotting
# Since the original r_peaks are indices, convert these to time in seconds
time_r_peaks <- filtered_peak_indices / fs  # Convert sample indices to time in seconds
time_r_peaks_filtered <- r_peaks_filtered / fs  # Convert sample indices to time in seconds

# Assuming ecg_signal represents the amplitude values of your ECG data
# Recalculate the time_points vector to accurately reflect time in seconds
time_points <- seq(from = 0, to = length(ecg_signal) / fs, length.out = length(ecg_signal))

# Plotting the ECG signal with original and filtered R-peaks corrected for time scale
plot(time_points, ecg_signal, type = "l", xlab = "Time (seconds)", ylab = "ECG Signal", main = "ECG Signal and R-Peaks")
points(time_r_peaks, ecg_signal[r_peaks], col = "red", pch = 19, cex = 1.5, lwd = 2, main = "Original vs. Filtered R-Peaks")
points(time_r_peaks_filtered, ecg_signal[r_peaks_filtered], col = "blue", pch = 17, cex = 1, lwd = 2)
legend("topright", legend = c("Original R-Peaks", "Filtered R-Peaks"), col = c("red", "blue"), pch = c(19, 17))

