#takes in a peak list and calculates absorbance values for each peak using Beer-Lambert's equation
calculate_peak_absorbance = function(peak_list){
  #Initialize data frame to store peak data and calculated absorbance
  peak_absorbance_df = data.frame(Sample = character(0), Peak = numeric(0), RT = numeric(0), Intensity = numeric(0), Absorbance = numeric(0))
  
  #Iterate through all chromatograms in the peak list
  for(i in 1:length(peak_list)){
  
    #Storing chromatogram intensity data for Beer-Lambert calculations
    intensity_280 = get(attr(peak_list, "chrom_list"))[[peak_list[[i]][[1]][1,1]]][, 4]
    
    #Baseline intensity for chromatogram (minimum positive non-zero intensity at 280)
    ref_intensity = min(intensity_280[intensity_280 > 0])
    
    #Iterating over all peaks in current chromatogram
    for(j in 1:length(peak_list[[i]][[1]][,1])){
      #Extracting and storing peak data
      peak_num = rownames(peak_list[[i]][[1]])[j]
      peak_rt = peak_list[[i]][[1]][j, 3]
        
      #Peak intensity (intensity at peak RT)
      pk_intensity = peak_list[[i]][[1]][j, 9]
      
      #Calculating absorbance using Beer-Lambert's equation
      calc_absorbance = log10(ref_intensity / pk_intensity)
      
      #Creating a data frame with current data
      n_row = data.frame(Sample = names(peak_list)[i], Peak = peak_num, RT = peak_rt, Intensity = pk_intensity, Absorbance = calc_absorbance)
      rownames(n_row) = NULL
      
      #Add new row to data frame
      peak_absorbance_df = rbind(peak_absorbance_df, n_row)
    }
  }
  return(peak_absorbance_df)
}

calculate_peak_absorbance(standard_mix_1_peaks)
calculate_peak_absorbance(doe_pre_corrected_75_pks_egh)
