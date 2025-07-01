cl = makeCluster(detectCores() - 1)

#Subset out standard mix runs from 2020 standards by using logical indexing on attributes of runs
standard_mixes_2020 = hplc_list_standards[sapply(hplc_list_standards, function(x) startsWith(attr(x, "sample_id"), "Standard"))]
print(length(standard_mixes_2020)) #should be 13 according to alignments file

#Subset out "check" runs from 2020 standards
standard_checks_2020 = hplc_list_standards[sapply(hplc_list_standards, function(x){
  startsWith(attr(x, "sample_id"), "Check") || startsWith(attr(x, "sample_id"), " Check")})] #quick fix because a check run had a space before its name
print(length(standard_checks_2020)) #should be 27 according to alignments file

#preprocess standard mixes and check runs
standard_mixes_2020_pre = preprocess(standard_mixes_2020, dim1 = seq(0, 54, 0.01), dim2 = seq(190, 400, 1), 
                                     remove.time.baseline = TRUE, cl = cl)
standard_checks_2020_pre = preprocess(standard_checks_2020, dim1 = seq(0, 54, 0.01), dim2 = seq(190, 400, 1), 
                                     remove.time.baseline = TRUE, cl = cl)
#plot standard mixes and checks
plot_chroms(standard_mixes_2020_pre, lambdas = c("280"), show_legend = FALSE)
plot_chroms(standard_checks_2020_pre, lambdas = c("280"), show_legend = FALSE)

#divide preprocessed and corrected 2020 standard mixes into 3 standard mix groupings
standard_mix_2020_1_pre = standard_mixes_2020_pre[1:4]
standard_mix_2020_2_pre = standard_mixes_2020_pre[5:8]
standard_mix_2020_3_pre = standard_mixes_2020_pre[10:13]

#correct each
standard_mix_2020_1_pre = correct_rt(standard_mix_2020_1_pre, lambdas = c('280'), what = c("corrected.values"), cl = cl)
standard_mix_2020_2_pre = correct_rt(standard_mix_2020_2_pre, lambdas = c('280'), what = c("corrected.values"), cl = cl)
standard_mix_2020_3_pre = correct_rt(standard_mix_2020_3_pre, lambdas = c('280'), what = c("corrected.values"), cl = cl)
#get peak lists for 3 2020 standard mixes
#ampthresh tweaked targetting 7 unique peaks
plot_chroms(standard_mix_2020_1_pre, lambdas = c('280'), show_legend = FALSE)

standard_mix_2020_1_peaks = get_peaks(standard_mix_2020_1_pre, lambdas = c('280'), fit = "egh", amp_thresh = 1900)
print(standard_mix_2020_1_peaks)
#peak finding on 2020 standard mix 1 at all concentrations is difficult due to increasing amplitudes at higher 
#concentrations preventing proper peak filtering with the amp_thresh parameter

#SOLUTION: split into individual concentrations from pre_processed chrom_list, and get peaks for each individual 
#concentration
standard_mix_2020_1_1mg = standard_mix_2020_1_pre[1]
standard_mix_2020_1_1mg_peaks = get_peaks(standard_mix_2020_1_1mg, lambdas = c('280'), fit = "egh",
                                          smooth_window = 0.001, amp_thresh = 10000)


standard_mix_2020_1_p5mg = standard_mix_2020_1_pre[2]
standard_mix_2020_1_p5mg_peaks = get_peaks(standard_mix_2020_1_p5mg, lambdas = c('280'), fit = "egh", 
                                           smooth_window = 0.001, amp_thresh = 7000)


standard_mix_2020_1_p25mg = standard_mix_2020_1_pre[3]
standard_mix_2020_1_p25mg_peaks = get_peaks(standard_mix_2020_1_p25mg, lambdas = c('280'), fit = "egh",
                                            smooth_window = 0.0025, amp_thresh = 3800)


standard_mix_2020_1_p125mg = standard_mix_2020_1_pre[4]
standard_mix_2020_1_p125mg_peaks = get_peaks(standard_mix_2020_1_p125mg, lambdas = c('280'), fit = "egh", 
                                           smooth_window = 0.004, sd.max = 10.4, amp_thresh = 1639)
rownames(standard_mix_2020_1_p125mg_peaks[[1]]$`280`) = c("1","2","3","4","5","6","7")


#for 2020 standard mix 2
#should find 5 peaks
plot_chroms(standard_mix_2020_2_pre, lambdas = c('280'), show_legend = FALSE)

standard_mix_2020_2_peaks = get_peaks(standard_mix_2020_2_pre, lambdas = c('280'), fit = "egh", amp_thresh = 100000)
print(standard_mix_2020_2_peaks)

#no issues with ampthresh here, so we can get all the data needed from standard_mix_2020_2_peaks


#manually compile concentration and absorbance data for each external standard in 2020 mix 1 to make linear regression 
#models for calculating unknown concentrations 
#store known concentrations as a vector
known_conc = c(1, 0.5, 0.25, 0.125)

#store absorbance data in order of known concentrations
#Gallic Acid (rt of 5.42 from external standard run)
smix_2020_gallic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[1], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[1], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[1],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[1])
#fit a linear model
smix_2020_gallic_acid_calibration = lm(smix_2020_gallic_acid_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_gallic_acid_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_gallic_acid_calibration, which = 3)

#Chlorogenic Acid (rt of 20.17 from external standard run)
smix_2020_chlorogenic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[2], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[2], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[2],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[2])
#fit a linear model
smix_2020_chlorogenic_acid_calibration = lm(smix_2020_chlorogenic_acid_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_chlorogenic_acid_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_chlorogenic_acid_calibration, which = 3)


#Vanillic Acid (rt of 22.29 from external standard run)
smix_2020_vanillic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[3], 
                                          standard_mix_2020_2_peaks[[2]]$'280'$height[3], 
                                          standard_mix_2020_2_peaks[[3]]$'280'$height[3],
                                          standard_mix_2020_2_peaks[[4]]$'280'$height[3])
#fit a linear model
smix_2020_vanillic_acid_calibration = lm(smix_2020_vanillic_acid_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_vanillic_acid_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_vanillic_acid_calibration, which = 3)


#Eluetheroside B (rt of 23.87 from external standard run)
smix_2020_eleuth_b_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[1], 
                                    standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[1], 
                                    standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[1],
                                    standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[1])
#fit a linear model
smix_2020_eleuth_b_calibration = lm(smix_2020_eleuth_b_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_eleuth_b_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_eleuth_b_calibration, which = 3)


#Coumaric Acid (rt of 32.16 from external standard run)
smix_2020_coumaric_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[4], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[4], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[4],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[4])
#fit a linear model
smix_2020_coumaric_acid_calibration = lm(smix_2020_coumaric_acid_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_coumaric_acid_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_coumaric_acid_calibration, which = 3)


#Ferulic Acid (rt of 35.98 from external standard run)
smix_2020_ferulic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[5], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[5], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[5],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[5])
#fit a linear model
smix_2020_ferulic_acid_calibration = lm(smix_2020_ferulic_acid_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_ferulic_acid_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_ferulic_acid_calibration, which = 3)


#Verbascoside (rt of 41.79 from external standard run)
smix_2020_verbascoside_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[2], 
                                    standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[2], 
                                    standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[2],
                                    standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[2])
#fit a linear model
smix_2020_verbascoside_calibration = lm(smix_2020_verbascoside_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_verbascoside_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_verbascoside_calibration, which = 3)


#Oleuropein (rt of 45.95 from external standard run)
smix_2020_oleuropein_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[3], 
                                        standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[3], 
                                        standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[3],
                                        standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[3])
#fit a linear model
smix_2020_oleuropein_calibration = lm(smix_2020_oleuropein_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_oleuropein_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_oleuropein_calibration, which = 3)


#Pinoresinol (rt of 46.92 from external standard run)
smix_2020_pinoresinol_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[4], 
                                      standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[4], 
                                      standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[4],
                                      standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[4])
#fit a linear model
smix_2020_pinoresinol_calibration = lm(smix_2020_pinoresinol_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_pinoresinol_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_pinoresinol_calibration, which = 3)


#Ligustrocide (rt of 48.43 from external standard run)
smix_2020_ligustrocide_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[5], 
                                       standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[5], 
                                       standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[5],
                                       standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[5])
#fit a linear model
smix_2020_ligustrocide_calibration = lm(smix_2020_ligustrocide_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_ligustrocide_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_ligustrocide_calibration, which = 3)


#Luteolin (rt of 50.23 from external standard run)
smix_2020_luteolin_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[6], 
                                        standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[6], 
                                        standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[6],
                                        standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[6])
#fit a linear model
smix_2020_luteolin_calibration = lm(smix_2020_luteolin_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_luteolin_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_luteolin_calibration, which = 3)


#Apagenin (rt of 52.57 from external standard run)
smix_2020_apagenin_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[7], 
                                    standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[7], 
                                    standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[7],
                                    standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[7])
#fit a linear model
smix_2020_apagenin_calibration = lm(smix_2020_apagenin_absorbance ~ known_conc)
#plot lm residuals vs fitted to check for non-linearity
plot(smix_2020_apagenin_calibration, which = 1)
#plot lm scale-location to check for heteroscedasticity
plot(smix_2020_apagenin_calibration, which = 3)


#compile linear models into vector
standard_calibration_model_list = list(smix_2020_gallic_acid_calibration, 
  smix_2020_chlorogenic_acid_calibration,
  smix_2020_vanillic_acid_calibration,
  smix_2020_eleuth_b_calibration,
  smix_2020_coumaric_acid_calibration,
  smix_2020_ferulic_acid_calibration,
  smix_2020_verbascoside_calibration,
  smix_2020_oleuropein_calibration,
  smix_2020_pinoresinol_calibration,
  smix_2020_ligustrocide_calibration,
  smix_2020_luteolin_calibration,
  smix_2020_apagenin_calibration)



#modify algorithmic detected table method to calculate unknown concentrations
algorithmic_detection_and_concentration_calculation = function(peak_list, stds_data_table, stds_calibration_models) {
  #Create an empty nrow(peak_tab) x 12 numeric data frame to store algorithmic detection data
  alg_det_pk_list_df = data.frame(matrix(NA, length(peak_list), 12))
  rownames(alg_det_pk_list_df) = names(peak_list)
  colnames(alg_det_pk_list_df) = stds_data_table[1, ]
  
  #Iterate over external standard data with walk
  for(col_index in seq_len(ncol(stds_data_table))) {
    #Use col_index to access individual column of the standards data table 
    std_data_table_col = stds_data_table[, col_index]
    
    #store data for external standard from selected column
    std_name = std_data_table_col[1]
    std_peak_rt = as.numeric(std_data_table_col[2])
    std_peak_sd = as.numeric(std_data_table_col[3])
    std_peak_spectra = as.numeric(std_data_table_col[4])
    
    #for current external standard, iterate over all subset runs
    for(sample_index in seq_along(peak_list)) {
      #Use name and positional indexing on sample_index from walk function to 
      #access peak data of individual runs 
      sample_peak_data = peak_list[[sample_index]]$`280`
      
      sample_name = sample_peak_data$sample[1]
      
      #Pull current run RTs as a vector 
      sample_peak_rts = sample_peak_data$rt
      
      #Use multiple RT search to get vector of RTs from sample that fall within
      #sd range, then use logical indexing to select RT that varies least from
      #the RT of the current external standard
      sample_det_rt_list = get_rt_in_sd_list(sample_peak_rts, std_peak_rt, std_peak_sd)
      sample_det_rt = sample_det_rt_list[which.min(abs(sample_det_rt_list - std_peak_rt))]
 
      if(length(sample_det_rt) == 1) {
        #use a mix of numerical and logical indexing to get rts corresponding absorbance
        unknown_absorbance = as.numeric(peak_list[[sample_index]]$`280`[[9]][which(peak_list[[sample_index]]$`280`[[3]] == sample_det_rt[1])])
        print(unknown_absorbance)
        predicted_conc = inverse.predict(stds_calibration_models[[col_index]], unknown_absorbance)
        alg_det_pk_list_df[which(rownames(alg_det_pk_list_df) == sample_name), 
                           which(colnames(alg_det_pk_list_df) == std_name)] = predicted_conc$Prediction
      }
    }
  }
  return(alg_det_pk_list_df)
}

algorithmic_detection_and_concentration_calculation(doe_pre_corrected_75_pks_egh, standards_data_table, standard_calibration_model_list)
algorithmic_detection_and_concentration_calculation(fre_pre_corrected_75_pks_egh, standards_data_table, standard_calibration_model_list)

standard_mix_1 = hplc_list_standards[c(12:15)]
standard_mix_2 = hplc_list_standards[c(16:19)]
standard_mix_3 = hplc_list_standards[c(35:38)]
standard_mix_check = hplc_list_standards[c(20,27:29,31:34)]
standard_mix_check_standard_1 = hplc_list_standards[c(39:47)]
standard_mix_check_standard_2 = hplc_list_standards[c(48,52:60)]
apagenin = hplc_list_standards[c(49:51)]


pre_standard_mix_1 =  preprocess(standard_mix_1, 
                                 dim1=seq(0,54,.01),
                                 dim2=seq(190,400,1),
                                 remove.time.baseline = TRUE,
                                 cl=TRUE)

corrected_rt_standard_mix_1 = correct_rt(pre_standard_mix_1, what = c("corrected.values"), lambdas=c("280"), 
                                         models = NULL, reference = "best", alg = c("vpdtw"), init.coef = c(0, 1, 0), 
                                         n.traces = NULL, n.zeros = 0, plot_it = FALSE, penalty = 1, 
                                         maxshift = 75, verbose = FALSE, progress_bar = FALSE)

plot_chroms(corrected_rt_standard_mix_1, lambdas = '280')

#splitting up peak detection in subset by treatment for better peak detection
standard_mix_1_peaks_split_1 = get_peaks(corrected_rt_standard_mix_1[1], lambdas = c('280'), fit = "egh", amp_thresh = 15000) 
print(standard_mix_1_peaks_split_1)

standard_mix_1_peaks_split_2 = get_peaks(corrected_rt_standard_mix_1[2], lambdas = c('280'), fit = "egh", amp_thresh = 7000) 
print(standard_mix_1_peaks_split_2)

standard_mix_1_peaks_split_3 = get_peaks(corrected_rt_standard_mix_1[3], lambdas = c('280'), fit = "egh", amp_thresh = 3500) 
#Extra fiddling to standardize data
standard_mix_1_peaks_split_3$`ASCIIData_Run#1_Comma013`$'280' = standard_mix_1_peaks_split_3$`ASCIIData_Run#1_Comma013`$'280'[-7,]
rownames(standard_mix_1_peaks_split_3$`ASCIIData_Run#1_Comma013`$'280')[7] = "7"
print(standard_mix_1_peaks_split_3)

standard_mix_1_peaks_split_4 = get_peaks(corrected_rt_standard_mix_1[4], lambdas = c('280'), fit = "egh", amp_thresh = 1700) 
#Extra fiddling to standardize data
standard_mix_1_peaks_split_4$`ASCIIData_Run#1_Comma014`$'280' = standard_mix_1_peaks_split_4$`ASCIIData_Run#1_Comma014`$'280'[-6,]
rownames(standard_mix_1_peaks_split_4$`ASCIIData_Run#1_Comma014`$'280')[6] = "6"
standard_mix_1_peaks_split_4$`ASCIIData_Run#1_Comma014`$'280' = standard_mix_1_peaks_split_4$`ASCIIData_Run#1_Comma014`$'280'[-7,]
rownames(standard_mix_1_peaks_split_4$`ASCIIData_Run#1_Comma014`$'280')[7] = "7"
print(standard_mix_1_peaks_split_4)

#Concatenating split peak lists to make a peak table
standard_mix_1_peaks = c(standard_mix_1_peaks_split_1, standard_mix_1_peaks_split_2, standard_mix_1_peaks_split_3, standard_mix_1_peaks_split_4)
print(standard_mix_1_peaks)
for (attr_name in names(attributes(standard_mix_1_peaks_split_1))) {
  if(attr_name != "names"){
    attr_value <- attributes(standard_mix_1_peaks_split_1)[[attr_name]]
    attr(standard_mix_1_peaks, attr_name) <- attr_value
  }
}
attr(standard_mix_1_peaks, "chrom_list") = "pre_standard_mix_1"
print(standard_mix_1_peaks)

standard_mix_1_peak_table = get_peaktable(standard_mix_1_peaks, response = 'area')
print(standard_mix_1_peak_table$pk_meta)




#Doing the same for standard mix 2
pre_standard_mix_2 =  preprocess(standard_mix_2, 
                                 dim1=seq(0,54,.01),
                                 dim2=seq(190,400,1),
                                 remove.time.baseline = TRUE,
                                 cl=TRUE)

plot_chroms(pre_standard_mix_2, lambdas = c('280'), show_legend = FALSE)

corrected_rt_standard_mix_2 = correct_rt(pre_standard_mix_2, what = c("corrected.values"), lambdas=c("280"), 
                                         models = NULL, reference = "best", alg = c("vpdtw"), init.coef = c(0, 1, 0), 
                                         n.traces = NULL, n.zeros = 0, plot_it = FALSE, penalty = 1, 
                                         maxshift = 75, verbose = FALSE, progress_bar = FALSE)

plot_chroms(corrected_rt_standard_mix_2, lambdas = '280')

#splitting up peak detection in subset by treatment for better peak detection
standard_mix_2_peaks_split_1 = get_peaks(corrected_rt_standard_mix_2[1], lambdas = c('280'), fit = "egh", amp_thresh = 800000) 
print(standard_mix_2_peaks_split_1)

standard_mix_2_peaks_split_2 = get_peaks(corrected_rt_standard_mix_2[2], lambdas = c('280'), fit = "egh", amp_thresh = 150000) 
print(standard_mix_2_peaks_split_2)

standard_mix_2_peaks_split_3 = get_peaks(corrected_rt_standard_mix_2[3], lambdas = c('280'), fit = "egh", amp_thresh = 8000)
print(standard_mix_2_peaks_split_3)

standard_mix_2_peaks_split_4 = get_peaks(corrected_rt_standard_mix_2[4], lambdas = c('280'), fit = "egh", amp_thresh = 5000) 
print(standard_mix_2_peaks_split_4)

#Concatenating split peak lists to make a peak table
standard_mix_2_peaks = c(standard_mix_2_peaks_split_1, standard_mix_2_peaks_split_2, standard_mix_2_peaks_split_3, standard_mix_2_peaks_split_4)
print(standard_mix_2_peaks)
for (attr_name in names(attributes(standard_mix_2_peaks_split_1))) {
  if(attr_name != "names"){
    attr_value <- attributes(standard_mix_2_peaks_split_1)[[attr_name]]
    attr(standard_mix_2_peaks, attr_name) <- attr_value
  }
}
attr(standard_mix_2_peaks, "chrom_list") = "pre_standard_mix_2"
print(standard_mix_2_peaks)

standard_mix_2_peak_table = get_peaktable(standard_mix_2_peaks, response = 'area')
print(standard_mix_2_peak_table$pk_meta)




pre_standard_mix_3 =  preprocess(standard_mix_3, 
                                 dim1=seq(0,54,.01),
                                 dim2=seq(274,284,2),
                                 remove.time.baseline = TRUE,
                                 cl=TRUE)

corrected_rt_standard_mix_3 = correct_rt(pre_standard_mix_3, what = c("corrected.values"), lambdas=c("280"), 
                                         models = NULL, reference = "best", alg = c("vpdtw"), init.coef = c(0, 1, 0), 
                                         n.traces = NULL, n.zeros = 0, scale = FALSE, plot_it = FALSE, penalty = 1, 
                                         maxshift = 75, verbose = FALSE, progress_bar = FALSE)

plot_chroms(corrected_rt_standard_mix_3, lambdas = '280')

standard_mix_3_peaks = get_peaks(pre_standard_mix_3, lambdas = c('280'), fit = "egh", amp_thresh = 10000) 
standard_mix_3_peak_table = get_peaktable(standard_mix_3_peaks, response = 'area')
print(standard_mix_3_peak_table$pk_meta)








pre_standard_mix_check = preprocess(standard_mix_check, 
                                    dim1=seq(0,54,.05),
                                    dim2=seq(274,284,2),
                                    remove.time.baseline = TRUE,
                                    cl=TRUE)
plot_chroms(pre_standard_mix_check, lambdas = '280', show_legend = FALSE)

pre_standard_mix_check_standard_1 = preprocess(standard_mix_check_standard_1, 
                                    dim1=seq(0,54,.05),
                                    dim2=seq(274,284,2),
                                    remove.time.baseline = TRUE,
                                    cl=TRUE)
plot_chroms(pre_standard_mix_check_standard_1, lambdas = '280', show_legend = FALSE)

pre_standard_mix_check_standard_2 = preprocess(standard_mix_check_standard_2, 
                                    dim1=seq(0,54,.05),
                                    dim2=seq(274,284,2),
                                    remove.time.baseline = TRUE,
                                    cl=TRUE)
plot_chroms(pre_standard_mix_check_standard_2, lambdas = '280', show_legend = FALSE)



get_subset_time_range = function(pre_subset) {
  pre_subset_check_time_range = data.frame(Run_Name = character(0), Run_DateTime = character(0))
  
  for(i in 1:length(pre_subset)) {
    new_row = data.frame(Run_Name = attr(pre_subset, "names")[i], Run_DateTime = attr(pre_subset[[i]], "run_datetime"))
    pre_subset_check_time_range = rbind(pre_subset_check_time_range, new_row)
    
  }
  return(pre_subset_check_time_range)
}

get_subset_time_range(standard_mix_check)
get_subset_time_range(standard_mix_check_standard_1)
get_subset_time_range(standard_mix_check_standard_2)


pre_apagenin =  preprocess(apagenin, 
                                 dim1=seq(0,54,.01),
                                 dim2=seq(274,284,2),
                                 remove.time.baseline = TRUE,
                                 cl=TRUE)

plot_chroms(pre_apagenin, lambdas = '280')

apagenin_peaks = get_peaks(pre_apagenin, lambdas = c('280'), fit = "egh", amp_thresh = 5000) 
apagenin_peak_table = get_peaktable(apagenin_peaks, response = 'area')

standard_mix21_1 = standards21_pre_standard_mix[c(1,2)]
standard_mix21_2 = standards21_pre_standard_mix[c(3,4)]
standard_mix21_3 = standards21_pre_standard_mix[5]
standard_mix21_check_standard = standards21_pre_standard_mix[6:26]
standard_mix21_check_standard_1 = standard_mix21_check_standard[1:6]
standard_mix21_check_standard_2 = standard_mix21_check_standard[7:21]


pre_standard_mix21_1 =  preprocess(standard_mix21_1, 
                                 dim1=seq(0,45,.01),
                                 dim2=seq(274,284,2),
                                 remove.time.baseline = TRUE,
                                 cl=TRUE)
plot_chroms(pre_standard_mix21_1, lambdas = '280')

standard_mix21_1_peaks = get_peaks(pre_standard_mix21_1, lambdas = c('280'), fit = "egh", amp_thresh = 3700) 
standard_mix21_1_peak_table = get_peaktable(standard_mix21_1_peaks, response = 'area')



pre_standard_mix21_2 =  preprocess(standard_mix21_2, 
                                   dim1=seq(0,45,.01),
                                   dim2=seq(274,284,2),
                                   remove.time.baseline = TRUE,
                                   cl=TRUE)
plot_chroms(pre_standard_mix21_2, lambdas = '280')

standard_mix21_2_peaks = get_peaks(pre_standard_mix21_2, lambdas = c('280'), fit = "egh", amp_thresh = 2000) 
standard_mix21_2_peak_table = get_peaktable(standard_mix21_2_peaks, response = 'area')


pre_standard_mix21_3 =  preprocess(standard_mix21_3, 
                                   dim1=seq(0,45,.01),
                                   dim2=seq(274,284,2),
                                   remove.time.baseline = TRUE,
                                   cl=TRUE)
plot_chroms(pre_standard_mix21_3, lambdas = '280')

standard_mix21_3_peaks = get_peaks(pre_standard_mix21_3, lambdas = c('280'), fit = "egh", amp_thresh = 3600) 
standard_mix21_3_peak_table = get_peaktable(standard_mix21_3_peaks, response = 'area')


pre_standard_mix21_check_standard_1 =  preprocess(standard_mix21_check_standard_1, 
                                   dim1=seq(0,45,.01),
                                   dim2=seq(274,284,2),
                                   remove.time.baseline = TRUE,
                                   cl=TRUE)
plot_chroms(pre_standard_mix21_check_standard_1, lambdas = '280', show_legend = FALSE)

standard_mix21_check_standard_1_peaks = get_peaks(pre_standard_mix21_check_standard_1, lambdas = c('280'), fit = "egh", amp_thresh = 2000) 
standard_mix21_check_standard_1_peak_table = get_peaktable(standard_mix21_check_standard_1_peaks, response = 'area')


pre_standard_mix21_check_standard_2 =  preprocess(standard_mix21_check_standard_2, 
                                                  dim1=seq(0,45,.01),
                                                  dim2=seq(274,284,2),
                                                  remove.time.baseline = TRUE,
                                                  cl=TRUE)
plot_chroms(pre_standard_mix21_check_standard_2, lambdas = '280', show_legend = FALSE)

standard_mix21_check_standard_2_peaks = get_peaks(pre_standard_mix21_check_standard_2, lambdas = c('280'), fit = "egh", amp_thresh = 2000) 
standard_mix21_check_standard_2_peak_table = get_peaktable(standard_mix21_check_standard_2_peaks, response = 'area')

stopCluster(cl)

