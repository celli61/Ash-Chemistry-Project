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

#Chlorogenic Acid (rt of 20.17 from external standard run)
smix_2020_chlorogenic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[2], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[2], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[2],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[2])

#Vanillic Acid (rt of 22.29 from external standard run)
smix_2020_vanillic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[3], 
                                          standard_mix_2020_2_peaks[[2]]$'280'$height[3], 
                                          standard_mix_2020_2_peaks[[3]]$'280'$height[3],
                                          standard_mix_2020_2_peaks[[4]]$'280'$height[3])

#Eluetheroside B (rt of 23.87 from external standard run)
smix_2020_eleuth_b_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[1], 
                                    standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[1], 
                                    standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[1],
                                    standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[1])

#Coumaric Acid (rt of 32.16 from external standard run)
smix_2020_coumaric_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[4], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[4], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[4],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[4])

#Ferulic Acid (rt of 35.98 from external standard run)
smix_2020_ferulic_acid_absorbance = c(standard_mix_2020_2_peaks[[1]]$'280'$height[5], 
                                       standard_mix_2020_2_peaks[[2]]$'280'$height[5], 
                                       standard_mix_2020_2_peaks[[3]]$'280'$height[5],
                                       standard_mix_2020_2_peaks[[4]]$'280'$height[5])

#Verbascoside (rt of 41.79 from external standard run)
smix_2020_verbascoside_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[2], 
                                    standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[2], 
                                    standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[2],
                                    standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[2])

#Oleuropein (rt of 45.95 from external standard run)
smix_2020_oleuropein_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[3], 
                                        standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[3], 
                                        standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[3],
                                        standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[3])

#Pinoresinol (rt of 46.92 from external standard run)
smix_2020_pinoresinol_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[4], 
                                      standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[4], 
                                      standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[4],
                                      standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[4])

#Ligustrocide (rt of 48.43 from external standard run)
smix_2020_ligustrocide_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[5], 
                                       standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[5], 
                                       standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[5],
                                       standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[5])

#Luteolin (rt of 50.23 from external standard run)
smix_2020_luteolin_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[6], 
                                        standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[6], 
                                        standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[6],
                                        standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[6])

#Apagenin (rt of 52.57 from external standard run)
smix_2020_apagenin_absorbance = c(standard_mix_2020_1_1mg_peaks[[1]]$'280'$height[7], 
                                    standard_mix_2020_1_p5mg_peaks[[1]]$'280'$height[7], 
                                    standard_mix_2020_1_p25mg_peaks[[1]]$'280'$height[7],
                                    standard_mix_2020_1_p125mg_peaks[[1]]$'280'$height[7])



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
                                 cl = cl)

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
                                 cl = cl)

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
                                 cl = cl)

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
                                    cl = cl)
plot_chroms(pre_standard_mix_check, lambdas = '280', show_legend = FALSE)

pre_standard_mix_check_standard_1 = preprocess(standard_mix_check_standard_1, 
                                    dim1=seq(0,54,.05),
                                    dim2=seq(274,284,2),
                                    remove.time.baseline = TRUE,
                                    cl = cl)
plot_chroms(pre_standard_mix_check_standard_1, lambdas = '280', show_legend = FALSE)

pre_standard_mix_check_standard_2 = preprocess(standard_mix_check_standard_2, 
                                    dim1=seq(0,54,.05),
                                    dim2=seq(274,284,2),
                                    remove.time.baseline = TRUE,
                                    cl = cl)
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
                                 cl = cl)

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
                                 cl = cl)
plot_chroms(pre_standard_mix21_1, lambdas = '280')

standard_mix21_1_peaks = get_peaks(pre_standard_mix21_1, lambdas = c('280'), fit = "egh", amp_thresh = 3700) 
standard_mix21_1_peak_table = get_peaktable(standard_mix21_1_peaks, response = 'area')



pre_standard_mix21_2 =  preprocess(standard_mix21_2, 
                                   dim1=seq(0,45,.01),
                                   dim2=seq(274,284,2),
                                   remove.time.baseline = TRUE,
                                   cl = cl)
plot_chroms(pre_standard_mix21_2, lambdas = '280')

standard_mix21_2_peaks = get_peaks(pre_standard_mix21_2, lambdas = c('280'), fit = "egh", amp_thresh = 2000) 
standard_mix21_2_peak_table = get_peaktable(standard_mix21_2_peaks, response = 'area')


pre_standard_mix21_3 =  preprocess(standard_mix21_3, 
                                   dim1=seq(0,45,.01),
                                   dim2=seq(274,284,2),
                                   remove.time.baseline = TRUE,
                                   cl = cl)
plot_chroms(pre_standard_mix21_3, lambdas = '280')

standard_mix21_3_peaks = get_peaks(pre_standard_mix21_3, lambdas = c('280'), fit = "egh", amp_thresh = 3600) 



pre_standard_mix21_check_standard_1 =  preprocess(standard_mix21_check_standard_1, 
                                   dim1=seq(0,45,.01),
                                   dim2=seq(274,284,2),
                                   remove.time.baseline = TRUE,
                                   cl = cl)
plot_chroms(pre_standard_mix21_check_standard_1, lambdas = '280', show_legend = FALSE)

standard_mix21_check_standard_1_peaks = get_peaks(pre_standard_mix21_check_standard_1, lambdas = c('280'), fit = "egh", amp_thresh = 2000) 
standard_mix21_check_standard_1_peak_table = get_peaktable(standard_mix21_check_standard_1_peaks, response = 'area')


pre_standard_mix21_check_standard_2 =  preprocess(standard_mix21_check_standard_2, 
                                                  dim1=seq(0,45,.01),
                                                  dim2=seq(274,284,2),
                                                  remove.time.baseline = TRUE,
                                                  cl = cl)
plot_chroms(pre_standard_mix21_check_standard_2, lambdas = '280', show_legend = FALSE)

standard_mix21_check_standard_2_peaks = get_peaks(pre_standard_mix21_check_standard_2, lambdas = c('280'), fit = "egh", amp_thresh = 2000) 
standard_mix21_check_standard_2_peak_table = get_peaktable(standard_mix21_check_standard_2_peaks, response = 'area')

stopCluster(cl)




###---2020(Standard Mix)---###

#Gallic Acid
#SD
standard_mix_gallic_acid_sd = ((standard_mix_2_peak_table$pk_meta[3, 1]-standard_mix_2_peak_table$pk_meta[4, 1])+
                                 (standard_mix_2_peak_table$pk_meta[5, 1]-standard_mix_2_peak_table$pk_meta[3, 1])/2)+standard_mix_2_peak_table$pk_meta[6,1]

#RT
standard_mix_gallic_acid_peak_rt = standard_mix_2_peak_table$pk_meta[3, 1]


#Chlorogenic Acid
#SD
standard_mix_chlorogenic_acid_sd = ((standard_mix_2_peak_table$pk_meta[3, 2]-standard_mix_2_peak_table$pk_meta[4, 2])+
                                      (standard_mix_2_peak_table$pk_meta[5, 2]-standard_mix_2_peak_table$pk_meta[3, 2])/2)+standard_mix_2_peak_table$pk_meta[6, 2]

#RT
standard_mix_chlorogenic_acid_peak_rt = standard_mix_2_peak_table$pk_meta[3, 2]


#Vanillic Acid
#SD
standard_mix_vanillic_acid_sd = ((standard_mix_2_peak_table$pk_meta[3, 3]-standard_mix_2_peak_table$pk_meta[4, 3])+
                                   (standard_mix_2_peak_table$pk_meta[5, 3]-standard_mix_2_peak_table$pk_meta[3, 3])/2)+standard_mix_2_peak_table$pk_meta[6, 3]

#RT
standard_mix_vanillic_acid_peak_rt = standard_mix_2_peak_table$pk_meta[3, 3]


#Eleutheroside B
#SD
standard_mix_eleutheroside_b_sd = ((standard_mix_1_peak_table$pk_meta[3, 1]-standard_mix_1_peak_table$pk_meta[4, 1])+
                                     (standard_mix_1_peak_table$pk_meta[5, 1]-standard_mix_1_peak_table$pk_meta[3, 1])/2)+standard_mix_1_peak_table$pk_meta[6, 1]

#RT
standard_mix_eleutheroside_b_peak_rt = standard_mix_1_peak_table$pk_meta[3, 1]


#Coumaric Acid
#SD
standard_mix_coumaric_acid_sd = ((standard_mix_2_peak_table$pk_meta[3, 4]-standard_mix_2_peak_table$pk_meta[4, 4])+
                                   (standard_mix_2_peak_table$pk_meta[5, 4]-standard_mix_2_peak_table$pk_meta[3, 4])/2)+standard_mix_2_peak_table$pk_meta[6, 4]

#RT
standard_mix_coumaric_acid_peak_rt = standard_mix_2_peak_table$pk_meta[3, 4]


#Ferulic Acid
#SD
standard_mix_ferulic_acid_sd = ((standard_mix_2_peak_table$pk_meta[3, 5]-standard_mix_2_peak_table$pk_meta[4, 5])+
                                  (standard_mix_2_peak_table$pk_meta[5, 5]-standard_mix_2_peak_table$pk_meta[3, 5])/2)+standard_mix_2_peak_table$pk_meta[6, 5]

#RT
standard_mix_ferulic_acid_peak_rt = standard_mix_2_peak_table$pk_meta[3, 5]


#Verbascoside
#SD
standard_mix_verbascoside_sd = ((standard_mix_1_peak_table$pk_meta[3, 2]-standard_mix_1_peak_table$pk_meta[4, 2])+
                                  (standard_mix_1_peak_table$pk_meta[5, 2]-standard_mix_1_peak_table$pk_meta[3, 2])/2)+standard_mix_1_peak_table$pk_meta[6, 2]

#RT
standard_mix_verbascoside_peak_rt = standard_mix_1_peak_table$pk_meta[3, 2]


#Oleuropein
#SD
standard_mix_oleuropein_sd = ((standard_mix_1_peak_table$pk_meta[3, 3]-standard_mix_1_peak_table$pk_meta[4, 3])+
                                (standard_mix_1_peak_table$pk_meta[5, 3]-standard_mix_1_peak_table$pk_meta[3, 3])/2)+standard_mix_1_peak_table$pk_meta[6, 3]

#RT
standard_mix_oleuropein_peak_rt = standard_mix_1_peak_table$pk_meta[3, 3]


#Pinoresinol
#SD
standard_mix_pinoresinol_sd = ((standard_mix_1_peak_table$pk_meta[3, 4]-standard_mix_1_peak_table$pk_meta[4, 4])+
                                 (standard_mix_1_peak_table$pk_meta[5, 4]-standard_mix_1_peak_table$pk_meta[3, 4])/2)+standard_mix_1_peak_table$pk_meta[6, 4]

#RT
standard_mix_pinoresinol_peak_rt = standard_mix_1_peak_table$pk_meta[3, 4]


#Ligustrocide
#SD
standard_mix_ligustrocide_sd = ((standard_mix_1_peak_table$pk_meta[3, 5]-standard_mix_1_peak_table$pk_meta[4, 5])+
                                  (standard_mix_1_peak_table$pk_meta[5, 5]-standard_mix_1_peak_table$pk_meta[3, 5])/2)+standard_mix_1_peak_table$pk_meta[6, 5]

#RT
standard_mix_ligustrocide_peak_rt = standard_mix_1_peak_table$pk_meta[3, 5]


#Luteolin
#SD
standard_mix_luteolin_sd = ((standard_mix_1_peak_table$pk_meta[3, 6]-standard_mix_1_peak_table$pk_meta[4, 6])+
                              (standard_mix_1_peak_table$pk_meta[5, 6]-standard_mix_1_peak_table$pk_meta[3, 6])/2)+standard_mix_1_peak_table$pk_meta[6, 6]

#RT
standard_mix_luteolin_peak_rt = standard_mix_1_peak_table$pk_meta[3, 6]


#Apagenin
#SD
standard_mix_apagenin_sd = ((standard_mix_1_peak_table$pk_meta[3, 7]-standard_mix_1_peak_table$pk_meta[4, 7])+
                              (standard_mix_1_peak_table$pk_meta[5, 7]-standard_mix_1_peak_table$pk_meta[3, 7])/2)+standard_mix_1_peak_table$pk_meta[6, 7]

#RT
standard_mix_apagenin_peak_rt = standard_mix_1_peak_table$pk_meta[3, 7]


#Creating the standard data table for iterative access

#storing like data in lists of the identical ordering
standard_name_list_standard_mix = list("Gallic Acid", "Chlorogenic Acid", "Vanillic Acid", "Eleutheroside B", "Coumaric Acid", 
                                       "Ferulic Acid", "Verbascoside", "Oleuropein", "Pinoresinol", "Ligustrocide", 
                                       "Luteolin", "Apagenin")
#creating a data frame from the list
standard_name_standard_mix = data.frame(standard_name_list_standard_mix)
#setting the column names to empty strings to prevent errors with rbind
colnames(standard_name_standard_mix) = rep("", 12)
#print to see if it looks right
print(standard_name_standard_mix)


#storing like data in lists of identical ordering
rt_list_standard_mix = list(gallic_acid_peak_rt, chlorogenic_acid_peak_rt, vanillic_acid_peak_rt, 
                            eleutheroside_b_peak_rt, coumaric_acid_peak_rt, ferulic_acid_peak_rt, 
                            verbascoside_peak_rt, oleuropein_peak_rt, pinoresinol_peak_rt, ligustrocide_peak_rt, 
                            luteolin_peak_rt, apagenin_peak_rt)
#creating a data frame from the list
standard_peak_rt_standard_mix = data.frame(rt_list_standard_mix)
#setting the column names to empty strings to prevent errors with rbind
colnames(standard_peak_rt_standard_mix) = rep("", 12)
#print to see if it looks right
print(standard_peak_rt_standard_mix)

#storing like data in lists of identical ordering
sd_list_standard_mix = list(gallic_acid_sd, chlorogenic_acid_sd, vanillic_acid_sd, eleutheroside_b_sd, 
                            coumaric_acid_sd, ferulic_acid_sd, verbascoside_sd, oleuropein_sd, pinoresinol_sd, 
                            ligustrocide_sd, luteolin_sd, apagenin_sd)
#creating a data frame from the list
standard_peak_sd_standard_mix = data.frame(sd_list_standard_mix)
#using lapply to round all data to 3 sig figs
standard_peak_sd_standard_mix = as.data.frame(lapply(standard_peak_sd_standard_mix, function(x) round(x, 3)))
#setting column names to empty strings to prevent errors with rbind
colnames(standard_peak_sd_standard_mix) = rep("", 12)
#print to see if it looks right
print(standard_peak_sd_standard_mix)


#Binding all individual data frames into a table with rbind
standards_data_table_standard_mix = rbind(standard_name_standard_mix, standard_peak_rt_standard_mix, standard_peak_sd_standard_mix)
#printing to show data structure
print(standards_data_table_standard_mix)


