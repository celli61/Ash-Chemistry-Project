##############################################################
# Analysis and construction of data into relevant structures #
##############################################################


###---creating the linear search method---###
# linear_search(list, key, sd)
# 
# A simple iterative linear search method that iterates through a list searching for a specified element, 
# optionally within a standard deviation, and returns the index of the element in the list if found as a numeric.
# If the element is not found in the specified list, the method returns -1. 
# 
# @param list - The list object to search through.
# @param key - The target value to search for.
# @param sd - [OPTIONAL] The standard deviation from the key to define a range of values an element of list can fall in 
#             to be considered an instance of the key value. defaults to 0 if unspecified.
#
# @return - A numeric value indicating the index of the key in the list, or a -1 indicating the value is not in the list.
#
linear_search = function(list, key, sd = 0){ 
  #Loop through each element in the list
  for(i in 1:length(list)){
    #If the object at index i in list falls in the range of the value of key +/- the standard deviation sd
    if(list[i] >= key-sd & list[i] <= key+sd){ 
      return(i) #Return the current index, i
    }
  }
  return(-1) #Return -1 if the key is not found within the defined range
}

###---creating the retention time within standard deviation range getter method---###
#2/27/25 - changed problematic behavior that resulted in premature return when 
#          encountering a value outside the range.
get_rt_in_sd_list = function(rt_list, rt_key, rt_sd){
  rt_in_sd_list = c() #instantiate an empty vector to contain RTs
  
  for(rt in rt_list){ #loop through RTs in the list
    if(rt >= rt_key - rt_sd && rt <= rt_key + rt_sd){ #if the current RT is within range
      rt_in_sd_list = c(rt_in_sd_list, rt) #append RT to list
    }
  }
  return(rt_in_sd_list) #return the list
}

doe_rt_list = get_rt_in_sd_list(doe_pk_tab$pk_meta[3,], 20.17, 1.14)
get_rt_in_sd_list(doe_pk_tab$pk_meta[3,], 5.42, 0.735)



###---Creating a max absorbance spectra table method---###
# get_max_absorbance_spectra_table(pk_tab)
#
# A method that searches through a peak list's reference spectra data and creates a table displaying each peak, and it's 
# corresponding maximum absorbance spectrum.
#
# @param pk_tab - A chromatographR peak table object with or without attached reference spectra data.
#
# @return - A data frame with each peak from the table and it's corresponding maximum absorbance spectrum.
#
get_max_absorbance_spectra_table = function(pk_tab){
  #Check if the peak table does not have ref_spectra attached
  if(is.na(pk_tab$ref_spectra)){
    #Attach ref_spectra to the peak table
    pk_tab = attach_ref_spectra(pk_tab, ref = "max.int")
    
  }
  #Initialize max_absorbance_spectra_df with all peaks in the Peaks column and and empty Max_Spectra column of the same length
  max_absorbance_spectra_df = data.frame(Peaks = colnames(pk_tab$ref_spectra), 
                                         Max_Spectra = character(length(colnames(pk_tab$ref_spectra))))
  
  #Iterate through the peaks
  for(i in 1:length(max_absorbance_spectra_df$Peaks)){
    #Call the linear_search method to find the index of the current peak's max spectrum from ref_spectra
    max_spectrum_index = linear_search(pk_tab$ref_spectra[,i], 1)
    
    #Check that the current peak's max spectrum is present in the given spectral range
    if(max_spectrum_index != -1){
      #Set the Max_Spectra for the current peak to the spectrum indexed by max_spectrum_index
      max_absorbance_spectra_df[i,2] = names(pk_tab$ref_spectra[,i])[max_spectrum_index]
      
    } else {
      #Set the Max_Spectra for the current peak to NA
      max_absorbance_spectra_df[i,2] = "NA"
    }
  }
  return(max_absorbance_spectra_df)
}


#Printing the max absorbance spectra table for doe 2020
get_max_absorbance_spectra_table(doe_pk_tab)

#Calling and storing max absorbance spectra for all 2020 standards
gallic_acid_max_spectra = get_max_absorbance_spectra_table(standards_gallic_acid_peak_table)$'Max_Spectra'
chlorogenic_acid_max_spectra = get_max_absorbance_spectra_table(standards_chlorogenic_acid_peak_table)$'Max_Spectra'
vanillic_acid_max_spectra = get_max_absorbance_spectra_table(standards_vanillic_acid_peak_table)$'Max_Spectra'
eleutheroside_b_max_spectra = get_max_absorbance_spectra_table(standards_eleutheroside_b_peak_table)$'Max_Spectra'
coumaric_acid_max_spectra = get_max_absorbance_spectra_table(standards_coumaric_acid_peak_table)$'Max_Spectra'
ferulic_acid_max_spectra = get_max_absorbance_spectra_table(standards_ferulic_acid_peak_table)$'Max_Spectra'
verbascoside_max_spectra = get_max_absorbance_spectra_table(standards_verbascoside_peak_table)$'Max_Spectra'
oleuropein_max_spectra = get_max_absorbance_spectra_table(standards_oleuropein_peak_table)$'Max_Spectra'
pinoresinol_max_spectra = get_max_absorbance_spectra_table(standards_pinoresinol_peak_table)$'Max_Spectra'
ligustrocide_max_spectra = get_max_absorbance_spectra_table(standards_ligustrocide_peak_table)$'Max_Spectra'
luteolin_max_spectra = get_max_absorbance_spectra_table(standards_luteolin_peak_table)$'Max_Spectra'
apagenin_max_spectra = get_max_absorbance_spectra_table(standards_apagenin_peak_table)$'Max_Spectra'


#Binding all max spectra data together in a vector
max_spectra_list = list(gallic_acid_max_spectra, chlorogenic_acid_max_spectra, vanillic_acid_max_spectra, eleutheroside_b_max_spectra, 
                coumaric_acid_max_spectra, ferulic_acid_max_spectra, verbascoside_max_spectra, oleuropein_max_spectra, 
                pinoresinol_max_spectra, ligustrocide_max_spectra, luteolin_max_spectra, apagenin_max_spectra)

#Storing the vector as a data frame
Standard_Max_Spectra = data.frame(max_spectra_list)
#removing column names
colnames(Standard_Max_Spectra) = rep("", 12)

print(Standard_Max_Spectra)

standards_data_table = rbind(standards_data_table, Standard_Max_Spectra)

print(standards_data_table)

#Printing max absorbance for other standards 2019 (in standard mix)
standard_mix_1_ref_spec = get_max_absorbance_spectra_table(standard_mix_1_peak_table)
standard_mix_2_ref_spec = get_max_absorbance_spectra_table(standard_mix_2_peak_table)


#Adding standard mix max spectra to standard mix data table 
standard_mix_ref_spec = c(standard_mix_2_ref_spec$Max_Spectra[1:3], standard_mix_1_ref_spec$Max_Spectra[1], 
                          standard_mix_2_ref_spec$Max_Spectra[4:5], standard_mix_1_ref_spec$Max_Spectra[2:7])

standards_data_table_standard_mix = rbind(standards_data_table_standard_mix, unname(standard_mix_ref_spec))


#Calling and storing max absorbance spectra for all 2021 standards
gallic_acid21_max_spectra = get_max_absorbance_spectra_table(standards21_gallic_acid_peak_table)$'Max_Spectra'
chlorogenic_acid21_max_spectra = get_max_absorbance_spectra_table(standards21_chlorogenic_acid_peak_table)$'Max_Spectra'
vanillic_acid21_max_spectra = get_max_absorbance_spectra_table(standards21_vanillic_acid_peak_table)$'Max_Spectra'
syringin21_max_spectra = get_max_absorbance_spectra_table(standards21_syringin_peak_table)$'Max_Spectra'
coumaric_acid21_max_spectra = get_max_absorbance_spectra_table(standards21_coumaric_acid_peak_table)$'Max_Spectra'
ferulic_acid21_max_spectra = get_max_absorbance_spectra_table(standards21_ferulic_acid_peak_table)$'Max_Spectra'
verbascoside21_max_spectra = get_max_absorbance_spectra_table(standards21_verbascoside_peak_table)$'Max_Spectra'
oleuropein21_max_spectra = get_max_absorbance_spectra_table(standards21_oleuropein_peak_table)$'Max_Spectra'
pinoresinol21_max_spectra = get_max_absorbance_spectra_table(standards21_pinoresinol_peak_table)$'Max_Spectra'
ligustrocide21_max_spectra = get_max_absorbance_spectra_table(standards21_ligustrocide_peak_table)$'Max_Spectra'
luteolin21_max_spectra = get_max_absorbance_spectra_table(standards21_luteolin_peak_table)$'Max_Spectra'
apagenin21_max_spectra = get_max_absorbance_spectra_table(standards21_apagenin_peak_table)$'Max_Spectra'


#Binding all max spectra data together in a vector
max_spectra_list_21 = list(gallic_acid21_max_spectra, chlorogenic_acid21_max_spectra, vanillic_acid21_max_spectra, syringin21_max_spectra, 
                        coumaric_acid21_max_spectra, ferulic_acid21_max_spectra, verbascoside21_max_spectra, oleuropein21_max_spectra, 
                        pinoresinol21_max_spectra, ligustrocide21_max_spectra, luteolin21_max_spectra, apagenin21_max_spectra)

#Storing the vector as a data frame
standard_max_spectra_21 = data.frame(max_spectra_list_21)
#removing column names
colnames(standard_max_spectra_21) = rep("", 12)

print(standard_max_spectra_21)

standards_data_table_21 = rbind(standards_data_table_21, standard_max_spectra_21)

print(standards_data_table_21)






###---creating a standard presence table method---###
# get_standard_presence_table(peak_tab, data_table)
#
# This method generates a Boolean matrix to indicate the presence of each of the 12 standards in all treatments in a 
# subset defined by peak_tab and data_table objects.
#
# @param peak_tab - A chromatographR peak_table object for a subset.
# @param data_table - A data frame containing the peak data for all 12 standards.
# 
# @return -  A Boolean matrix with dimensions: number-of-samples-in-subset x 12, where each element indicates the
#            presence (TRUE) or absence (FALSE) of a standard for each treatment in a subset.
#
# @details
#   - Iterates through each subset in std_pres_tab_list and retrieves the corresponding metadata from subset_meta_list.
#   - Computes the number of true and false entries for each standard using the table function
#   - Constructs a new row for the presence rate table for each standard, and appends it to the data frame
#
get_standard_presence_table <- function(peak_tab, data_table){
  #Create an empty nrow(peak_tab) x 12 numeric data frame to store algorithmic detection data 
  standard_presence_table = data.frame(matrix(FALSE, length(peak_tab$tab$V1), 12))
  rownames(standard_presence_table) = rownames(peak_tab$tab)
  colnames(standard_presence_table) = data_table[1, ]
  
  #Retrieve and store the RTs of the all peaks in the peak table as a vector
  peaks_rts = peak_tab$pk_meta[3,]
  
  #Iterate over each standard in data_table
  for(i in 1:length(standard_presence_table)){ 
    #Retrieve and store the RT and SD of the current standard from the data table as a numeric
    standard_rt = as.numeric(data_table[2,i])
    standard_sd = as.numeric(data_table[3,i])
    
    #Call linear_search on peak_rts for the current standard_rt with the current standard_sd and store the return
    standard_index = linear_search(peaks_rts, standard_rt, standard_sd)
    
    #If a valid index (> -1) is found
    if(standard_index != -1){ 
      #Iterate over every treatment in the subset
      for(j in 1:length(peak_tab$tab$V1)){
        
        #Check if the PDA data at the peak corresponding to standard_index is non-zero in the current treatment
        if(peak_tab$tab[j, standard_index] > 0){ 
          #Update spot [j,i] in standard_presence_table to indicate the detection of the current standard
          standard_presence_table[j,i] = TRUE 
        }
      }
    }
  }
  #Return the matrix to be stored in a variable
  return(standard_presence_table)
}


###---Standard presence tables - 2020 subsets---###
#Doe farm treatments 
doe_standard_presence_table = get_standard_presence_table(doe_pk_tab, standards_data_table)
#print(doe_standard_presence_table)

#French Conservation farm treatments
fre_standard_presence_table = get_standard_presence_table(fre_pk_tab, standards_data_table)
#print(fre_standard_presence_table)

#East Foss farm treatments
eas_standard_presence_table = get_standard_presence_table(eas_pk_tab, standards_data_table)
#print(eas_standard_presence_table)

#Jennings forest treatments
jen_standard_presence_table = get_standard_presence_table(jen_pk_tab, standards_data_table)
#print(jen_standard_presence_table)


###---Standard presence tables - 2021 subsets---###
#Tuttle Swamp treatments 
tut_standard_presence_table = get_standard_presence_table(tut_pk_tab, standards_data_table_21)
#print(tut_standard_presence_table)

#Powder Major's forest treatments
pow_standard_presence_table = get_standard_presence_table(pow_pk_tab, standards_data_table_21)
#print(pow_standard_presence_table)

#Town of Lee forest treatments
lee_standard_presence_table = get_standard_presence_table(lee_pk_tab, standards_data_table_21)
#print(lee_standard_presence_table)

#Little River Conservation Area treatments
lit_standard_presence_table = get_standard_presence_table(lit_pk_tab, standards_data_table_21)
#print(lit_standard_presence_table)

#putting all the standard presence tables in a list
standard_presence_table_list = list(doe_standard_presence_table, fre_standard_presence_table, eas_standard_presence_table, jen_standard_presence_table,
                                    tut_standard_presence_table, pow_standard_presence_table, lee_standard_presence_table, lit_standard_presence_table)

#subsetting 2021 meta data to match other subsets
pow_meta = pow[, c(1, 2, 3, 4, 5, 11, 12, 13)]
tut_meta = tut[, c(1, 2, 3, 4, 5, 11, 12, 13)]
lee_meta = lee[, c(1, 2, 3, 4, 5, 11, 12, 13)]
lit_meta = lit[, c(1, 2, 3, 4, 5, 11, 12, 13)]
#putting each subset's data in a list with index corresponding to the subsets position in std_pres_tab_list
subset_metadata_list = list(doe_meta, fre_meta, eas_meta, jen_meta, tut_meta, pow_meta, lee_meta, lit_meta)


###---creating a standard presence table method---###
# get_standard_presence_table(peak_tab, data_table)
#
# This method generates a Boolean matrix to indicate the presence of each of the 12 standards in all treatments in a 
# subset defined by peak_tab and data_table objects.
#
# @param peak_tab - A chromatographR peak_table object for a subset.
# @param data_table - A data frame containing the peak data for all 12 standards.
# 
# @return -  A Boolean matrix with dimensions: number-of-samples-in-subset x 12, where each element indicates the
#            presence (TRUE) or absence (FALSE) of a standard for each treatment in a subset.
#
# @details
#   - Iterates through each subset in std_pres_tab_list and retrieves the corresponding metadata from subset_meta_list.
#   - Computes the number of true and false entries for each standard using the table function
#   - Constructs a new row for the presence rate table for each standard, and appends it to the data frame
#
get_standard_presence_table_spec_check <- function(peak_tab, data_table, spec_sd = 20){
  #Create an empty nrow(peak_tab) x 12 numeric data frame to store algorithmic detection data 
  standard_presence_table_spec_check = data.frame(matrix(FALSE, length(peak_tab$tab$V1), 12))
  rownames(standard_presence_table_spec_check) = rownames(peak_tab$tab)
  colnames(standard_presence_table_spec_check) = data_table[1, ]
  
  #Create and store the max spectra table for the peak tab
  peaks_max_spectra = get_max_absorbance_spectra_table(peak_tab)
  
  #Retrieve and store the RTs of the all peaks in the peak table as a vector
  peaks_rts = peak_tab$pk_meta[3,]
  
  #Iterate over each standard in data_table
  for(i in 1:length(standard_presence_table_spec_check)){ 
    #Retrieve and store the RT and SD of the current standard from the data table as a numeric
    standard_rt = as.numeric(data_table[2,i])
    standard_sd = as.numeric(data_table[3,i])
    standard_max_spec = as.numeric(data_table[4,i])
    
    #Call linear_search on peak_rts for the current standard_rt with the current standard_sd and store the return
    standard_index = linear_search(peaks_rts, standard_rt, standard_sd)
    
    #If a valid index (> -1) is found
    if(standard_index != -1){ 
      
      #Retrieve and store the max absorbance spectra of the peak
      peak_max_spec = as.numeric(peaks_max_spectra[standard_index,2])
        
        #Iterate over every treatment in the subset
        for(j in 1:length(peak_tab$tab$V1)){
          
          #If the current peak's max absorbance spectra matches the current standard's
          if(abs(peak_max_spec - standard_max_spec) <= spec_sd){
          
            #Check if the PDA data at the peak corresponding to standard_index is non-zero in the current treatment
            if(peak_tab$tab[j, standard_index] > 0){ 
              
              #Update spot [j,i] in standard_presence_table to indicate the detection of the current standard
              standard_presence_table_spec_check[j,i] = TRUE 
            }
          }
        }
    }
  }
  #Return the matrix to be stored in a variable
  return(standard_presence_table_spec_check)
}

doe_standard_presence_table_ref_spec = get_standard_presence_table_spec_check(doe_pk_tab, standards_data_table)
tut_standard_presence_table_ref_spec = get_standard_presence_table_spec_check(tut_pk_tab, standards_data_table_21)



###---Algorithmic detection with alternative peak finding---###
get_standard_presence_table_multiple = function(peak_tab, data_table) {
  #Create an empty nrow(peak_tab) x 12 numeric data frame to store algorithmic detection data 
  standard_presence_table_multiple = data.frame(matrix(FALSE, length(peak_tab$tab$V1), 12))
  rownames(standard_presence_table_multiple) = rownames(peak_tab$tab)
  colnames(standard_presence_table_multiple) = data_table[1, ]
  
  #Create and store the max spectra table for the peak tab
  peaks_max_spectra = get_max_absorbance_spectra_table(peak_tab)
  
  #Retrieve and store the RTs of the all peaks in the peak table as a vector
  peaks_rts = peak_tab$pk_meta[3,]
  
  #Iterate over each standard in data_table
  for(i in 1:length(standard_presence_table_multiple)){
    #Retrieve and store the RT and SD of the current standard from the data table as a numeric
    standard_rt = as.numeric(data_table[2,i])
    standard_sd = as.numeric(data_table[3,i])
    standard_max_spec = as.numeric(data_table[4,i])
    
    #Get the list of peaks within the current standard SD
    detected_rts = get_rt_in_sd_list(peaks_rts, standard_rt, standard_sd)
    
    if(length(detected_rts) > 0) { 
      
      if(length(detected_rts) == 1) {
        
        detected_rt_index = which(peaks_rts == detected_rts[1])
        
      } else {
        #Calculate the absolute difference between the standard RT and all detected RTs to see variation from standard RT
        detected_rts_diff = abs(detected_rts - standard_rt)
        
        #Find the index of the detected RT with the smallest variation from the standard RT
        closest_detected_rt_index = which.min(detected_rts_diff)
        
        detected_rt_index = which(peaks_rts == detected_rts[closest_detected_rt_index])
      }
      
      detected_rt_max_spec = as.numeric(peaks_max_spectra[detected_rt_index, 2])
      
      for(j in 1:length(peak_tab$tab$V1)) {
        if(peak_tab$tab[j, detected_rt_index] > 0){
          #if(abs(standard_max_spec - detected_rt_max_spec) < 20) {
            standard_presence_table_multiple[j, i] = TRUE
          #}
        }
      }
    } 
  }
  return(standard_presence_table_multiple)
}

doe_standard_presence_table_multiple = get_standard_presence_table_multiple(doe_pk_tab, standards_data_table)






#2/27/25(CE) - I messed around with different vectorized access methods but ultimately
# ended up refactoring with for loops. Most vectorized operations had unintended 
# side effects (mostly unwanted console output) that made print debugging impossible.
# Also found out that since R 3.4, for loops are basically just as efficient.
# ***PICK-UP POINT***: about to add reference spectra conditional after RT conditional
algorithmic_detection_peak_list = function(peak_list, stds_data_table) {
  # Create an empty nrow(peak_tab) x 12 numeric data frame to store algorithmic detection data
  alg_det_pk_list_df = data.frame(matrix(FALSE, length(peak_list), 12))
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
    #debug print to make sure we are accessing the correct data
    #print(paste0(std_name, " ", std_peak_rt, " ", std_peak_sd, " ", std_peak_spectra))
    
    #for current external standard, iterate over all subset runs
    for(sample_index in seq_along(peak_list)) {
      #Use name and positional indexing on sample_index from walk function to 
      #access peak data of individual runs 
      sample_peak_data = peak_list[[sample_index]]$`280`
      #confirm we are isolating the correct data with a print
      #print(sample_peak_data)
      
      sample_name = sample_peak_data$sample[1]
      
      #Pull current run RTs as a vector 
      sample_peak_rts = sample_peak_data$rt
      #confirm we are getting a vector with a print
      #print(sample_peak_rts)
      
      #Use multiple RT search to get vector of RTs from sample that fall within
      #sd range, then use logical indexing to select RT that varies least from
      #the RT of the current external standard
      sample_det_rt_list = get_rt_in_sd_list(sample_peak_rts, std_peak_rt, std_peak_sd)
      sample_det_rt = sample_det_rt_list[which.min(abs(sample_det_rt_list - std_peak_rt))]
      #confirm we are getting at least some matches with a print
      #print(sample_det_rt)
      
      #print(str(alg_det_pk_list_df))
      if(length(sample_det_rt) == 1) {
        #!!!NOTE!!!: since we are looping over the peak list rather than the peak
        # tab, we do not need to confirm the intensity of the peak to be 
        # above baseline to say we have detected the current standard. In the 
        # peak table, peaks are generalized and represented on a by run basis in 
        # $tab, with a 0 at the intersection of the current run and current peak
        # representing a lack of the peak within that run. In this way, using the
        # peak list may actually be more effecient, as conditionals can be taxing
        alg_det_pk_list_df[which(rownames(alg_det_pk_list_df) == sample_name), 
                           which(colnames(alg_det_pk_list_df) == std_name)] = TRUE
      }
    }
  }
  return(alg_det_pk_list_df)
}

doe_algorithmic_detection_peak_list = algorithmic_detection_peak_list(doe_pre_corrected_75_pks_egh, standards_data_table)
fre_algorithmic_detection_peak_list = algorithmic_detection_peak_list(fre_pre_corrected_75_pks_egh, standards_data_table)

#tut_algorithmic_detection_peak_list = algorithmic_detection_peak_list(tut_pre_corrected_75_pks_egh, standards_data_table_21)
pow_algorithmic_detection_peak_list = algorithmic_detection_peak_list(pow_pre_corrected_75_pks_egh, standards_data_table_21)



###---Creating a Presence Rate Table---###
# get_presence_rate_table(std_pres_tab_list, subset_meta_list)
#
# This method iterates through a list of standard presence tables for different subsets, and a corresponding list of 
# subset meta data to calculate and store presence rates for each standard across subsets.
#
# @param std_pres_tab_list - A list of standard_presence_table objects for all subsets
# @param subset_meta_list - A corresponding list of all meta data for all subsets
#
# @return - A data frame, structured to accumulate data iteratively, with the following columns:
#   - Site: Site name of subset.
#   - Year: Year of subset.
#   - Species: Ash species of subset.
#   - Standard: Name of standard.
#   - False: Number of treatments in the subset where the standard was absent.
#   - True: Number of treatments in the subset where the standard was present.
#   - Algorithmic_Detection_Rate: Rate of algorithmic detection for each standard (True/(True+False))
#
# @details
#   - Iterates through each subset in std_pres_tab_list and retrieves the corresponding metadata from subset_meta_list.
#   - Computes the number of true and false entries for each standard using the table function
#   - Constructs a new row for the presence rate table for each standard, and appends it to the data frame
#
get_presence_rate_table = function(std_pres_tab_list, subset_meta_list){
  #Initializing presence_rate_table as a data frame with empty, named columns
  presence_rate_table = data.frame(Site = character(0), Year = character(0), Species = character(0), 
                                   Standard = character(0), False = numeric(0), True = numeric(0), 
                                   Algorithmic_Detection_Rate = numeric(0), stringsAsFactors = FALSE) 
  #Iterate through each subset in std_pres_tab_list
  for(i in 1:length(std_pres_tab_list)){ 
    #Retrieve and store subset metadata
    subset_sitename = subset_meta_list[[i]][1,8]
    subset_year = subset_meta_list[[i]][1,7]
    subset_species = subset_meta_list[[i]][1,3]
    
    #Iterate through each standard in the current subset
    for(j in 1:length(std_pres_tab_list[[i]][1,])){
      #Calculate and store the number of TRUE and FALSE entries there are for the current standard in the current susbet
      standard_table = table(std_pres_tab_list[[i]][,j])
      trues = numeric()
      falses = numeric()
      
      #Check if standard_table has less than 2 entries (to handle edge cases)
      if(length(standard_table) != 2) { 
        #Check if we only have FALSE entries
        if(names(standard_table)[1] == "FALSE") {
          #Set the variable for number of TRUE entries to 0
          trues = 0 
          #Set the variable for the number of FALSE entries to the value in standard_table
          falses = standard_table[1]
        } else { #If we only have TRUE entries
          #Set the variable for number of FALSE entries to 0
          falses = 0
          #Set the variable for the number of TRUE entries to the value in standard_table
          trues = standard_table[1]
        }
      } else { #If standard_table has 2 entries
        #Set the variables to store the number of TRUE and FALSE entries to the corresponding values from standard_table
        falses = standard_table[1]
        trues = standard_table[2]
      }
      #Calculate the algorithmic detection rate for the current standard
      algorithmic_rate = trues / (falses + trues)
      #Create the new row for the presence rate table
      new_row = data.frame(Site = subset_sitename, Year = subset_year, Species = subset_species, 
                           Standard = colnames(std_pres_tab_list[[i]])[j], False = falses, True = trues, 
                           Algorithmic_Detection_Rate = algorithmic_rate, stringsAsFactors = FALSE)
      #Reset row names to avoid issues with rbind
      rownames(new_row) = NULL
      #Append new row to the data frame using rbind
      presence_rate_table = rbind(presence_rate_table, new_row)
    }
  }
  #Return the data frame to be stored in a variable
  return(presence_rate_table)
}

standard_presence_rate_table = get_presence_rate_table(standard_presence_table_list, subset_metadata_list)
#print(standard_presence_rate_table)





std_detection_rt_table = function(peak_tab, data_table) {
  std_det_rt_tab = data.frame(Sample_ID = character(0), Standard = character(0), Detected_RT = numeric(0), Intensity = numeric(0))
  
  #Retrieve and store the RTs of the all peaks in the peak table as a vector
  peaks_rts = peak_tab$pk_meta[3,]
  
  for(i in 1:length(data_table)) {
  
    #Retrieve and store the RT and SD of the search standard from the data table as a numeric
    standard_rt = as.numeric(data_table[2, i])
    standard_sd = as.numeric(data_table[3, i])
    
    #Get the list of peaks within the search standard SD
    detected_rts = get_rt_in_sd_list(peaks_rts, standard_rt, standard_sd)
    
    if(length(detected_rts) > 0) {
      if(length(detected_rts) == 1) {
        
        detected_rt_index = which(peaks_rts == detected_rts[1])
        
      } else {
        #Calculate the absolute difference between the standard RT and all detected RTs to see variation from standard RT
        detected_rts_diff = abs(detected_rts - standard_rt)
        
        #Find the index of the detected RT with the smallest variation from the standard RT
        closest_detected_rt_index = which.min(detected_rts_diff)
        
        detected_rt_index = which(peaks_rts == detected_rts[closest_detected_rt_index])
      }
      for(j in 1:length(peak_tab$tab$V1)) {
        if(peak_tab$tab[j, detected_rt_index] > 0) {
          #Updated this row to actually use sample ID instead of hplc name
          new_row = data.frame(Sample_ID = peak_tab$sample_meta$sample_name[j], Standard = data_table[1, i], Detected_RT = unname(peaks_rts[detected_rt_index]), Intensity = unname(peak_tab$tab[j, detected_rt_index]))
          std_det_rt_tab = rbind(std_det_rt_tab, new_row)
        }
      }
    }
  }
  row.names(std_det_rt_tab) = NULL
  return(std_det_rt_tab)
}

doe_std_detection_rt = std_detection_rt_table(doe_pk_tab, standards_data_table)
#print(doe_std_detection_rt)
#print(std_detection_rt_table(fre_pk_tab, standards_data_table))
#print(std_detection_rt_table(eas_pk_tab, standards_data_table))
#print(std_detection_rt_table(jen_pk_tab, standards_data_table))




###---Quantifying Peak AUC Contributions---###
# quantify_peaks(pk_list)
#
# This method iterates through a list of peak data objects, extracting and quantifying 
# peak AUC (Area Under the Curve) values for each run in the list, while normalizing them 
# as percentages of the total AUC for each run.
#
# @param pk_list - A list of peak data objects, each containing AUC values for individual peaks 
#                  in different sample runs.
#
# @return - A data frame structured to accumulate peak quantification data iteratively, 
# with the following columns:
#   - Run: Name of the sample run.
#   - PK1, PK2, ..., PKn: Normalized AUC percentages for each peak
#     are marked as NA.
#
# @details
#   - Iterates through each run in `pk_list` and extracts the AUC values for all peaks.
#   - Calculates the total AUC for peaks and normalizes AUC by dividing by the total.
#   - Marks AUC values ≤ 0.02 as NA to signify non-relevant or absent peaks.
#   - Appends the normalized AUC values for each peak to an accumulating data frame.
#
quantify_peaks = function(pk_tab) {
  #store peak tab `tab` data as a data frame, *removing all peaks with less than 3 occurrances
  tab_data = bind_rows(pk_tab$tab) 
  #%>% select(where(~ sum(.x > 0, na.rm = TRUE) >= 3))
  
  #call and store peak list from peak tab attributes
  pk_list = get(pk_tab$args$peak_list)
  
  #instantiate the return data frame
  quantify_peaks_df = data.frame()
  
  #looping over all runs in a peak list
  for(run_index in seq_along(pk_list)) {
    #extract and store run peak AUC data vector
    run_auc_data = pk_list[[run_index]]$`280`$`area`
    
    #sum AUC for all peaks (no exclusion of AUC < 0.02)
    run_total_auc = sum(run_auc_data)  # Changed to sum all AUC values
    
    #store run name for later use as row name
    run_name = names(pk_list)[run_index]

    #create a named list to fill with peak quantification data 
    new_row_list = list()
    
    #looping over all peaks in `tab` data
    for(pk_index in seq_along(tab_data)) {
      #make a character variable to give each peak a unique name
      pk_rt = as.character(pk_tab$pk_meta[3, pk_index])
      
      #store nested list of current peak data from peak list
      pk_data = pk_list[[run_index]]$`280`
      
      # Quantify all peaks (no exclusion of AUC < 0.02)
      if(run_total_auc > 0 && tab_data[run_index, pk_index] > 0) {  # Avoid division by or of zero
        #get the index for the area value for the current peak
        match_idx = which(pk_data$`height` == tab_data[run_index, pk_index])
        if(length(match_idx) > 0) {
          pk_height = pk_data$`area`[match_idx[1]]
          new_row_list[[pk_rt]] = round(((pk_height / run_total_auc) * 100), 6)
        } else {
          new_row_list[[pk_rt]] = 0
        }
      } else {
        new_row_list[[pk_rt]] = 0  # If total AUC is zero, mark as 0
      }
    }
    
    #convert list to data frame and assign run name as row name
    new_row_df = as.data.frame(new_row_list, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(new_row_df) = run_name
    
    #bind new row to greater data frame
    quantify_peaks_df = rbind(quantify_peaks_df, new_row_df)
  }
  
  #attach the argument as an attribute to the output
  attr(quantify_peaks_df, "args") = list(pk_tab = substitute(pk_tab))
  
  return(quantify_peaks_df)
}

#store call on doe subset
doe_quant_tab = quantify_peaks(doe_pk_tab)
#print(doe_quant_tab)
#check that all rows sum to 100
#print(rowSums(doe_quant_tab))

eas_quant_tab = quantify_peaks(eas_pk_tab)
#print(eas_quant_tab)

fre_quant_tab = quantify_peaks(fre_pk_tab)
#print(fre_quant_tab)

jen_quant_tab = quantify_peaks(jen_pk_tab)
#print(jen_quant_tab)



tut_quant_tab = quantify_peaks(tut_pk_tab)
#print(tut_quant_tab)

pow_quant_tab = quantify_peaks(pow_pk_tab)
#print(pow_quant_tab)

lee_quant_tab = quantify_peaks(lee_pk_tab)
#print(lee_quant_tab)

lit_quant_tab = quantify_peaks(lit_pk_tab)
#print(lit_quant_tab)




combined_peak_meta_table = function(pk_tab) {
  tab_df = bind_rows(pk_tab$tab)
  meta_df = bind_rows(pk_tab$sample_meta)
  #print(tab_df)
  #print(meta_df)
  combined_df = bind_cols(meta_df, tab_df)
  return(combined_df)
}

doe_combined_df = combined_peak_meta_table(doe_pk_tab)
#print(doe_combined_df)
write.csv(doe_combined_df, "doe_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

fre_combined_df = combined_peak_meta_table(fre_pk_tab)
#print(fre_combined_df)
write.csv(fre_combined_df, "fre_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

eas_combined_df = combined_peak_meta_table(eas_pk_tab)
#print(eas_combined_df)
write.csv(eas_combined_df, "eas_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

jen_combined_df = combined_peak_meta_table(jen_pk_tab)
#print(jen_combined_df)
write.csv(jen_combined_df, "jen_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")


tut_combined_df = combined_peak_meta_table(tut_pk_tab)
#(tut_combined_df)
write.csv(tut_combined_df, "tut_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

pow_combined_df = combined_peak_meta_table(pow_pk_tab)
#print(pow_combined_df)
write.csv(pow_combined_df, "pow_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

lee_combined_df = combined_peak_meta_table(lee_pk_tab)
#print(lee_combined_df)
write.csv(lee_combined_df, "lee_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

lit_combined_df = combined_peak_meta_table(lit_pk_tab)
#print(lit_combined_df)
write.csv(lit_combined_df, "lit_peak_and_meta.csv", row.names = FALSE, fileEncoding = "UTF-8")

