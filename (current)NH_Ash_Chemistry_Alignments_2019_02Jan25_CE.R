
###---2019---###

#2019 meta subset from above
hplc_names_subset_2019_df = hplc_names_subset_df[-c(209:626),]
hplc_names_subset_2019_df

#Now, to remove the standards from the list I have to set them as NULL as below
hplc2019_pre_append = hplc_list 
length(hplc2019_pre_append) #626
hplc2019_pre_append[c(209:626)] = NULL #list of elements to be removed from our list

length(hplc2019_pre_append) #208
length(hplc_names_subset_2019_df[,1]) #208

#need to make some manual adjustments to the dataframe otherwise metadata is difficult to pull below
hplc_names_subset_2019_df[13,3] = "Deer W Post 4-2"
hplc_names_subset_2019_df[25,3] = "Staff Pre 3-10"
hplc_names_subset_2019_df[40,3] = "S Pre 1-11"
hplc_names_subset_2019_df[208,3] = "Tut Pre 1-12F"

#adjust dataframe to account for new row reference numbers
hplc_names_subset_2019_df$row_number = seq(1:length(hplc_names_subset_2019_df[,1])) 

#use the stri_replace_asll_regex function to replace multiple strings
hplc_names_subset_2019_df$hplc_new_sample_names_subset = 
  stri_replace_all_regex(hplc_names_subset_2019_df$hplc_sample_names_subset,
                         pattern=c('Tut', 'Turtle',
                                   'Deer W', 'Dee W', 'D W', 'Deep W',
                                   'Deer G',
                                   'Staff', 'S1-11', 'S Pre', 'SPre',
                                   'Pre','pre','Post', 'post',
                                   ' '),
                         replacement=c('TS', 'TS',
                                       'DW', 'DW', 'DW', 'DW',
                                       'DG',
                                       'SF', 'SF 1-11', 'SF Pre','SF Pre',
                                       'Pr','Pr','Po', 'Po',
                                       '_'),
                         vectorize=FALSE)

#fix remaining strings
hplc_names_subset_2019_df$hplc_new_sample_names_subset = 
  stri_replace_all_regex(hplc_names_subset_2019_df$hplc_new_sample_names_subset,
                         pattern=c('TSF'
                         ),
                         replacement=c('TS'
                         ),
                         vectorize=FALSE)

#     4c) update dataframe with species and treatment information

hplc_names_subset_2019_updated_df = hplc_names_subset_2019_df #new variable
length(hplc_names_subset_2019_updated_df$hplc_new_sample_names_subset) #208

#use grepl to find instances of sites and make changes accordingly # | is "or"
hplc_names_subset_2019_updated_df$species <- ifelse(grepl("DG|TS",hplc_names_subset_2019_updated_df$hplc_new_sample_names_subset),
                                                    'Green','White') #basically if it's not a green site, it's a white site

#do same thing for pre- and post-samples
hplc_names_subset_2019_updated_df$sample_type <- ifelse(grepl("Pr",hplc_names_subset_2019_updated_df$hplc_new_sample_names_subset),
                                                        'Pre','Post') #if it's not a pre-sample, it's a post-sample


#Manual update of all the hplc_sample_names_subset that are still missing information
hplc_names_subset_2019_final_df = hplc_names_subset_2019_updated_df

#add 20_ in front of the samples that don't have it
# Subset the relevant column of the dataframe
subset_col19 <- hplc_names_subset_2019_final_df[, 5]

# Use the paste() function to add "19_" in front of each value
new_col19 <- paste("19_", subset_col19, sep = "")

# Assign the new column to the original dataframe
hplc_names_subset_2019_final_df[, 5] <- new_col19

#extract tree number by taking _,-, and numeric characters from our sample names subset
hplc_names_subset_2019_final_df$tree_id <- gsub("[^_\\x2D0-9\\-]", "", #x2D is the hex code for "-"; which doesn't work here for some reason, so added \\- at end of gsub  
                                                substr(hplc_names_subset_2019_final_df$hplc_new_sample_names_subset,
                                                       10, 14))

# create a new column "dbh_cat" based on "tree_id"
#set tree size category according to the last two digits of the string we extracted and cleaned:
#NOTE: 9May23: NAs are sample ids with missing data to be addressed at a later date
#tree number 1-3 = 3-6 cm DBH
#            4-6 = 6.01-9 cm DBH
#            7-9- = 9.01-12 cm DBH
#            10-12 = 12.01-15 cm DBH
hplc_names_subset_2019_final_df$dbh_cat <- ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2019_final_df$tree_id), 2, 3) %in% c("1", "2", "3"), "3-6",
                                                  ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2019_final_df$tree_id), 2, 3) %in% c("4", "5", "6"), "6.01-9",
                                                         ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2019_final_df$tree_id), 2, 3) %in% c("7", "8", "9"), "9.01-12",
                                                                ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2019_final_df$tree_id), 2, 3) %in% c("10", "11", "12"), "12.01-15", NA))))

#assign year
hplc_names_subset_2019_final_df$sample_year <- ifelse(grepl("19", hplc_names_subset_2019_final_df$hplc_new_sample_names_subset), "2019", #there are no 2019 samples in this subset
                                                      ifelse(grepl("20", hplc_names_subset_2019_final_df$hplc_new_sample_names_subset), "2019", "2019")) #last value was NA, but because these are all 2020 samples, I wanted to make sure everything was assigned

#assign site
hplc_names_subset_2019_final_df$site <- ifelse(grepl("DW", hplc_names_subset_2019_final_df$hplc_new_sample_names_subset), "Deerfield White",
                                               ifelse(grepl("DG", hplc_names_subset_2019_final_df$hplc_new_sample_names_subset), "Deerfield Green",
                                                      ifelse(grepl("SF", hplc_names_subset_2019_final_df$hplc_new_sample_names_subset), "Strafford Forest",
                                                             ifelse(grepl("TS", hplc_names_subset_2019_final_df$hplc_new_sample_names_subset), "Tuttle Swamp", 
                                                                    NA))))
#final dataframe
hplc_2019_df = hplc_names_subset_2019_final_df[,c(1,3,5:11)]
length(hplc_2019_df$hplc_new_sample_names_subset) #208

#---alignment---#

#-subset by site-#

#create subsets
deer_g = subset(hplc_2019_df,site =="Deerfield Green")
deer_w = subset(hplc_2019_df,site =="Deerfield White")
straff =subset(hplc_2019_df,site =="Strafford Forest")
tut19 = subset(hplc_2019_df,site =="Tuttle Swamp")

#hplc data
rt_2019_deer_g = hplc2019_pre_append
rt_2019_deer_w = hplc2019_pre_append
rt_2019_straff = hplc2019_pre_append
rt_2019_tut19 = hplc2019_pre_append

#get row ids for green and white subsets / create variable
deer_g$row_number = as.numeric(rownames(deer_g)) #length(deer_g$row_number) = 51
deer_w$row_number = as.numeric(rownames(deer_w)) #length(deer_w$row_number) = 52
straff$row_number = as.numeric(rownames(straff)) #length(straff$row_number) = 54
tut19$row_number = as.numeric(rownames(tut19)) #length(tut19$row_number) = 51

#remove non-deer_g rows
rt_2019_deer_g[c(deer_w$row_number,straff$row_number,tut19$row_number)] = NULL #list of elements to be removed from our list
length(rt_2019_deer_g) #length is now 51

#check names to confirm we only have deer_g samples
rt_2019_deer_g_sample_names = sapply(rt_2019_deer_g, function(x) attr(x,"sample_name")) #confirmed

#remove non-deer_w rows
rt_2019_deer_w[c(deer_g$row_number,straff$row_number,tut19$row_number)] = NULL #list of elements to be removed from our list
length(rt_2019_deer_w) #length is now 52

#check names to confirm we only have deer_w samples
rt_2019_deer_w_sample_names = sapply(rt_2019_deer_w, function(x) attr(x,"sample_name")) #confirmed

#remove non-straff rows
rt_2019_straff[c(deer_w$row_number,deer_g$row_number,tut19$row_number)] = NULL #list of elements to be removed from our list
length(rt_2019_straff) #length is now 54

#check names to confirm we only have straff samples
rt_2019_straff_sample_names = sapply(rt_2019_straff, function(x) attr(x,"sample_name")) #confirmed

#remove non-tut19 rows
rt_2019_tut19[c(deer_w$row_number,straff$row_number,deer_g$row_number)] = NULL #list of elements to be removed from our list
length(rt_2019_tut19) #length is now 51

#check names to confirm we only have tut19 samples
rt_2019_tut19_sample_names = sapply(rt_2019_tut19, function(x) attr(x,"sample_name")) #confirmed


###---Align peaks---###

#10May23: I cut the old commentary here about the green ash chromatograms from 2020.
#that code/commentary can be found in previous versions of this file

#this suggests a simple way to drop NULLs from our list of lists: https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
#rt_2020_white_only_pre2 = compact(rt_2020_white_only_pre2)

#lets align chromatograms by site
#first, lets preprocess after subsetting, which is unlike what we did before (we previously aligned all by number in the list)?


#---pre-process data---#

#NOTE: There are duplicates in some of these subsets. Make sure that those that are removed have been retained or returned to 
#the correct position within their dataframe

#2020 - Doe Farm - Green Ash

#before we align, there are some duplicate and extra samples that were somehow retained. Lets remove them.
#There were some remove duplicates from processed data
rt_2020_doe_nodup <- rt_2020_doe[!duplicated(rt_2020_doe)]
length(rt_2020_doe_nodup) #106

#I could not figure out how to pull these out programmatically so I did it manually. The three extra observations
#in doe_pre_nodup are 18 (duplicate),24 (Fc), and 55 (Jf).  
rt_2020_doe_nodup[c(18,24,55)] = NULL 
rt_2020_doe_nodup = compact(rt_2020_doe_nodup)
length(rt_2020_doe_nodup) #103

#remove duplicates from metadata
hplc_2020_df_doe_nodup = doe[!duplicated(doe$hplc_sample_names_subset),]

#change column name to be consistent in list and dataframe
colnames(hplc_2020_df_doe_nodup) <- c("hplc_name","sample_name","tree_species","sample_type","treatment","tree_size","year","site")
length(hplc_2020_df_doe_nodup$sample_name) #same length as doe_pre_nodup

#according to the ?find_peaks help
#It is recommended to do pre-processing using the preprocess function before peak detection. Overly high chromatographic resolution 
#can sometimes lead to peaks being split into multiple segments. In this case, it is recommended to increase the smooth_window or 
#reduce the resolution along the time axis by adjusting the dim2 argument before proceeding with further analyses.

#I originally thought I wanted a very fine time resolution but perhaps I do not, I am going to try and pre-process at 0.05 and then 0.1 to see how this
#alters the number of peaks found. I want to see ~70-90 as the HPLC output originally suggested
doe_pre = preprocess(rt_2020_doe_nodup, 
                     dim1=seq(0,60,.05),
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE,
                     parallel=FALSE)#, #I don't know why I removed this previously
#mc.cores= 6,
#mc.cleanup=TRUE) #https://www.rdocumentation.org/packages/parallel/versions/3.4.0/topics/mclapply

#align data
doe_pre_05_corrected_75 = 
  correct_rt(
    doe_pre, 
    what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
    lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
    models = NULL,
    reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
    alg = c("vpdtw"),
    init.coef = c(0, 1, 0),
    n.traces = NULL,
    n.zeros = 0,
    scale = TRUE, #default is FALSE
    plot_it = FALSE,
    penalty = 1, #default 5
    maxshift = 75, #
    verbose = FALSE,
    progress_bar = FALSE)

#now I am wondering if I should be using ptw?

###---Peak finding---###
#adjusting the smooth window changes the number of peaks output/detected below
#0.75 = 55 peaks;0.05 = 60 peaks; 0.04 = 72 peaks (seems close to HPLC determination); 0.035 = 80 peaks (seems close to HPLC determination)  
#0.025 = 94 peaks; 0.02 = 103 peaks; 0.015 = 121 peaks

#use this subset for the report
doe_pre_corrected_75_pks_egh <- get_peaks(doe_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.035, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000) #works

#plot peak finding
plot(doe_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035) looks like garbage but is it actually?

###---Creating a peak table---###
doe_pre_corrected_75_egh_peak <- get_peaktable(doe_pre_corrected_75_pks_egh, response = 'height') #

#subset the metadata dataframe for doe
doe_meta = subset(hplc_2020_df_doe_nodup, site=="Doe Farm")
length(doe_meta$site) #103

#attach metadata
doe_pk_tab <- attach_metadata(peak_table = doe_pre_corrected_75_egh_peak, 
                              metadata = doe_meta, 
                              column = "hplc_name") #this works! 

#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
doe_pk_tab_df = as.data.frame(do.call(cbind, doe_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
doe_pk_tab_rt_df = as.data.frame(do.call(cbind, doe_pk_tab[2]))





