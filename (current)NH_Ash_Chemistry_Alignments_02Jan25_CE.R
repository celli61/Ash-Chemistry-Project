##################################################################################################################
### Project: Testing the ontogenetic resistance hypothesis in green and white ash across four size             ###
###          classes in New Hampshire                                                                          ###
###                                                                                                            ###
### Dates of sampling: Year 1 Presample: XXJun19   Postsample: XXJul19 CHECK THIS IN NOTEBOOK IN DESK IN OFFICE###
###                    Year 2 Presample: XXJun20   Postsample: XXJul20                                         ###
###                    Year 3 Presample: XXJun21   Postsample: XXJul21                                         ###
###                                                                                                            ###
### Goal of project: To determine how composition of phloem defensive metabolites changes after emerald        ###
###                  ash borer attack, or induction with the plant hormone methyl jasmonate (Meja)             ###
###                                                                                                            ###
### Predictions: 1) EAB attack will influence the identity and abundance of defensive metabolites              ###
###              2) Induction with Meja will influence the  identity and abundance of defensive metabolites    ###
###              3) Tree size class will influence pre- and post-induced identity and abundance of defensive   ###
###                 metabolites                                                                                ###
###              4) Tree species will influence pre- and post-induced identity and abundance of defensive      ###
###                                                                                                            ###
##################################################################################################################

#---NOTES---#

#26May23: I copied all the code from NH_Ash_Chemistry26May23 and put it here, except I cut the NMDS at the end
#         The goal is to put all the processing and alignment code in this separate R file, to make it easier
#         to find and work with the relevant code for each step of this analysis, which right now will be
#         1) preprocessing/alignment, 2) NMDS, 3) Random Forest Analysis

#load libraries
pacman::p_load(chromConverter, usethis, chromatographR, stringi, tidyverse, pbapply, patchwork, gtools, ggrepel, data.table, parallel)
#chromConverter - read data from propriatary shimadzu format
#chromatographR - post-process and analyze spectra converted from chromConverter
#stringi - to replace multiple patterns (better than multiple gsub())
#pbapply - to allow parallel processing on windows
#pCompare - to allow floating point comparison with tolerance

#initialize a cluster of size (#machine cpu cores) - 1 for parallel processing
cl = makeCluster(detectCores() - 1)
###---read data---###

#Separate comma format (note that the file folder has an typo, fix at some point)

#setwd("/Users/toddjohnson/Sync/Box\ Sync/OneDrive\ -\ LSU\ AgCenter/Box\ Sync\ Files/Data/UNH\ HPLC\ data/Ash\ Phenolics/2020/Seperate_Comma/")

#for running on my laptop
#run_2019_2020_separate_comma <- read_chroms("/Users/colto/Documents/AshPhenolics/2020/Seperate_Comma/",
#                                            read_metadata = TRUE,
#                                            format_in = "shimadzu_dad", cl = cl)

run_2019_2020_separate_comma <- read_chroms("/Users/colto/Documents/AshPhenolics/2020/Seperate_Comma/",
                                            read_metadata = TRUE,
                                            format_in = "shimadzu_dad", cl = cl)


###---Pre-process data---###
#Pre-preprocessing of data is done to remove baseline noise (i.e., peaks that are not related to our sample) in the chromatogram,
#as well as to look at relevant retention times

#exclude standards from peak alignment (changed to now removing before baseline correction 15Feb23) as it may cause issues and reduce the quality of the alignment
#list of standards to be removed

#Run 1 - 001 Ligustrocide
#Run 1 - 002 Eleutheroside B
#Run 1 - 003 Oleuropein
#Run 1 - 004 Verbascoside
#Run 1 - 005 Pinoresinol
#Run 1 - 006 Luteolin
#Run 1 - 007 Chlorogenic Acid
#Run 1 - 008 Vanillic Acid
#Run 1 - 009 Ferulic Acid
#Run 1 - 010 Coumaric Acid
#Run 1 - 011 Gallic Acid
#Run 1 - 012 Standard Mix 1 mg/ml
#Run 1 - 013 Standard Mix 0.5 mg/ml
#Run 1 - 014 Standard Mix 0.25 mg/ml
#Run 1 - 015 Standard Mix 0.125 mg/ml
#Run 1 - 016 Standard Mix 2 1 mg/ml
#Run 1 - 017 Standard Mix 2 0.5 mg/ml
#Run 1 - 018 Standard Mix 2 0.25 mg/ml
#Run 1 - 019 Standard Mix 2 0.125 mg/ml
#Run 1 - 045 Check 1
#Run 1 - 046 Null 1
#Run 1 - 047 Ligustrocide
#Run 1 - 048 Standards 0.06 mg/ml
#Run 1 - 059 Ligustroside Standard Third
#Run 1 - 060 Ligustroside 0.25 mg/ml
#Run 1 - 061 Ligustroside 0.06 mg/ml
#Run 1 - 077 Check 2
#Run 1 - 103 Check 3
#Run 1 - 129 Check 4
#Run 1 - 130 Eleutheroside B 2 mg/ml
#Run 1 - 155 Check 5
#Run 1 - 181 Check 6
#Run 1 - 207 Check 7
#Run 1 - 233 Check 8
#Run 1 - 239 Standard Mix 3 0.5 mg/ml
#Run 1 - 240 Standard Mix 3 0.25 mg/ml
#Run 1 - 241 Standard Mix 3 0.125 mg/ml
#Run 1 - 242 Standard Mix 3 0.06 mg/ml

#Run 2 - 002 Check Standard 1 #249 corrected
#Run 2 - 027 Check Standard 2 #274
#Run 2 - 052 Check Standard 3 #299
#Run 2 - 013 Check Standard 4 #260
#Run 2 - 039 Check Standard 5 #286
#Run 2 - 065 Check Standard 6 #312
#Run 2 - 091 Check Standard 7 #338
#Run 2 - 117 Check Standard 8 #364
#Run 2 - 143 Check Standard 9 #390

#Run 3 - 002 Check Standard 1 #470
#Run 3 - 003 Apigenin 0.01 mg/ml
#Run 3 - 004 Apigenin 0.1 mg/ml
#Run 3 - 005 Apigenin 0.2 mg/ml #473
#Run 3 - 028 Check Standard 2 #496
#Run 3 - 054 Check Standard 3 #522
#Run 3 - 080 Check Standard 4 #548
#Run 3 - 106 Check Standard 5 #574
#Run 3 - 132 Check Standard 6 #600
#Run 3 - 158 Check Standard 7 #626
#Run 3 - 184 Check Standard 8 #652
#Run 3 - 210 Check Standard 9 #678
#Run 3 - 220 Check Standard 10 #686

#make new variable to prevent modification raw data
hplc_list = run_2019_2020_separate_comma# hplc_pre_append

#pull the names of each list item (these are the file names of the chromatograms)
hplc_names = names(hplc_list)

#get a vector with a number for each element in the vector (this is so I know which rows to call to remove standards prior to pre-processing and alignment)
hplc_length = seq(1:length(hplc_names))

#get sample names - this pulls the actual sample id as was input into the HPLC
hplc_sample_names = sapply(run_2019_2020_separate_comma, function(x) attr(x,"sample_name")
)

#For some reason I was having trouble binding this to the rest of the dataframe so i set it as a vector
hplc_sample_names = as.vector(hplc_sample_names)

hplc_length = as.vector(hplc_length) #need to change this

#bind together in a dataframe
hplc_names_df = cbind.data.frame(hplc_sample_names,hplc_length)

#this allows me to add hplc_names to the dataframe
hplc_names_df$hplc_names = hplc_names

#Copy standards into list to process separately
hplc_list_standards = hplc_list[c(1:19,45:48,59:61,77,103,129,130,155,181,207,233,239:242,
                                  247,272,297,323,349,375,401,427,453,
                                  468:471,494,520,546,572,598,624,650,676,686)]

#Remove standards from the list by setting them as NULL
hplc_list[c(1:19,45:48,59:61,77,103,129,130,155,181,207,233,239:242,
            247,272,297,323,349,375,401,427,453,
            468:471,494,520,546,572,598,624,650,676,686)] = NULL #list of elements to be removed from our list

#pre-processing our lists without standards

#dim1 changes the x-axis and essentially slices out the portions you don't want to keep.
#hopefully this makes alignment easier; dim1=seq(start minute,end minute,interval in between minutes)
#dim2 selects the lambdas or wavelengths to include; dim2=seq(190,400,2) 190 nm to 400 nm by 2 nm

#this plot can be used to look at the peaks
#matplot(sapply(hplc_list_pre1_.5,function(x)x[,"280"]),type='l')
#legend("topright", legend="0.5", bty = "n")

#Note 1- as of 26May23 I removed a bunch of preprocessing trials were I adjusted the number of points in between minute
#        this is the third number in dim1. I found that 0.01 is the best slicing for these peaks
#
#.    2- I previously pre-processed everything according to location in the hplc_list. This has changed and I don't 
#        pre-process until later in the script. Commented below out for posterity.

#preprocess at 0.01. If I "zoom in" I see a lot of recovered smaller peaks compared to other slicings. I will now try
#to pre-process everything at 0.01; I also modified the number of cores to be used from 2 (default) to 6 to zoom zoom
#hplc_list_pre1 = preprocess(hplc_list[1:200], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)
#hplc_list_pre2 = preprocess(hplc_list[201:250], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)
#hplc_list_pre3 = preprocess(hplc_list[251:260], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)
#hplc_list_pre4 = preprocess(hplc_list[261:265], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)
#hplc_list_pre5 = preprocess(hplc_list[266:268], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)
#hplc_list_pre6 = preprocess(hplc_list[269], dim1=seq(10,60,.5),dim2=seq(190,400,2))#this one doesn't work, excluding for now
#hplc_list_pre6 = preprocess(hplc_list[270:300], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)
#hplc_list_pre7 = preprocess(hplc_list[301:626], dim1=seq(10,60,.01),dim2=seq(190,400,2),mc.cores=6)

#because I pre-processed chromatograms separately, I need to append them all back into a list
#because I can't figure out how to append them all at once, I have appended to a growing list
#hplc_pre_append = append(hplc_list_pre1,hplc_list_pre2)
#hplc_pre_append = append(hplc_pre_append,hplc_list_pre3)
#hplc_pre_append = append(hplc_pre_append, hplc_list_pre4)
#hplc_pre_append = append(hplc_pre_append,hplc_list_pre5)
#hplc_pre_append = append(hplc_pre_append,hplc_list_pre6)
#hplc_pre_append = append(hplc_pre_append,hplc_list_pre7)

#check length of list to confirm we appended everything correctly 626 - 1 chromatogram above = 625
#length(hplc_pre_append) #625

###---Align pre-processed standards---###

#06Feb24_CE: I basically copied exactly how the samples are subset, aligned, and pre-processed to produce a peak table of standards

#assign new variable new to prevent changes to old
hplc_pre_append_standards = hplc_list_standards #length 60

hplc_names_subset_standards = names(hplc_pre_append_standards) #get names of files

hplc_sample_names_subset_standards = sapply(hplc_pre_append_standards, function(x) attr(x,"sample_name")) #get sample_names attribute
hplc_sample_names_subset_standards = as.vector(hplc_sample_names_subset_standards) #set as a vector

hplc_length_subset_standards = seq(1:length(hplc_pre_append_standards)) #assign a number for each row
hplc_length_subset_standards = as.vector(hplc_length_subset_standards) #set as a vector

hplc_names_subset_standards_df = cbind.data.frame(hplc_names_subset_standards,hplc_length_subset_standards) #bind together in a dataframe
hplc_names_subset_standards_df$hplc_sample_names_subset_standards = hplc_sample_names_subset_standards #add names


#sub-setting hplc_list_standards by standards and standard mix to aid in alignment 
standards_pre_ligustrocide = hplc_list_standards[c(22,24:26)]
standards_pre_eleutheroside_b = hplc_list_standards[c(2,30)]
standards_pre_oleuropein = hplc_list_standards[c(3)]
standards_pre_verbascoside = hplc_list_standards[c(4)]
standards_pre_pinoresinol = hplc_list_standards[c(5)]
standards_pre_luteolin = hplc_list_standards[c(6)]
standards_pre_chlorogenic_acid = hplc_list_standards[c(7)]
standards_pre_vanillic_acid = hplc_list_standards[c(8)]
standards_pre_ferulic_acid = hplc_list_standards[c(9)]
standards_pre_coumaric_acid = hplc_list_standards[c(10)]
standards_pre_gallic_acid = hplc_list_standards[c(11)]
standards_pre_standard_mix = hplc_list_standards[c(12:20,27:29,31:48,52:60)]
standards_pre_apagenin = hplc_list_standards[c(49:51)]


#pre-processing subset standards
rt_standards_pre_ligustrocide = preprocess(standards_pre_ligustrocide, 
                                           dim1=seq(0,54,.01),
                                           dim2=seq(190,400,1),
                                           remove.time.baseline = TRUE,
                                           cl=cl, show_progress = FALSE)

rt_standards_pre_eleutheroside_b = preprocess(standards_pre_eleutheroside_b, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_oleuropein = preprocess(standards_pre_oleuropein, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_verbascoside = preprocess(standards_pre_verbascoside, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_pinoresinol = preprocess(standards_pre_pinoresinol, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_luteolin = preprocess(standards_pre_luteolin, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_chlorogenic_acid = preprocess(standards_pre_chlorogenic_acid, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_vanillic_acid = preprocess(standards_pre_vanillic_acid, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_ferulic_acid = preprocess(standards_pre_ferulic_acid, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_coumaric_acid = preprocess(standards_pre_coumaric_acid, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

rt_standards_pre_gallic_acid = preprocess(standards_pre_gallic_acid, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)

#rt_standards_pre_standard_mix = preprocess(standards_pre_standard_mix, 
#                                          dim1=seq(0,54,.01),
#                                          dim2=seq(190,400,1),
#                                          remove.time.baseline = TRUE,
#                                          cl=cl, show_progress = FALSE)

rt_standards_pre_apagenin = preprocess(standards_pre_apagenin, 
                                          dim1=seq(0,54,.01),
                                          dim2=seq(190,400,1),
                                          remove.time.baseline = TRUE,
                                          cl=cl, show_progress = FALSE)



#---alignment---#


#plotting ligustrocide 
plot_chroms(rt_standards_pre_ligustrocide, lambdas = c('280'), show_legend = TRUE)

#corrected_rt_standards_pre_ligustrocide = rt_standards_pre_ligustrocide
corrected_rt_standards_pre_ligustrocide = correct_rt(
  rt_standards_pre_ligustrocide, 
  what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
  lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
  models = NULL,
  reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
  alg = c("vpdtw"),
  init.coef = c(0, 1, 0),
  n.traces = NULL,
  n.zeros = 0,
  scale = FALSE, #default is FALSE
  plot_it = FALSE,
  penalty = 1, #default 5
  maxshift = 75, #
  verbose = FALSE,
  progress_bar = FALSE, 
  cl = cl)

#plotting ligustrocide after alignment 
plot_chroms(corrected_rt_standards_pre_ligustrocide, lambdas = c('280'), show_legend = FALSE) #looks perfect


#plotting eleutheroside b before alignment 
plot_chroms(rt_standards_pre_eleutheroside_b, lambdas = c('280'), show_legend = FALSE) #pre-processing looks good

corrected_rt_standards_pre_eleutheroside_b = correct_rt(
  rt_standards_pre_eleutheroside_b, 
  what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
  lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
  models = NULL,
  reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
  alg = c("vpdtw"),
  init.coef = c(0, 1, 0),
  n.traces = NULL,
  n.zeros = 0,
  scale = FALSE, #default is FALSE
  plot_it = FALSE,
  penalty = 1, #default 5
  maxshift = 75, #
  verbose = FALSE,
  progress_bar = FALSE,
  cl = cl)

#plotting eleutheroside b after alignment 
plot_chroms(corrected_rt_standards_pre_eleutheroside_b, lambdas = c('280'), show_legend = FALSE) #looks perfect


#plotting oleuropein
plot_chroms(rt_standards_pre_oleuropein, lambdas = c('280'), show_legend = FALSE) #pre-processing looks good

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_oleuropein = c(rt_standards_pre_oleuropein, copy(rt_standards_pre_oleuropein))


#plotting verbascoside 
plot_chroms(rt_standards_pre_verbascoside, lambdas = c('280'), show_legend = FALSE) #pre-processing looks good

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_verbascoside = c(rt_standards_pre_verbascoside, copy(rt_standards_pre_verbascoside))


#plotting pinoresinol 
plot_chroms(rt_standards_pre_pinoresinol, lambdas = c('280'), show_legend = FALSE) #pre-processing looks good

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_pinoresinol = c(rt_standards_pre_pinoresinol, copy(rt_standards_pre_pinoresinol))


#plotting luteolin 
plot_chroms(rt_standards_pre_luteolin, lambdas = c('280'), show_legend = FALSE) 

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_luteolin = c(rt_standards_pre_luteolin, copy(rt_standards_pre_luteolin))


#plotting chlorogenic acid 
plot_chroms(rt_standards_pre_chlorogenic_acid, lambdas = c('280'), show_legend = FALSE) 

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_chlorogenic_acid = c(rt_standards_pre_chlorogenic_acid, copy(rt_standards_pre_chlorogenic_acid))


#plotting vanillic acid
plot_chroms(rt_standards_pre_vanillic_acid, lambdas = c('280'), show_legend = FALSE) 

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_vanillic_acid = c(rt_standards_pre_vanillic_acid, copy(rt_standards_pre_vanillic_acid))


#plotting ferulic acid before alignment 
plot_chroms(rt_standards_pre_ferulic_acid, lambdas = c('280'), show_legend = FALSE) 

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_ferulic_acid = c(rt_standards_pre_ferulic_acid, copy(rt_standards_pre_ferulic_acid))


#plotting coumaric acid 
plot_chroms(rt_standards_pre_coumaric_acid, lambdas = c('280')) 

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_coumaric_acid = c(rt_standards_pre_coumaric_acid, copy(rt_standards_pre_coumaric_acid))


#plotting gallic acid 
plot_chroms(rt_standards_pre_gallic_acid, lambdas = c('280')) 

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_gallic_acid = c(rt_standards_pre_gallic_acid, copy(rt_standards_pre_gallic_acid))


#plotting standard mix before alignment 
#plot_chroms(rt_standards_pre_standard_mix, lambdas = c('280'), show_legend = FALSE)

#corrected_rt_standards_pre_standard_mix = rt_standards_pre_standard_mix
#corrected_rt_standards_pre_standard_mix = correct_rt(
#  rt_standards_pre_standard_mix, 
#  what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
#  lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
#  models = NULL,
#  reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
#  alg = c("vpdtw"),
#  init.coef = c(0, 1, 0),
#  n.traces = NULL,
#  n.zeros = 0,
#  scale = TRUE, #default is FALSE
#  plot_it = FALSE,
#  penalty = 1, #default 5
#  maxshift = 75, #
#  verbose = FALSE,
#  progress_bar = FALSE)

#plotting standard mix after alignment
#plot_chroms(corrected_rt_standards_pre_standard_mix, lambdas = c('280'), show_legend = FALSE) 


#plotting apagenin before alignment
plot_chroms(rt_standards_pre_apagenin, lambdas = c('280'), show_legend = FALSE)

#storing in a new variable and appending a deep copy of the single chromatogram using copy()
corrected_rt_standards_pre_apagenin = c(rt_standards_pre_apagenin, copy(rt_standards_pre_apagenin))



###---Peak finding---###

#Ligustrocide
standards_ligustrocide_peaks = get_peaks(corrected_rt_standards_pre_ligustrocide, lambdas = c('280'), fit = "egh", amp_thresh = 5000, show_progress = FALSE, cl = cl) 

standards_ligustrocide_peak_table = get_peaktable(standards_ligustrocide_peaks, response = "height")

standards_ligustrocide_peak_table_df = standards_ligustrocide_peak_table$tab

#the common peak in all 4 ligustrocide treatments is peak v4, with peak retention time of 48.4, and range of 48.12 to 48.90 and a standard deviation of 0.07
print(standards_ligustrocide_peak_table$pk_meta)
print(standards_ligustrocide_peak_table_df)


#Eletheroside B
standards_eleutheroside_b_peaks = get_peaks(corrected_rt_standards_pre_eleutheroside_b, lambdas = c('280'), fit = "egh", amp_thresh = 21000, show_progress = FALSE, cl = cl) 
standards_eleutheroside_b_peak_table = get_peaktable(standards_eleutheroside_b_peaks, response = 'height')
standards_eleutheroside_b_peak_table_df = standards_eleutheroside_b_peak_table$tab

print(standards_eleutheroside_b_peak_table$pk_meta)
print(standards_eleutheroside_b_peak_table_df)

#Oleuropein
standards_oleuropein_peaks = get_peaks(corrected_rt_standards_pre_oleuropein, lambdas = c('280'), fit="egh", amp_thresh = 50000, show_progress = FALSE, cl = cl)
standards_oleuropein_peak_table = get_peaktable(standards_oleuropein_peaks, response = 'height')

print(standards_oleuropein_peak_table$pk_meta)

#Verbascoside
standards_verbascoside_peaks = get_peaks(corrected_rt_standards_pre_verbascoside, lambdas = c('280'), fit="egh", amp_thresh = 200000, show_progress = FALSE, cl = cl)
standards_verbascoside_peak_table = get_peaktable(standards_verbascoside_peaks, response = 'height')

print(standards_verbascoside_peak_table$pk_meta)

#Pinoresinol
standards_pinoresinol_peaks = get_peaks(corrected_rt_standards_pre_pinoresinol, lambdas = c('280'), fit="egh", amp_thresh = 200000, show_progress = FALSE, cl = cl)
standards_pinoresinol_peak_table = get_peaktable(standards_pinoresinol_peaks, response = 'height')

print(standards_pinoresinol_peak_table$pk_meta)

#Luteolin
standards_luteolin_peaks = get_peaks(corrected_rt_standards_pre_luteolin, lambdas = c('280'), fit="egh", amp_thresh = 17000, show_progress = FALSE, cl = cl)
standards_luteolin_peak_table = get_peaktable(standards_luteolin_peaks, response = 'height')

print(standards_luteolin_peak_table$pk_meta)

#Chlorogenic acid
standards_chlorogenic_acid_peaks = get_peaks(corrected_rt_standards_pre_chlorogenic_acid, lambdas = c('280'), fit="egh", amp_thresh = 600000, show_progress = FALSE, cl = cl)
standards_chlorogenic_acid_peak_table = get_peaktable(standards_chlorogenic_acid_peaks, response = 'height')

print(standards_chlorogenic_acid_peak_table$pk_meta)

#Vanillic acid
standards_vanillic_acid_peaks = get_peaks(corrected_rt_standards_pre_vanillic_acid, lambdas = c('280'), fit="egh", amp_thresh = 800000, show_progress = FALSE, cl = cl)
standards_vanillic_acid_peak_table = get_peaktable(standards_vanillic_acid_peaks, response = 'height')

print(standards_vanillic_acid_peak_table$pk_meta)

#Ferulic acid
standards_ferulic_acid_peaks = get_peaks(corrected_rt_standards_pre_ferulic_acid, lambdas = c('280'), fit="egh", amp_thresh = 1600000, show_progress = FALSE, cl = cl)
standards_ferulic_acid_peak_table = get_peaktable(standards_ferulic_acid_peaks, response = 'height')

print(standards_ferulic_acid_peak_table$pk_meta)

#Coumaric acid
standards_coumaric_acid_peaks = get_peaks(corrected_rt_standards_pre_coumaric_acid, lambdas = c('280'), fit="egh", amp_thresh = 2000000, show_progress = FALSE, cl = cl)
standards_coumaric_acid_peak_table = get_peaktable(standards_coumaric_acid_peaks, response = 'height')

print(standards_coumaric_acid_peak_table$pk_meta)

#Gallic acid
standards_gallic_acid_peaks = get_peaks(corrected_rt_standards_pre_gallic_acid, lambdas = c('280'), fit="egh", amp_thresh = 470000, show_progress = FALSE, cl = cl)
standards_gallic_acid_peak_table = get_peaktable(standards_gallic_acid_peaks, response = 'height')

print(standards_gallic_acid_peak_table$pk_meta)

#Apagenin
standards_apagenin_peaks = get_peaks(corrected_rt_standards_pre_apagenin, lambdas = c('280'), fit="egh", amp_thresh = 11000, show_progress = FALSE, cl = cl)
standards_apagenin_peak_table = get_peaktable(standards_apagenin_peaks, response = 'height')
standards_apagenin_peak_table_df = standards_apagenin_peak_table$tab

print(standards_apagenin_peak_table$pk_meta)
print(standards_apagenin_peak_table_df)










#Combine all standards peak lists
standards_peaks = c(standards_gallic_acid_peaks, standards_chlorogenic_acid_peaks, standards_vanillic_acid_peaks, 
                    standards_eleutheroside_b_peaks, standards_coumaric_acid_peaks, standards_ferulic_acid_peaks,
                    standards_verbascoside_peaks, standards_oleuropein_peaks, standards_pinoresinol_peaks, 
                    standards_ligustrocide_peaks, standards_luteolin_peaks, standards_apagenin_peaks)
print(standards_peaks)
for (attr_name in names(attributes(standards_gallic_acid_peaks))) {
  if(attr_name != "names"){
    attr_value <- attributes(standards_gallic_acid_peaks)[[attr_name]]
    attr(standards_peaks, attr_name) <- attr_value
  }
}
attr(standards_peaks, "chrom_list") = "pre_standards"
print(standards_peaks)

#Make a peak table out of the combined standard peak list's
standards_peak_table = get_peaktable(standards_peaks, response = 'height') #not all too useful, attributes get lost so no 
                                                                         #further analysis


#combining external standard subsets again for generalized peak finding and analysis 
standards_pre = c(standards_pre_gallic_acid, standards_pre_chlorogenic_acid, 
                  standards_pre_vanillic_acid, standards_pre_eleutheroside_b, 
                  standards_pre_coumaric_acid, standards_pre_ferulic_acid, 
                  standards_pre_verbascoside, standards_pre_oleuropein, 
                  standards_pre_pinoresinol, standards_pre_ligustrocide, 
                  standards_pre_luteolin, standards_pre_apagenin)

rt_standards_pre = preprocess(standards_pre, 
           dim1=seq(0,54,.05),
           dim2=seq(274,284,2),
           remove.time.baseline = TRUE,
           cl=cl)

plot_chroms(rt_standards_pre, lambdas = c('280'), show_legend = FALSE)












###---Align pre-processed peaks---###

#26May23 UPDATE THIS AND THEN REMOVE THIS COMMENT
#     3) Create new columns with relevant sample information (e.g., site, species of ash) and attempt to subset my data further
#        it is possible that chemistry is highly variable between species or other experimental factors.
#        I will be moving forward with subsetting data from run 3 

#Because the only relevant sample information is contained in the "sample_name" attribute, I will:
#
#     4a) Create a new data frame as above with row number, sample name, file name


#assign new variable new to prevent changes to old
hplc_pre_append = hplc_list #length 626

hplc_names_subset = names(hplc_pre_append) #get names of files

hplc_sample_names_subset = sapply(hplc_pre_append, function(x) attr(x,"sample_name")) #get sample_names attribute
hplc_sample_names_subset = as.vector(hplc_sample_names_subset) #set as a vector

hplc_length_subset = seq(1:length(hplc_pre_append)) #assign a number for each row
hplc_length_subset = as.vector(hplc_length_subset) #set as a vector

hplc_names_subset_df = cbind.data.frame(hplc_names_subset,hplc_length_subset) #bind together in a dataframe; interestingly, on my Mac there is no error with as.data.frame(); when I run this on the PC I get an error, Stackoverflow suggested this cbind.data.frame() here: https://stackoverflow.com/questions/48845309/row-names-is-not-a-character-vector-of-length
hplc_names_subset_df$hplc_sample_names_subset = hplc_sample_names_subset #add names

#     4b) Standardize site and sample type abbreviations

#Four sites were used in 2020:
#Doe Farm - Green Ash
#French Conservation Area - Green Ash
#East Foss Farm - White Ash
#Jennings Forest - White Ash

#Upon looking at "hplc_sample_names_subset" I realized that the naming conventions shifted halfway through run 3.
#[1:271] is named in a somewhat haphazard fashion with sites named multiple things in different numbers of characters
#[272:625] sample names are standardized so that they are easier to refer to and understand.

#Thus to process these samples into something usable, I am going to process these samples semi-separately.

#Another observation here (and may be influential to how well these samples have been able to be corrected) is that a number of samples from
#2019 were apparently re-ran (I forgot about this) since Chad improved the runs. I need to first remove these similar to how I removed
#the standards from the samples. Remove everything except for the four sites listed above


###---2019---###

#non-2020 sites: 1:208

#create a 2019 subset
hplc_names_subset_2019_df = hplc_names_subset_df[-c(209:626),]

#Need to do

###---2020---###

#remove 2019 sites to create a dataframe that is just 2020
hplc_names_subset_2020_df = hplc_names_subset_df[-c(1:208),]

#Now, to remove the standards from the list I have to set them as NULL as below
hplc_pre_append[c(1:208)] = NULL #list of elements to be removed from our list

length(hplc_pre_append) #length is now 418
length(hplc_names_subset_2020_df[,1]) #418

#adjust dataframe to account for new row reference numbers
hplc_names_subset_2020_df$row_number = seq(1:length(hplc_names_subset_2020_df[,1])) 

#use the stri_replace_asll_regex function to replace multiple strings
hplc_names_subset_2020_df$hplc_sample_names_subset = 
  stri_replace_all_regex(hplc_names_subset_2020_df$hplc_sample_names_subset,
                         pattern=c('Jen', 
                                   'Frenem', 'French',
                                   'Does','Doe',
                                   'Pre','Post','post',
                                   ' '),
                         replacement=c('JF', 
                                       'FC', 'FC',
                                       'DF','DF',
                                       'Pr','Po','Po',
                                       '_'),
                         vectorize=FALSE)

#     4c) update dataframe with species and treatment information

#Performing according to: https://stackoverflow.com/questions/25372082/create-column-based-on-presence-of-string-pattern-and-ifelse

hplc_names_subset_2020_updated_df = hplc_names_subset_2020_df #new variable
length(hplc_names_subset_2020_updated_df$hplc_sample_names_subset) #418; length of just the subset is 4 for some reason

#use grepl to find instances of sites and make changes accordingly # | is "or"
hplc_names_subset_2020_updated_df$species <- ifelse(grepl("DF|FC",hplc_names_subset_2020_updated_df$hplc_sample_names_subset),
                                                    'Green','White') #basically if it's not a green site, it's a white site

#do same thing for pre- and post-samples
hplc_names_subset_2020_updated_df$sample_type <- ifelse(grepl("Pr",hplc_names_subset_2020_updated_df$hplc_sample_names_subset),
                                                        'Pre','Post') #if it's not a pre-sample, it's a post-sample

#Manual update of all the hplc_sample_names_subset that are still missing information
hplc_names_subset_2020_final_df = hplc_names_subset_2020_updated_df

#add 20_ in front of the samples that don't have it
# Subset the relevant column of the dataframe
subset_col <- hplc_names_subset_2020_final_df[c(1:64,163,165,175), 3]

# Use the paste() function to add "19_" in front of each value
new_col <- paste("20_", subset_col, sep = "")

# Assign the new column to the original dataframe
hplc_names_subset_2020_final_df[c(1:64,163,165,175), 3] <- new_col

#manually assign trts; note I made adjustments to this for the unpreprocessed data. I probably should have commented this out and copied it but oops
#If I use the preprocessed data I will need to go back in and adjust
hplc_names_subset_2020_final_df[1,3] = "20_JF_Po_2-2_Meja0"
hplc_names_subset_2020_final_df[2,3] = "20_JF_Po_2-1_EAB00"
hplc_names_subset_2020_final_df[3,3] = "20_FC_Po_4-11_EAB00"
hplc_names_subset_2020_final_df[4,3] = "20_FC_Po_1-11_Meja0"
hplc_names_subset_2020_final_df[5,3] = "20_EF_Po_2-12_Meja0"
hplc_names_subset_2020_final_df[6,3] = "20_DF_Pr_2-2_EAB00"
hplc_names_subset_2020_final_df[7,3] = "20_EF_Po_2-1_Meja0"
hplc_names_subset_2020_final_df[8,3] = "20_DF_Pr_4-2_Meja0"
hplc_names_subset_2020_final_df[9,3] = "20_DF_Po_3-10_EAB00"
hplc_names_subset_2020_final_df[10,3] = "20_DF_Po_1-2_Meja0"
hplc_names_subset_2020_final_df[11,3] = "20_EF_Pr_3-3_Meja0"
hplc_names_subset_2020_final_df[12,3] = "20_EF_Pr-2-12_Meja0"
hplc_names_subset_2020_final_df[13,3] = "20_JF_Pr_1-1_Contr"
hplc_names_subset_2020_final_df[14,3] = "20_JF_Pr_4-2_EAB00"
hplc_names_subset_2020_final_df[15,3] = "20_FC_Pr_3-10_Meja0"
hplc_names_subset_2020_final_df[16,3] = "20_FC_Pr_1-9_Contr"
hplc_names_subset_2020_final_df[17,3] = "20_EF_Po_3-11_Contr"
hplc_names_subset_2020_final_df[18,3] = "20_EF_Pr_4-1_Contr"
hplc_names_subset_2020_final_df[19,3] = "20_FC_Pr_1-12_EAB00"
hplc_names_subset_2020_final_df[20,3] = "20_DF_pre_1-10_Contr"
hplc_names_subset_2020_final_df[21,3] = "20_FC_Po_3-10_Contr"
hplc_names_subset_2020_final_df[22,3] = "20_JF_Pr_1-11_Meja0"
hplc_names_subset_2020_final_df[23,3] = "20_EF_Pr_2-2_EAB00"
hplc_names_subset_2020_final_df[24,3] = "20_JF_Po_4-11_Meja0"
hplc_names_subset_2020_final_df[25,3] = "20_DF_Pr_3-2_Contr"
hplc_names_subset_2020_final_df[26,3] = "20_JF_Pr_2-12_EAB00"
hplc_names_subset_2020_final_df[27,3] = "20_DF_Po_2-11_Contr"
hplc_names_subset_2020_final_df[28,3] = "20_JF_Po_4-1_Contr"
hplc_names_subset_2020_final_df[29,3] = "20_FC_Pr_2-3_Contr"
hplc_names_subset_2020_final_df[30,3] = "20_DF_Po_1-1_EAB00"
hplc_names_subset_2020_final_df[31,3] = "20_FC_Po_1_Unkno" #unknown
hplc_names_subset_2020_final_df[32,3] = "20_EF_Po_4-3_EAB00"
hplc_names_subset_2020_final_df[33,3] = "20_EF_Po_8_Unkno" #unknown
hplc_names_subset_2020_final_df[34,3] = "20_FC_Po_2-3_Contr"
hplc_names_subset_2020_final_df[35,3] = "20_DF_Pr_3-11_Meja0"
hplc_names_subset_2020_final_df[36,3] = "20_EF_Pr_3-10_Meja0"
hplc_names_subset_2020_final_df[37,3] = "20_EF_Po_2-10_EAB00"
hplc_names_subset_2020_final_df[38,3] = "20_DF_Po_4-12_EAB00"
hplc_names_subset_2020_final_df[39,3] = "20_DF_Po_2-12_EAB00"
hplc_names_subset_2020_final_df[40,3] = "20_FC_Pr_4-10_Contr"
hplc_names_subset_2020_final_df[41,3] = "20_JF_Po_1-2_Meja0"
hplc_names_subset_2020_final_df[42,3] = "20_EF_Pr_4-2_Meja0"
hplc_names_subset_2020_final_df[43,3] = "20_DF_Pr_4-10_Meja0"
hplc_names_subset_2020_final_df[44,3] = "20_JF_Pr_2-1_EAB00"
hplc_names_subset_2020_final_df[45,3] = "20_FC_Po_1-10_Contr"
hplc_names_subset_2020_final_df[46,3] = "20_JF_Po_3-12_EAB00"
hplc_names_subset_2020_final_df[47,3] = "20_JF_Pr_4-10_EAB00"
hplc_names_subset_2020_final_df[48,3] = "20_FC_Pr_3-5_EAB00"
hplc_names_subset_2020_final_df[49,3] = "20_DF_Po_3-12_Contr"
hplc_names_subset_2020_final_df[50,3] = "20_DF_Pr_2-3_Meja0"
hplc_names_subset_2020_final_df[51,3] = "20_EF_Pr_1-10_Meja0"
hplc_names_subset_2020_final_df[52,3] = "20_JF_Pr_4-12_Contr"
hplc_names_subset_2020_final_df[53,3] = "20_DF_Po_3-1_EAB00"
hplc_names_subset_2020_final_df[54,3] = "20_EF_Po_2-3_Meja0"
hplc_names_subset_2020_final_df[55,3] = "20_FC_Pr_1-3_Meja0"
hplc_names_subset_2020_final_df[56,3] = "20_JF_Po_1-11_Meja0"
hplc_names_subset_2020_final_df[57,3] = "20_JF_Pr_2-11_Contr"
hplc_names_subset_2020_final_df[58,3] = "20_FC_Po_4-12_Meja0"
hplc_names_subset_2020_final_df[59,3] = "20_JF_Po_8_Unkno" #unknown
hplc_names_subset_2020_final_df[60,3] = "20_EF_Po_1_Unkno" #unknown
hplc_names_subset_2020_final_df[61,3] = "20_EF_Pr_4-3_EAB00"
hplc_names_subset_2020_final_df[62,3] = "20_DF_Pr_4-11_Contr"
hplc_names_subset_2020_final_df[63,3] = "20_FC_Pr_4-9_EAB00"
hplc_names_subset_2020_final_df[64,3] = "20_FC_Po_3-3_Meja0"
hplc_names_subset_2020_final_df[163,3] = "20_JF_Po_8-10_Unkno" #unknown
hplc_names_subset_2020_final_df[165,3] = "20_JF_Po_8-1_Unkno" #unknown; maybe the 8 was a 3?
hplc_names_subset_2020_final_df[175,3] = "20_DF_Pr_8-12_Unkno" #unknown; maybe the 8 was a 3?

#create a treatment column from the last 5 characters in the sample ID
hplc_names_subset_2020_final_df$treatment <- substr(hplc_names_subset_2020_final_df$hplc_sample_names_subset, 
                                                    nchar(hplc_names_subset_2020_final_df$hplc_sample_names_subset)-4, 
                                                    nchar(hplc_names_subset_2020_final_df$hplc_sample_names_subset))

#replace treatment codes with actual treatment names
hplc_names_subset_2020_final_df$treatment = 
  stri_replace_all_regex(hplc_names_subset_2020_final_df$treatment,
                         pattern=c('Meja0', 'MeJA0', 
                                   'Contr',
                                   'EAB00'),
                         replacement=c('Methyl jasmonate', 'Methyl jasmonate', 
                                       'Control',
                                       'Emerald ash borer'),
                         vectorize=FALSE)

#extract tree number by taking _,-, and numeric characters from our sample names subset
hplc_names_subset_2020_final_df$tree_id <- gsub("[^_\\x2D0-9]", "", #x2D is the hex code for "-"  
                                                substr(hplc_names_subset_2020_final_df$hplc_sample_names_subset, 
                                                       4, 13))

#not great looking but closer, lets pull the last 4 characters
hplc_names_subset_2020_final_df$tree_id_new <- substr(hplc_names_subset_2020_final_df$tree_id, 
                                                      nchar(hplc_names_subset_2020_final_df$tree_id)-3, 
                                                      nchar(hplc_names_subset_2020_final_df$tree_id))

#lets replace "D", "0", and "_" with ""
hplc_names_subset_2020_final_df$tree_id_new = 
  stri_replace_all_regex(hplc_names_subset_2020_final_df$tree_id_new,
                         pattern=c('D','_', '0' 
                         ),
                         replacement=c('', '', ''
                         ),
                         vectorize=FALSE)

#insert a hyphen between first and second character
hplc_names_subset_2020_final_df$tree_id_new2 <- paste0(substr(hplc_names_subset_2020_final_df$tree_id_new, 1, 1),
                                                       "-", 
                                                       substr(hplc_names_subset_2020_final_df$tree_id_new, 2, 
                                                              nchar(hplc_names_subset_2020_final_df$tree_id_new)))

# create a new column "dbh_cat" based on "tree_id_new2"
#set tree size category according to the last two digits of the string we extracted and cleaned:
#NOTE: 9May23: NAs are sample ids with missing data to be addressed at a later date
#tree number 1-3 = 3-6 cm DBH
#            4-6 = 6.01-9 cm DBH
#            7-9- = 9.01-12 cm DBH
#            10-12 = 12.01-15 cm DBH
hplc_names_subset_2020_final_df$dbh_cat <- ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2020_final_df$tree_id_new2), 2, 3) %in% c("1", "2", "3"), "3-6",
                                                  ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2020_final_df$tree_id_new2), 2, 3) %in% c("4", "5", "6"), "6.01-9",
                                                         ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2020_final_df$tree_id_new2), 2, 3) %in% c("7", "8", "9"), "9.01-12",
                                                                ifelse(substr(gsub("[^0-9]", "", hplc_names_subset_2020_final_df$tree_id_new2), 2, 3) %in% c("10", "11", "12"), "12.01-15", NA))))

# print the updated columns
print(hplc_names_subset_2020_final_df[,c(1,3,5:7,11)]) #looks good; prior to 18May started wiht 3. need 1 for metadata attach below

#assign year
hplc_names_subset_2020_final_df$sample_year <- ifelse(grepl("19", hplc_names_subset_2020_final_df$hplc_sample_names_subset), "2020", #there are no 2019 samples in this subset
                                                      ifelse(grepl("20", hplc_names_subset_2020_final_df$hplc_sample_names_subset), "2020", "2020")) #last value was NA, but because these are all 2020 samples, I wanted to make sure everything was assigned

#assign site
hplc_names_subset_2020_final_df$site <- ifelse(grepl("DF", hplc_names_subset_2020_final_df$hplc_sample_names_subset), "Doe Farm",
                                               ifelse(grepl("EF", hplc_names_subset_2020_final_df$hplc_sample_names_subset), "East Foss Farm",
                                                      ifelse(grepl("JF", hplc_names_subset_2020_final_df$hplc_sample_names_subset), "Jennings Forest",
                                                             ifelse(grepl("FC", hplc_names_subset_2020_final_df$hplc_sample_names_subset), "French Conservation Area", 
                                                                    NA))))
#final dataframe
hplc_2020_df = hplc_names_subset_2020_final_df[,c(1,3,5:7,11:13)]
length(hplc_2020_df$hplc_sample_names_subset) #418; length() without is 7 which is the number of columns?

#---alignment---#

#I spoke with Casey Philbin (UNR Chem Dept) about aligning peaks and he seemed to think aligning by species would be a good idea
#I'd like to try and take this further and:
# 1) align all pre-samples according to tree species
# 2) align post-samples by treatment

#12May23- aligning white ash by itself causes R to crash, attempting to align by site instead

#-subset by site-#

#create subsets
doe = subset(hplc_2020_df,site =="Doe Farm")
eas = subset(hplc_2020_df,site =="East Foss Farm")
fre =subset(hplc_2020_df,site =="French Conservation Area")
jen = subset(hplc_2020_df,site =="Jennings Forest")

#hplc data length = 417; 418
rt_2020_doe = hplc_pre_append 
rt_2020_eas = hplc_pre_append 
rt_2020_fre = hplc_pre_append 
rt_2020_jen = hplc_pre_append 

#get row ids for green and white subsets / create variable
doe$row_number = as.numeric(rownames(doe)) #length(doe$row_number) = 104
eas$row_number = as.numeric(rownames(eas)) #length(eas$row_number) = 104
fre$row_number = as.numeric(rownames(fre)) #length(fre$row_number) = 105
jen$row_number = as.numeric(rownames(jen)) #length(jen$row_number) = 102

#correct row number by subtracting 208 from them all
doe$row_number_new = doe$row_number-208
eas$row_number_new = eas$row_number-208
fre$row_number_new = fre$row_number-208
jen$row_number_new = jen$row_number-208

#remove non-doe rows
rt_2020_doe[c(eas$row_number_new,fre$row_number_new,jen$row_number_new)] = NULL #list of elements to be removed from our list
length(rt_2020_doe) #length is now 106 (not sure why 2 more than length above)

#check names to confirm we only have doe samples
rt_2020_doe_sample_names = sapply(rt_2020_doe, function(x) attr(x,"sample_name")) #confirmed

#remove non-eas rows
rt_2020_eas[c(doe$row_number_new,fre$row_number_new,jen$row_number_new)] = NULL #list of elements to be removed from our list
length(rt_2020_eas) #length is now 106 (not sure why 2 more than length above); 107 21Jun23

#check names to confirm we only have doe samples
rt_2020_eas_sample_names = sapply(rt_2020_eas, function(x) attr(x,"sample_name")) #confirmed

#remove non-fre rows
rt_2020_fre[c(eas$row_number_new,doe$row_number_new,jen$row_number_new)] = NULL #list of elements to be removed from our list
length(rt_2020_fre) #length is now 107 (not sure why 2 more than length above)

#check names to confirm we only have doe samples
rt_2020_fre_sample_names = sapply(rt_2020_fre, function(x) attr(x,"sample_name")) #confirmed

#remove non-jen rows
rt_2020_jen[c(eas$row_number_new,fre$row_number_new,doe$row_number_new)] = NULL #list of elements to be removed from our list
length(rt_2020_jen) #length is now 104 (not sure why 2 more than length above)

#check names to confirm we only have doe samples
rt_2020_jen_sample_names = sapply(rt_2020_jen, function(x) attr(x,"sample_name")) #confirmed

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
                     dim2=seq(200,340,1),
                     remove.time.baseline = TRUE,
                     cl = cl) #I don't know why I removed this previously
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
    progress_bar = TRUE, 
    cl = cl)

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
                                           max.iter=2000,
                                           cl = cl) #works

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

#looking at this current alignment makes it seem like there are quite a few peaks that have been artificially split, for example
#looking at V46 and V47 we see most of the absorbance in V46 and when there is a 0 in V47 there is an absorbance value there, and when
#there are absorbance values in V46 they are 0 in V47. Looks like this has happened a number of places
#I am saving a 19.5May23 file and I'm going to go back up and see if I can make some of the alignments a bit more stringent to see if
#we can lump this together more. All the old alignments will be in 19May23 files and earlier. See egh2 above

#2020 - French Conservation Area - Green

#---pre-process data---#

#NOTE: There are duplicates in some of these subsets. Make sure that those that are removed have been retained or returned to 
#the correct position within their dataframe

#remove duplicates from metadata
hplc_2020_df_fre_nodup = fre[!duplicated(fre$hplc_sample_names_subset),] #no duplicates
fre_meta_df = cbind.data.frame(fre$hplc_names_subset,fre$hplc_sample_names_subset)

#before we align, there are some duplicate and extra samples that were somehow retained. Lets remove them.
#There were some remove duplicates from processed data
rt_2020_fre_nodup <- rt_2020_fre[!duplicated(rt_2020_fre)] #no duplicates
length(rt_2020_fre_nodup) #107
length(hplc_2020_df_fre_nodup$site) #103 - where is the difference here?

#check names of french samples
rt_2020_fre_sample_names = sapply(rt_2020_fre, function(x) attr(x,"sample_name")) #confirmed
fre_sample_df = as.data.frame(rt_2020_fre_sample_names)

#remove JF sample from list of lists
rt_2020_fre_nodup[60] = NULL 
length(rt_2020_fre_nodup) #106


#add the missing Comma088 to the metadata
hplc_2020_df_fre_nodup = hplc_2020_df_fre_nodup[-c(9:10)] #remove these extraneous columns that idk where they came from

colnames(hplc_2020_df_fre_nodup) <- c("hplc_name","sample_name","tree_species","sample_type","treatment","tree_size","year","site")


#add row to end of meta dataframe
hplc_2020_df_fre_nodup = hplc_2020_df_fre_nodup %>% add_row(hplc_name="ASCIIData_Run#2_Comma088",
                                                            sample_name="20_FC_Po_05_0_Contr",
                                                            tree_species="Green",
                                                            sample_type="Post",
                                                            treatment="Control",
                                                            tree_size="6.01-9",
                                                            year="2020",
                                                            site="French Conservation Area")

hplc_2020_df_fre_nodup = hplc_2020_df_fre_nodup %>% add_row(hplc_name="ASCIIData_Run#2_Comma097",
                                                            sample_name="20_FC_Po_08_0_Contr",
                                                            tree_species="Green",
                                                            sample_type="Post",
                                                            treatment="Control",
                                                            tree_size="9.01-12",
                                                            year="2020",
                                                            site="French Conservation Area")

hplc_2020_df_fre_nodup = hplc_2020_df_fre_nodup %>% add_row(hplc_name="ASCIIData_Run#2_Comma216",
                                                            sample_name="20_FC_Po_02_1_EAB00",
                                                            tree_species="Green",
                                                            sample_type="Post",
                                                            treatment="Emerald Ash Borer",
                                                            tree_size="3-6",
                                                            year="2020",
                                                            site="French Conservation Area")

#this fixed the missing rows
length(hplc_2020_df_fre_nodup$site) #106 this is equal to the list of lists now
length(rt_2020_fre_nodup) #106

#according to the ?find_peaks help
#It is recommended to do pre-processing using the preprocess function before peak detection. Overly high chromatographic resolution 
#can sometimes lead to peaks being split into multiple segments. In this case, it is recommended to increase the smooth_window or 
#reduce the resolution along the time axis by adjusting the dim2 argument before proceeding with further analyses.

#I originally thought I wanted a very fine time resolution but perhaps I do not, I am going to try and pre-process at 0.05 and then 0.1 to see how this
#alters the number of peaks found. I want to see ~70-90 as the HPLC output originally suggested
fre_pre = preprocess(rt_2020_fre_nodup, 
                     dim1=seq(0,60,.05),
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE,
                     cl = cl)#, #I don't know why I removed this previously

#align data
fre_pre_05_corrected_75 = 
  correct_rt(
    fre_pre, 
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
    progress_bar = TRUE,
    cl = cl)

#now I am wondering if I should be using ptw?

###---Peak finding---###
#adjusting the smooth window changes the number of peaks output/detected below
#0.75 = 55 peaks;0.05 = 60 peaks; 0.04 = 72 peaks (seems close to HPLC determination); 0.035 = 80 peaks (seems close to HPLC determination)  
#0.025 = 94 peaks; 0.02 = 103 peaks; 0.015 = 121 peaks

#use this subset for the report
fre_pre_corrected_75_pks_egh <- get_peaks(fre_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.035, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000,
                                          cl = cl) #works


#plot peak finding
plot(fre_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035) looks like garbage but is it actually?

###---Creating a peak table---###
fre_pre_corrected_75_egh_peak <- get_peaktable(fre_pre_corrected_75_pks_egh, response = 'height') #

#subset the metadata dataframe for fre
fre_meta = subset(hplc_2020_df_fre_nodup, site=="French Conservation Area")
length(fre_meta$site) #106

#fix messed up treatment/ID
fre_meta[72,2] = "20-FC-Pr-03-5-EAB00"
fre_meta[72,5] = "Emerald Ash Borer"

#attach metadata
#fre_pk_tab <- attach_metadata(peak_table = fre_pre_corrected_75_egh_peak, 
#                              metadata = fre_meta, 
#                              column = "hplc_name") 

#Error in attach_metadata(peak_table = fre_pre_corrected_75_egh_peak, metadata = fre_meta,  : 
#Sample names must be unique. Please check column ‘hplc_name’ for duplicates.

#correct data for analysis
table(fre_meta$sample_name) #2 of: 20_FC_Po_02_1_EAB00 - bizarre that both this and the following were missing before and are now duplicats?
table(fre_meta$hplc_name) #no dups so I guess the sample above was ran twice?

table(rownames(fre_pre_corrected_75_egh_peak$tab))

peak_names = rownames(fre_pre_corrected_75_egh_peak$tab)
peak_names_alp = peak_names[order(peak_names)]
fre_meta_alp = fre_meta[order(fre_meta$hplc_name),]

fre_check = cbind(peak = peak_names_alp, meta = fre_meta_alp[,1]) #23 and 24 for meta are dups

#remove row 23 from meta but that will make it unequal?
fre_meta = fre_meta[-23,] #remov 20_FC_Po_05_Control
fre_meta = fre_meta[-53,] #remove 20_FC_Po_02_1EAB00 dup
length(fre_meta$sample_name) #104
fre_meta_alp = fre_meta[order(fre_meta$hplc_name),]

length(fre_pre_corrected_75_egh_peak$tab$V1) #106

#212 is missing in meta and present in the peak table, so add in
fre_meta = fre_meta %>% add_row(hplc_name="ASCIIData_Run#3_Comma212",
                                sample_name="20_FC_Po_04_7_Meja0",
                                tree_species="Green",
                                sample_type="Post",
                                treatment="Methyl Jasmonate",
                                tree_size="9.01-12",
                                year="2020",
                                site="French Conservation Area")

length(fre_meta$sample_name) #105
#add something to meta that is missing?

meta_hplc_names = table(fre_meta$hplc_name)
peak_hplc_names = table(rownames(fre_pre_corrected_75_egh_peak$tab))

#data frame with a list of ids and the number of times they occurred.
n_occur <- data.frame(table(fre_meta$hplc_name))

#which ids occurred more than once
n_occur[n_occur$Freq > 1,] #in meta: ASCIIData_Run#2_Comma216 occurs twice

#data frame with a list of ids and the number of times they occurred.
n_occur <- data.frame(table(rownames(fre_pre_corrected_75_egh_peak$tab)))

#which ids occurred more than once
n_occur[n_occur$Freq > 1,] #in meta: Nothing occurs twice

meta_hplc_names = fre_meta$hplc_name
meta_hplc_names = as.data.frame(meta_hplc_names)
peak_hplc_names = rownames(fre_pre_corrected_75_egh_peak$tab)

peak_hplc_names = peak_hplc_names[order(peak_hplc_names)]
meta_hplc_names = meta_hplc_names %>% add_row(meta_hplc_names="test"
                                              )
fre_check = cbind(peak = peak_hplc_names, meta = meta_hplc_names) #23 and 24 for meta are dups

#values not in meta_hplc_names (see: https://stackoverflow.com/questions/13774773/check-whether-values-in-one-data-frame-column-exist-in-a-second-data-frame)
fre_check$peak[!fre_check$peak %in% fre_check$meta_hplc_names] #"ASCIIData_Run#3_Comma216"
fre_check$meta_hplc_names[!fre_check$meta_hplc_names %in% fre_check$peak] #test

#add 20-FC-Po-03-8-EAB00 to fre_meta
fre_meta = fre_meta %>% add_row(hplc_name="ASCIIData_Run#3_Comma216",
                                sample_name="20_FC_Po_03_8_EAB00",
                                tree_species="Green",
                                sample_type="Post",
                                treatment="Emerald Ash Borer",
                                tree_size="9.01-12",
                                year="2020",
                                site="French Conservation Area")

length(fre_meta$sample_name)

#attach metadata - finally works
fre_pk_tab <- attach_metadata(peak_table = fre_pre_corrected_75_egh_peak, 
                              metadata = fre_meta, 
                              column = "hplc_name")


#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
fre_pk_tab_df = as.data.frame(do.call(cbind, fre_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
fre_pk_tab_rt_df = as.data.frame(do.call(cbind, fre_pk_tab[2]))

#2020- East Foss Farm - White

#before we align, there are some duplicate and extra samples that were somehow retained. Lets remove them.
#There were some remove duplicates from processed data
rt_2020_eas_nodup <- rt_2020_eas[!duplicated(rt_2020_eas)]
length(rt_2020_eas) #107
length(rt_2020_eas_nodup) #107

eas_dup_check = sapply(rt_2020_eas_nodup, function(x) attr(x,"sample_name"))
eas_dup_check = as.data.frame(eas_dup_check)
eas_dup_check$row_number = seq(1:length(eas_dup_check[,1])) 

rt_2020_eas_nodup[c(26,53)] = NULL 
eas_dup_null_check = sapply(rt_2020_eas_nodup, function(x) attr(x,"sample_name"))
length(eas_dup_null_check) #105

#remove duplicates from metadata
hplc_2020_df_eas_nodup = eas[!duplicated(eas$hplc_sample_names_subset),]

#change column name to be consistent in list and dataframe
colnames(hplc_2020_df_eas_nodup) <- c("hplc_name","sample_name","tree_species","sample_type","treatment","tree_size","year","site")
 length(hplc_2020_df_eas_nodup$sample_name) #same length as eas_pre_nodup #105

#pre-process
#matplot(sapply(rt_2020_eas_nodup,function(x)x[,"280.18"]),type='l')
#legend("topright", legend="0.05", bty = "n")

rt_2020_eas_nodup_end = sapply(rt_2020_eas_nodup, function(x) attr(x,"end_time"))

eas_pre = preprocess(rt_2020_eas_nodup, #causes an error, but okay I think: https://github.com/facebookexperimental/Robyn/issues/576
                     dim1=seq(0,59,.05), #there was an error because one of the times ended at 59.648 not 60
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE, 
                     cl = cl) 

#align data
eas_pre_05_corrected_75 = 
  correct_rt(
    eas_pre, 
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
    progress_bar = TRUE,
    cl = cl)

###---Peak finding---###
#adjusting the smooth window changes the number of peaks output/detected below
#using 0.035 here as it's what I used above = 73 peaks (seems close to HPLC determination)  


#use this subset for the report
eas_pre_corrected_75_pks_egh <- get_peaks(eas_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.035, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000,
                                          cl = cl) #works

#plot peak finding
plot(eas_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035) looks like garbage but is it actually?

###---Creating a peak table---###
eas_pre_corrected_75_egh_peak <- get_peaktable(eas_pre_corrected_75_pks_egh, response = 'height') #

#subset the metadata dataframe for eas
eas_meta = subset(hplc_2020_df_eas_nodup, site=="East Foss Farm")
length(eas_meta$site) #105

#attach metadata
eas_pk_tab <- attach_metadata(peak_table = eas_pre_corrected_75_egh_peak, 
                              metadata = eas_meta[c(1:8)], 
                              column = "hplc_name") #this works! 

#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
eas_pk_tab_df = as.data.frame(do.call(cbind, eas_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
eas_pk_tab_rt_df = as.data.frame(do.call(cbind, eas_pk_tab[2]))

###---Jennings Forest---###

#before we align, there are some duplicate and extra samples that were somehow retained. Lets remove them.
#There were some remove duplicates from processed data
rt_2020_jen_nodup <- rt_2020_jen[!duplicated(rt_2020_jen)]
length(rt_2020_jen) #104
length(rt_2020_jen_nodup) #104

jen_dup_check = sapply(rt_2020_jen_nodup, function(x) attr(x,"sample_name"))
jen_dup_check = as.data.frame(jen_dup_check)
jen_dup_check$row_number = seq(1:length(jen_dup_check[,1])) 
table(jen_dup_check$jen_dup_check) #2 20_JF_Pr_04_8_Contr

rt_2020_jen_nodup[64] = NULL #remove first instance of duplicate
length(rt_2020_jen_nodup) #103

jen_dup_null_check = sapply(rt_2020_jen_nodup, function(x) attr(x,"sample_name"))
jen_dup_null_check2 = names(rt_2020_jen_nodup)
length(jen_dup_null_check) #103
length(jen$hplc_sample_names_subset) #102

#remove duplicates from metadata
hplc_2020_df_jen_nodup = jen[!duplicated(jen$hplc_sample_names_subset),]
length(hplc_2020_df_jen_nodup$hplc_names_subset) #100

table(hplc_2020_df_jen_nodup$hplc_sample_names_subset)
length(hplc_2020_df_jen_nodup$hplc_sample_names_subset) #100, which are missing

#values not in jen_dup_null_check (see: https://stackoverflow.com/questions/13774773/check-whether-values-in-one-data-frame-column-exist-in-a-second-data-frame)
jen_dup_null_check[!jen_dup_null_check %in% hplc_2020_df_jen_nodup$hplc_sample_names_subset]

#this is a wild number of things not in the other? Until line 942 these are all just named differently
#ASCIIData_Run#2_Comma001 ASCIIData_Run#2_Comma002 ASCIIData_Run#2_Comma013 ASCIIData_Run#2_Comma014 ASCIIData_Run#2_Comma022 
#"Jen Post 2-2F"          "Jen Post 2-1F"           "Jen Pre 1-1F"           "Jen Pre 4-2F"          "Jen Pre 1-11F" 
#ASCIIData_Run#2_Comma024 ASCIIData_Run#2_Comma027 ASCIIData_Run#2_Comma029 ASCIIData_Run#2_Comma042 ASCIIData_Run#2_Comma045 
#"Jen Post 4-11F"          "Jen Pre 2-12F"          "Jen Post 4-1F"          "Jen Post 1-2F"           "Jen Pre 2-1F" 
#ASCIIData_Run#2_Comma047 ASCIIData_Run#2_Comma048 ASCIIData_Run#2_Comma054 ASCIIData_Run#2_Comma058 ASCIIData_Run#2_Comma059 
#"Jen Post 3-12F"          "Jen Pre 4-10F"          "Jen Pre 4-12F"         "Jen Post 1-11F"          "Jen Pre 2-11F" 
#ASCIIData_Run#2_Comma061 ASCIIData_Run#2_Comma097 ASCIIData_Run#2_Comma169 ASCIIData_Run#2_Comma171 ASCIIData_Run#3_Comma018 
#"Jen Post 8F"    "20_Fc_Po_08_0_Contr"         "Jen Post 8-10F"          "Jen Post 8-1F"    "20-Jf-Pr-02-5-MeJA0" 

#ASCIIData_Run#2_Comma169 ASCIIData_Run#2_Comma171
#"Jen Post 8-10F"          "Jen Post 8-1F"    

#ASCIIData_Run#2_Comma097 ASCIIData_Run#3_Comma018 
#"20_Fc_Po_08_0_Contr".   "20-Jf-Pr-02-5-MeJA0"

#lets try and look at the actual hplc names
jen_dup_null_check2[!jen_dup_null_check2 %in% hplc_2020_df_jen_nodup$hplc_names_subset]
# "ASCIIData_Run#2_Comma097" "ASCIIData_Run#2_Comma175" "ASCIIData_Run#3_Comma018" "ASCIIData_Run#3_Comma055"

#change column name to be consistent in list and dataframe
colnames(hplc_2020_df_jen_nodup) <- c("hplc_name","sample_name","tree_species","sample_type","treatment","tree_size","year","site")
jen_meta = hplc_2020_df_jen_nodup[-c(9:10)]

jen_meta = jen_meta %>% add_row(hplc_name="ASCIIData_Run#2_Comma175",
                                sample_name="20_JF_Po_4_12_Contr",
                                tree_species="White",
                                sample_type="Post",
                                treatment="Control",
                                tree_size="12.01-15",
                                year="2020",
                                site="Jennings Forest")

jen_meta = jen_meta %>% add_row(hplc_name="ASCIIData_Run#3_Comma018",
                                sample_name="20_JF_Pr_02_5_MeJA0",
                                tree_species="White",
                                sample_type="Pre",
                                treatment="Methyl Jasmonate",
                                tree_size="6.01-9",
                                year="2020",
                                site="Jennings Forest")

jen_meta = jen_meta %>% add_row(hplc_name="ASCIIData_Run#3_Comma055",
                                sample_name="20_JF_Pr_04_8_Contr",
                                tree_species="White",
                                sample_type="Pre",
                                treatment="Control",
                                tree_size="9.01-12",
                                year="2020",
                                site="Jennings Forest")

length(jen_meta$hplc_name) #103

#remove French sample
rt_2020_jen_nodup[23] = NULL
length(rt_2020_jen_nodup) #102

#data frame with a list of ids and the number of times they occurred.
n_occur <- data.frame(table(jen_meta$hplc_name))

#which ids occurred more than once
n_occur[n_occur$Freq > 1,]

#data frame with a list of ids and the number of times they occurred.
#n_occur <- data.frame(table(rownames(jen_pre_corrected_75_egh_peak$tab)))

#which ids occurred more than once
n_occur[n_occur$Freq > 1,] #in meta: Nothing occurs twice

length(jen_meta$hplc_name) #same length as jen_pre_nodup #103
length(rt_2020_jen_nodup) #102

jen_dup_null_check3 = names(rt_2020_jen_nodup)
jen_meta$hplc_name[!jen_meta$hplc_name %in% jen_dup_null_check3]
#[1] "ASCIIData_Run#3_Comma017" #for sake of ease (at least now, just drop this from the data to make it easier)

jen_meta = jen_meta[-62,]

jen_meta$hplc_name[!jen_meta$hplc_name %in% jen_dup_null_check3] #none
jen_dup_null_check3[!jen_dup_null_check3 %in% jen_meta$hplc_name] #none 

length(jen_meta$hplc_name) #102
length(jen_dup_null_check3) #102

#pre-process
jen_pre = preprocess(rt_2020_jen_nodup,
                     dim1=seq(0,60,.05),
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE,
                     cl = cl) 

#align data
jen_pre_05_corrected_75 = 
  correct_rt(
    jen_pre, 
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
    progress_bar = TRUE,
    cl = cl)

###---Peak finding---###
#adjusting the smooth window changes the number of peaks output/detected below
#using 0.035 here as it's what I used above = 73 peaks (seems close to HPLC determination)  

#use this subset for the report
jen_pre_corrected_75_pks_egh <- get_peaks(jen_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.035, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000,
                                          cl = cl) #works

#plot peak finding
plot(jen_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035) looks like garbage but is it actually?

###---Creating a peak table---###
jen_pre_corrected_75_egh_peak <- get_peaktable(jen_pre_corrected_75_pks_egh, response = 'height') #

#attach metadata
jen_pk_tab <- attach_metadata(peak_table = jen_pre_corrected_75_egh_peak, 
                              metadata = jen_meta, 
                              column = "hplc_name") #this works! 

#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
jen_pk_tab_df = as.data.frame(do.call(cbind, jen_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
jen_pk_tab_rt_df = as.data.frame(do.call(cbind, jen_pk_tab[2]))









###---2021---###
#setwd("/Users/toddjoh/Ash Phenolics/2021/Seperate_Comma")

###---read data---###

#for running on my laptop
#run_2021_separate_comma <- read_chroms("/Users/colto/Documents/AshPhenolics/2021/Seperate_Comma", 
#                                       parser="chromconverter",
#                                       read_metadata = TRUE,
#                                       format_in = "shimadzu_dad", cl = cl)

#Separate comma format (note that the file folder has an typo, fix at some point)
run_2021_separate_comma <- read_chroms("/Users/colto/Documents/AshPhenolics/2021/Seperate_Comma", 
                                       parser="chromconverter",
                                       read_metadata = TRUE,
                                       format_in = "shimadzu_dad", cl = cl)

###---Pre-process data---###

#exclude standards from peak alignment
#first get sample names to know which to remove
hplc_2021_sample_names = sapply(run_2021_separate_comma, function(x) attr(x,"sample_name")) #this pulls the standards, the sample names are in another column for some reason
hplc_2021_sample_ids = sapply(run_2021_separate_comma, function(x) attr(x,"sample_id"))

#list of standards to be removed

#ASCIIData_2021_Comma "Luteolin_0.5"
#ASCIIData_2021_Comma001 "Oleuropein_0.5"
#ASCIIData_2021_Comma002 "Syringin_0.5"
#ASCIIData_2021_Comma003 "Verbascoside_0.5"
#ASCIIData_2021_Comma004 "Pinoresinol_0.5"
#ASCIIData_2021_Comma005 "Ligustrocide_0.5"
#ASCIIData_2021_Comma006 "Chlorogenic Acid_0.5"
#ASCIIData_2021_Comma007 "Vanillic Acid_0.5"
#ASCIIData_2021_Comma008 "Ferulic Acid_0.5"
#ASCIIData_2021_Comma009 "Coumaric Acid_0.5"
#ASCIIData_2021_Comma010 "Gallic Acid_0.5"
#ASCIIData_2021_Comma011 "Apigenin_0.5"
#ASCIIData_2021_Comma012 "Apigenin_0.2"
#ASCIIData_2021_Comma013 "Apigenin_0.01"
#ASCIIData_2021_Comma014 "Standard #1_0.25"
#ASCIIData_2021_Comma015 "Standard #1_0.125"
#ASCIIData_2021_Comma016 "Standard #2_0.25"
#ASCIIData_2021_Comma017 "Standard #2_0.125"
#ASCIIData_2021_Comma018 "Syringin_2"
#ASCIIData_2021_Comma019 "Syringin_1"
#ASCIIData_2021_Comma020 "Standard #3_0.25"
#ASCIIData_2021_Comma021 "Ligustrocide_0.125"
#ASCIIData_2021_Comma022 "Check Standard #1"

#ASCIIData_2021_Comma048 "Check Standard #2"
#ASCIIData_2021_Comma074 "Check Standard #3"
#ASCIIData_2021_Comma100 "Check Standard #4"
#ASCIIData_2021_Comma126 "Check Standard #5"
#ASCIIData_2021_Comma152 "Check Standard #6"
#ASCIIData_2021_Comma178 "Check Standard #7"
#ASCIIData_2021_Comma203 "Check Standard #8"
#ASCIIData_2021_Comma229 "Check Standard #9"
#ASCIIData_2021_Comma255 "Check Standard #10"
#ASCIIData_2021_Comma281 "Check Standard #11"
#ASCIIData_2021_Comma307 "Check Standard #12"
#ASCIIData_2021_Comma332 "Check Standard #13"
#ASCIIData_2021_Comma358 "Check Standard #14"
#ASCIIData_2021_Comma384 "Check Standard #15" 
#ASCIIData_2021_Comma410 "Check Standard #16"
#ASCIIData_2021_Comma436 "Check Standard #17"
#ASCIIData_2021_Comma462 "Check Standard #18"
#ASCIIData_2021_Comma488 "Check Standard #19"
#ASCIIData_2021_Comma514 "Check Standard #20"
#ASCIIData_2021_Comma540 "Check Standard #21" 

#make new variable to prevent modification raw data
hplc21_list = run_2021_separate_comma# hplc_pre_append

#pull the names of each list item (these are the file names of the chromatograms)
hplc21_names = names(hplc21_list)

#get a vector with a number for each element in the vector (this is so I know which rows to call to remove standards prior to pre-processing and alignment)
hplc21_length = seq(1:length(hplc21_names))

#get sample names - this pulls the actual sample id as was input into the HPLC
hplc_2021_sample_ids #from above

#For some reason I was having trouble binding this to the rest of the dataframe so i set it as a vector
hplc_2021_sample_ids = as.vector(hplc_2021_sample_ids)

hplc21_length = as.vector(hplc21_length) #need to change this

#bind together in a dataframe
hplc21_names_df = cbind.data.frame(hplc_2021_sample_ids,hplc21_length)

#this allows me to add hplc_names to the dataframe
hplc21_names_df$hplc21_names = hplc21_names

#Remove standards from the list by setting them as NULL; 53 removed
hplc21_list_standards = hplc21_list[c(1:23,49,75,101,127,153,179,204,230,256,282,308,333,359,385,411,437,463,489,515,541)]
hplc21_list[c(1:23,49,75,101,127,153,179,204,230,256,282,308,333,359,385,411,437,463,489,515,541)] = NULL #list of elements to be removed from our list

hplc21_list_updated_sample_id = sapply(hplc21_list, function(x) attr(x,"sample_id"))
hplc21_list_updated_sample_id = cbind(hplc21_list_updated_sample_id,(1:length(hplc21_list_updated_sample_id)))


hplc21_pre_append_standards = hplc21_list_standards 

hplc21_names_subset_standards = names(hplc21_pre_append_standards) #get names of files

hplc21_sample_names_subset_standards = sapply(hplc21_pre_append_standards, function(x) attr(x,"sample_name")) #get sample_names attribute
hplc21_sample_names_subset_standards = as.vector(hplc21_sample_names_subset_standards) #set as a vector

hplc21_length_subset_standards = seq(1:length(hplc21_pre_append_standards)) #assign a number for each row
hplc21_length_subset_standards = as.vector(hplc21_length_subset_standards) #set as a vector

hplc21_names_subset_standards_df = cbind.data.frame(hplc21_names_subset_standards,hplc21_length_subset_standards) #bind together in a dataframe
hplc21_names_subset_standards_df$hplc21_sample_names_subset_standards = hplc21_sample_names_subset_standards #add names

#sub-setting hplc21_list_standards by standards and standard mix to aid in alignment 
standards21_pre_luteolin = hplc21_list_standards[c(1)]
standards21_pre_oleuropein = hplc21_list_standards[c(2)]
standards21_pre_syringin = hplc21_list_standards[c(3,19,20)]
standards21_pre_verbascoside = hplc21_list_standards[c(4)]
standards21_pre_pinoresinol = hplc21_list_standards[c(5)]
standards21_pre_ligustrocide = hplc21_list_standards[c(6,22)]
standards21_pre_chlorogenic_acid = hplc21_list_standards[c(7)]
standards21_pre_vanillic_acid = hplc21_list_standards[c(8)]
standards21_pre_ferulic_acid = hplc21_list_standards[c(9)]
standards21_pre_coumaric_acid = hplc21_list_standards[c(10)]
standards21_pre_gallic_acid = hplc21_list_standards[c(11)]
standards21_pre_standard_mix = hplc21_list_standards[c(15:18,21,23:43)]
standards21_pre_apagenin = hplc21_list_standards[c(12:14)]


#pre-processing subset standards and peak finding
rt_standards21_pre_luteolin = preprocess(standards21_pre_luteolin, 
                                           dim1=seq(0,45,.05),
                                           dim2=seq(210,400,1),
                                           remove.time.baseline = TRUE,
                                           cl=cl)
plot_chroms(rt_standards21_pre_luteolin, lambdas = '280')

rt_standards21_pre_luteolin = c(rt_standards21_pre_luteolin, copy(rt_standards21_pre_luteolin))
standards21_luteolin_peaks = get_peaks(rt_standards21_pre_luteolin, lambdas = c('280'), fit = "egh", amp_thresh = 9000, cl = cl)
standards21_luteolin_peak_table = get_peaktable(standards21_luteolin_peaks, response = 'height')
print(standards21_luteolin_peak_table$pk_meta)



rt_standards21_pre_oleuropein = preprocess(standards21_pre_oleuropein, 
                                         dim1=seq(0,45,.05),
                                         dim2=seq(210,400,1),
                                         remove.time.baseline = TRUE,
                                         cl=cl)
plot_chroms(rt_standards21_pre_oleuropein, lambdas = '280')

rt_standards21_pre_oleuropein = c(rt_standards21_pre_oleuropein, copy(rt_standards21_pre_oleuropein))
standards21_oleuropein_peaks = get_peaks(rt_standards21_pre_oleuropein, lambdas = c('280'), fit = "egh", amp_thresh = 3000, cl = cl)
standards21_oleuropein_peak_table = get_peaktable(standards21_oleuropein_peaks, response = 'height')
print(standards21_oleuropein_peak_table$pk_meta)



rt_standards21_pre_syringin = preprocess(standards21_pre_syringin, 
                                         dim1=seq(0,45,.05),
                                         dim2=seq(210,400,1),
                                         remove.time.baseline = TRUE,
                                         cl=cl)
plot_chroms(rt_standards21_pre_syringin, lambdas = '280') #needs alignment 
rt_standards21_pre_syringin = correct_rt(
  rt_standards21_pre_syringin, 
  what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
  lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
  models = NULL,
  reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
  alg = c("vpdtw"),
  init.coef = c(0, 1, 0),
  n.traces = NULL,
  n.zeros = 0,
  scale = FALSE, 
  plot_it = FALSE,
  penalty = 1, #default 5
  maxshift = 75, #
  verbose = FALSE,
  progress_bar = FALSE,
  cl = cl)
plot_chroms(rt_standards21_pre_syringin, lambdas = '280') #looks better
standards21_syringin_peaks = get_peaks(rt_standards21_pre_syringin, lambdas = c('280'), fit = "egh", amp_thresh = 5000, cl = cl)
standards21_syringin_peak_table = get_peaktable(standards21_syringin_peaks, response = 'height')
print(standards21_syringin_peak_table$pk_meta)



rt_standards21_pre_verbascoside = preprocess(standards21_pre_verbascoside, 
                                         dim1=seq(0,45,.05),
                                         dim2=seq(210,400,1),
                                         remove.time.baseline = TRUE,
                                         cl=cl)
plot_chroms(rt_standards21_pre_verbascoside, lambdas = '280')

rt_standards21_pre_verbascoside = c(rt_standards21_pre_verbascoside, copy(rt_standards21_pre_verbascoside))
standards21_verbascoside_peaks = get_peaks(rt_standards21_pre_verbascoside, lambdas = c('280'), fit = "egh", amp_thresh = 3300, cl = cl)
standards21_verbascoside_peak_table = get_peaktable(standards21_verbascoside_peaks, response = 'height')
print(standards21_verbascoside_peak_table$pk_meta)



rt_standards21_pre_pinoresinol = preprocess(standards21_pre_pinoresinol, 
                                         dim1=seq(0,45,.05),
                                         dim2=seq(210,400,1),
                                         remove.time.baseline = TRUE,
                                         cl=cl)
plot_chroms(rt_standards21_pre_pinoresinol, lambdas = '280')

rt_standards21_pre_pinoresinol = c(rt_standards21_pre_pinoresinol, copy(rt_standards21_pre_pinoresinol))
standards21_pinoresinol_peaks = get_peaks(rt_standards21_pre_pinoresinol, lambdas = c('280'), fit = "egh", amp_thresh = 9000, cl = cl)
standards21_pinoresinol_peak_table = get_peaktable(standards21_pinoresinol_peaks, response = 'height')
print(standards21_pinoresinol_peak_table$pk_meta)



rt_standards21_pre_ligustrocide = preprocess(standards21_pre_ligustrocide, 
                                         dim1=seq(0,45,.05),
                                         dim2=seq(210,400,1),
                                         remove.time.baseline = TRUE,
                                         cl=cl)
plot_chroms(rt_standards21_pre_ligustrocide, lambdas = '280')
rt_standards21_pre_ligustrocide = correct_rt(
  rt_standards21_pre_ligustrocide, 
  what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
  lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
  models = NULL,
  reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
  alg = c("vpdtw"),
  init.coef = c(0, 1, 0),
  n.traces = NULL,
  n.zeros = 0,
  scale = FALSE, 
  plot_it = FALSE,
  penalty = 1, #default 5
  maxshift = 75, #
  verbose = FALSE,
  progress_bar = FALSE, 
  cl = cl)
plot_chroms(rt_standards21_pre_ligustrocide, lambdas = '280')
standards21_ligustrocide_peaks = get_peaks(rt_standards21_pre_ligustrocide, lambdas = c('280'), fit = "egh", amp_thresh = 1800, cl = cl)
standards21_ligustrocide_peak_table = get_peaktable(standards21_ligustrocide_peaks, response = 'height')
print(standards21_ligustrocide_peak_table$pk_meta)



rt_standards21_pre_chlorogenic_acid = preprocess(standards21_pre_chlorogenic_acid, 
                                         dim1=seq(0,45,.05),
                                         dim2=seq(210,400,1),
                                         remove.time.baseline = TRUE,
                                         cl=cl)
plot_chroms(rt_standards21_pre_chlorogenic_acid, lambdas = '280')

rt_standards21_pre_chlorogenic_acid = c(rt_standards21_pre_chlorogenic_acid, copy(rt_standards21_pre_chlorogenic_acid))
standards21_chlorogenic_acid_peaks = get_peaks(rt_standards21_pre_chlorogenic_acid, lambdas = c('280'), fit = "egh", amp_thresh = 309000, cl = cl)
standards21_chlorogenic_acid_peak_table = get_peaktable(standards21_chlorogenic_acid_peaks, response = 'height')
print(standards21_chlorogenic_acid_peak_table$pk_meta)



rt_standards21_pre_vanillic_acid = preprocess(standards21_pre_vanillic_acid, 
                                                 dim1=seq(0,45,.05),
                                                 dim2=seq(210,400,1),
                                                 remove.time.baseline = TRUE,
                                                 cl=cl)
plot_chroms(rt_standards21_pre_vanillic_acid, lambdas = '280')

rt_standards21_pre_vanillic_acid = c(rt_standards21_pre_vanillic_acid, copy(rt_standards21_pre_vanillic_acid))
standards21_vanillic_acid_peaks = get_peaks(rt_standards21_pre_vanillic_acid, lambdas = c('280'), fit = "egh", amp_thresh = 3600, cl = cl)
standards21_vanillic_acid_peak_table = get_peaktable(standards21_vanillic_acid_peaks, response = 'height')
print(standards21_vanillic_acid_peak_table$pk_meta)



rt_standards21_pre_ferulic_acid = preprocess(standards21_pre_ferulic_acid, 
                                                 dim1=seq(0,45,.05),
                                                 dim2=seq(210,400,1),
                                                 remove.time.baseline = TRUE,
                                                 cl=cl)
plot_chroms(rt_standards21_pre_ferulic_acid, lambdas = '280')

rt_standards21_pre_ferulic_acid = c(rt_standards21_pre_ferulic_acid, copy(rt_standards21_pre_ferulic_acid))
standards21_ferulic_acid_peaks = get_peaks(rt_standards21_pre_ferulic_acid, lambdas = c('280'), fit = "egh", amp_thresh = 4700, cl = cl)
standards21_ferulic_acid_peak_table = get_peaktable(standards21_ferulic_acid_peaks, response = 'height')
print(standards21_ferulic_acid_peak_table$pk_meta)



rt_standards21_pre_coumaric_acid = preprocess(standards21_pre_coumaric_acid, 
                                                 dim1=seq(0,45,.05),
                                                 dim2=seq(210,400,1),
                                                 remove.time.baseline = TRUE,
                                                 cl=cl)
plot_chroms(rt_standards21_pre_coumaric_acid, lambdas = '280')

rt_standards21_pre_coumaric_acid = c(rt_standards21_pre_coumaric_acid, copy(rt_standards21_pre_coumaric_acid))
standards21_coumaric_acid_peaks = get_peaks(rt_standards21_pre_coumaric_acid, lambdas = c('280'), fit = "egh", amp_thresh = 4700, cl = cl)
standards21_coumaric_acid_peak_table = get_peaktable(standards21_coumaric_acid_peaks, response = 'height')
print(standards21_coumaric_acid_peak_table$pk_meta)



rt_standards21_pre_gallic_acid = preprocess(standards21_pre_gallic_acid, 
                                                 dim1=seq(0,45,.05),
                                                 dim2=seq(210,400,1),
                                                 remove.time.baseline = TRUE,
                                                 cl=cl)
plot_chroms(rt_standards21_pre_gallic_acid, lambdas = '280')

rt_standards21_pre_gallic_acid = c(rt_standards21_pre_gallic_acid, copy(rt_standards21_pre_gallic_acid))
standards21_gallic_acid_peaks = get_peaks(rt_standards21_pre_gallic_acid, lambdas = c('280'), fit = "egh", amp_thresh = 690000, cl = cl)
standards21_gallic_acid_peak_table = get_peaktable(standards21_gallic_acid_peaks, response = 'height')
print(standards21_gallic_acid_peak_table$pk_meta)



#rt_standards21_pre_standard_mix = preprocess(standards21_pre_standard_mix, 
#                                                 dim1=seq(0,45,.01),
#                                                 dim2=seq(210,400,1),
#                                                 remove.time.baseline = TRUE,
#                                                 cl=TRUE)
#plot_chroms(rt_standards21_pre_standard_mix, lambdas = '280', show_legend = FALSE)


rt_standards21_pre_apagenin = preprocess(standards21_pre_apagenin, 
                                                 dim1=seq(0,45,.05),
                                                 dim2=seq(210,400,1),
                                                 remove.time.baseline = TRUE,
                                                 cl=cl)
plot_chroms(rt_standards21_pre_apagenin, lambdas = '280')
#rt_standards21_pre_apagenin = correct_rt(
#  rt_standards21_pre_apagenin, 
#  what = c("corrected.values"), #What to return: either the 'corrected.values' (useful for visual inspection) or the warping 'models' (for further programmatic use).
#  lambdas=c("280"), #can only select one wavelength to align, will use 280 which was what we read at
#  models = NULL,
#  reference = "best", #"best" is default, I can pick different chromatograms based upon hplc_pre_no_std[[chromatogram number]]# hplc_pre_no_std[[1]][,46]  #Index of the sample that is to be considered the reference sample; not sure how "best" is determined
#  alg = c("vpdtw"),
#  init.coef = c(0, 1, 0),
#  n.traces = NULL,
#  n.zeros = 0,
#  scale = FALSE, 
#  plot_it = FALSE,
#  penalty = 1, #default 5
#  maxshift = 75, #
#  verbose = FALSE,
#  progress_bar = FALSE)
plot_chroms(rt_standards21_pre_apagenin, lambdas = '280')
#rt_standards21_pre_apagenin = c(rt_standards21_pre_apagenin, copy(rt_standards21_pre_apagenin))
standards21_apagenin_peaks = get_peaks(rt_standards21_pre_apagenin, lambdas = c('280'), fit = "egh", amp_thresh = 2100, cl = cl)
standards21_apagenin_peak_table = get_peaktable(standards21_apagenin_peaks, response = 'height')
print(standards21_apagenin_peak_table$pk_meta)


###---Align pre-processed peaks---###

#hplc21_names_df #metadata
#hplc21_list #data

#Standardize site and sample type abbreviations

#Four sites were used in 2021:
#Powder Major's Forest - Green Ash
#Tuttle Swamp - Green Ash
#Town of Lee Forest - White Ash
#Little River Conservation Area - White Ash

length(hplc21_list) #length is now 498
length(hplc21_names_df[,1]) #541; #remove the same rows from the metadata as in the list of chromatograms

hplc21_meta = hplc21_names_df[-c(1:23,49,75,101,127,153,179,204,230,256,282,308,333,359,385,411,437,463,489,515,541),]
length(hplc21_meta[,1]) #498; now equal to the hplc21_list

#create new columns to fix sample_id names
hplc21_meta$length = seq(1:length(hplc21_meta$hplc_2021_sample_ids))
hplc21_meta$sample_id_corrected = hplc21_meta$hplc_2021_sample_ids

#Because we were better with sample naming, we do not need to use the stri_replace_asll_regex function to replace multiple strings
#Instead, I will make some manual updates to sample_ids that have some typos
hplc21_meta[1,5] = "21-TS-Pr-01-5-NEYM2" #according to the field datasheets, this is a 2
hplc21_meta[3,5] = "21-PM-Pr-03-6-YEYM1"
hplc21_meta[4,5] = "21-LF-Pr-02-3-YEYM2"
hplc21_meta[9,5] = "21-TS-P0-02-2-YEYM2"
hplc21_meta[76,5] = "21-PM-Pr-04-7-YEYM1"
hplc21_meta[80,5] = "21-PM-Pr-02-1-NEYM2"
hplc21_meta[91,5] = "21-PM-Pr-1-14-YENMX"
hplc21_meta[98,5] = "21-PM-Pr-1-13-NEYM2"
hplc21_meta[99,5] = "221-PM-Pr-02-5-YENMX"
hplc21_meta[145,5] = "21-LF-Po-03-4-YENMX"
hplc21_meta[263,5] = "21-PM-Po-4-15-NEYM0" #it looks like there are a number of samples that were ran twice. This will show up in the original name column. I think I would probably use the second instance of these as I imagine there was a problem with the first run. I will tabulate the corrected ids to find the duplicates then drop the first instance of the duplicates and drop it in the meta df and chromatogram list
hplc21_meta[269,5] = "21-PM-Po-4-13-YEYM0"
hplc21_meta[271,5] = "21-PM-Pr-02-7-YEYM1"
hplc21_meta[272,5] = "21-PM-Po-2-10-NEYM2"
hplc21_meta[273,5] = "21-PM-Po-04-6-YEYM0"
hplc21_meta[289,5] = "21-PM-Po-01-5-NEYM0"
hplc21_meta[312,5] = "21-PM-Po-04-3-NENMX"
hplc21_meta[334,5] = "21-LR-Pr-04-6-YENMX"
hplc21_meta[351,5] = "21-LR-Pr-02-6-NEYM1"
hplc21_meta[374,5] = "21-LR-Pr-1-15-YEYM0" #1 and 2 have the same trt; I just guessed and put 1, if there is a duplicate it must be 2.
hplc21_meta[381,5] = "21-LR-Po-2-?-YEYM0" #I guess we can figure this one out by looking at which number is missing in the 2-x set of replicates
hplc21_meta[382,5] = "21-LR-Po-04-4-NEYM2"
hplc21_meta[417,5] = "21-LR-Po-2-10-NEYM2" #had an asterisk in the sample_id for some reason; maybe another duplicate
hplc21_meta[430,5] = "21-LR-Po-03-3-NEYM1"
hplc21_meta[458,5] = "21-TS-Po-2-14-YENMX"
hplc21_meta[460,5] = "21-TS-Po-02-5-NENMX" #sample_id was NEYM0; field datasheet says NENMX so that's what I put; make sure to remove 21-TS-Po-02-5-NEYM0
hplc21_meta[465,5] = "21-TS-Po-4-13-YEYM1"
hplc21_meta[467,5] = "21-TS-Po-01-4-YEYM2"
hplc21_meta[484,5] = "21-TS-Po-01-8-YEYM0"
hplc21_meta[487,5] = "21-TS-Po-2-15-NENMX"

#update dataframe with species and treatment information; see: https://stackoverflow.com/questions/25372082/create-column-based-on-presence-of-string-pattern-and-ifelse
hplc21_meta2 = hplc21_meta #new variable

#use grepl to find instances of sites and make changes accordingly # | is "or"
hplc21_meta2$species <- ifelse(grepl("TS|PM",hplc21_meta2$sample_id_corrected),
                               'Green','White') #basically if it's not a green site, it's a white site

#do same thing for pre- and post-samples
hplc21_meta2$sample_type <- ifelse(grepl("Pr",hplc21_meta2$sample_id_corrected),
                                   'Pre','Post') #if it's not a pre-sample, it's a post-sample


#create a treatment column from the last 5 characters in the sample ID
hplc21_meta2$treatment <- substr(hplc21_meta2$sample_id_corrected, 
                                                    nchar(hplc21_meta2$sample_id_corrected)-4, 
                                                    nchar(hplc21_meta2$sample_id_corrected))

#check and make sure correct number of treatments with no misspellings
table(hplc21_meta2$treatment) 

#NENMX NEYM0 NEYM1 NEYM2 YENMX YEYM0 YEYM1 YEYM2 - Looks good; there are some duplicates to be removed though
# 62    64    63    66    62    59    62    60 

#lets find and remove the duplicates
table(hplc21_meta2$sample_id_corrected) 

#duplicate sample_ids
#21-LF-Po-04-1-NENMX; 183 and 184; there is no Pr for this so I imagine one of these is the presample; because this is unclear about the number of the sample (which is listed as below in the spreadsheet), I am going to remove both of these
#   there is a 21-LF-Pr-04-1-NEYM1 which makes me wonder if one of these is labeled incorrectly; There is an associated Po sample for this id

#21-LF-Po-4-11-YEYM0; 8 and 194; dropping the first one
#21-LF-Pr-01-6-NEYM2; 30 and 55; dropping the first one
#21-LF-Pr-02-5-NEYM0; 26 and 60; dropping the first one
#21-LF-Pr-4-15-NENMX; 21 and 41; dropping the first one
#21-LR-Po-02-2-YENMX; 10 and 432; dropping the first one
#21-LR-Pr-02-1-NEYM2; 357 and 371; there is no post for this sample; it is possible one of these is Pr and other is Po, but to be the safest, I will drop both of these
#21-LR-Pr-3-13-NEYM0; 331 and 337; dropping the first one
#21-PM-Po-01-5-NEYM0; 289 and 307; there is no post for this sample; it is possible one of these is Pr and other is Po, but to be the safest, I will drop both of these
#21-PM-Po-03-8-YEYM2; 310 and 315; dropping the first one
#21-PM-Po-04-6-YEYM0; 273 and 275; dropping the first one
#21-PM-Po-2-10-NEYM2; 257 and 272; dropping the first one
#21-PM-Po-4-13-YEYM0; 269 and 296; dropping the first one
#21-PM-Pr-02-2-NEYM0; 83, 85, and 123; dropping the first two
#21-PM-Pr-02-7-YEYM1; 113 and 271; dropping the first one
#21-PM-Pr-1-10-YEYM0; 84 and 88; it is possible one of these is Pr and other is Po, but to be the safest, I will drop both of these
#21-TS-Pr-4-12-YEYM2; 220 and 234, also the Po had two 456, 475; dropping 220 and 475
#21-TS-Pr-4-15-YEYM0; 214 and 239; it is possible one of these is Pr and other is Po, but to be the safest, I will drop both of these

#new dataframe before removing rows
hplc21_meta3 = hplc21_meta2
hplc21_list2 = hplc21_list

#-remove duplicate samples-#

#hplc meta
hplc21_meta3 = hplc21_meta3[-c(8,10,21,26,30,83,84,85,88,113,183,184,214,220,239,257,273,289,307,310,331,357,371,475),]
length(hplc21_meta3$sample_id_corrected) #474

#hplc list
hplc21_list2[c(8,10,21,26,30,83,84,85,88,113,183,184,214,220,239,257,273,289,307,310,331,357,371,475)] = NULL
length(hplc21_list2) #474

#replace treatment codes with actual treatment names
hplc21_meta3$treatment_full = 
  stri_replace_all_regex(hplc21_meta3$treatment,
                         pattern=c('NENMX',
                                   'YENMX',
                                   'NEYM2',
                                   'NEYM1',
                                   'NEYM0',
                                   'YEYM2',
                                   'YEYM1',
                                   'YEYM0'
                                   ),
                         replacement=c('No eggs No methyl Jasmonate',
                                       'Yes eggs No methyl Jasmonate',
                                       'No eggs Yes methyl Jasmonate 2 weeks prior',
                                       'No eggs Yes methyl Jasmonate 1 weeks prior',
                                       'No eggs Yes methyl Jasmonate 0 weeks prior',
                                       'Yes eggs Yes methyl Jasmonate 2 weeks prior',
                                       'Yes eggs Yes methyl Jasmonate 1 weeks prior',
                                       'Yes eggs Yes methyl Jasmonate 0 weeks prior'
                                       ),
                         vectorize=FALSE)


#extract EAB Y/N, MeJA Y/N, and MeJA trt 

#eab
hplc21_meta3$eab = 
  stri_replace_all_regex(hplc21_meta3$treatment,
                         pattern=c('NENMX',
                                   'YENMX',
                                   'NEYM2',
                                   'NEYM1',
                                   'NEYM0',
                                   'YEYM2',
                                   'YEYM1',
                                   'YEYM0'
                         ),
                         replacement=c('No',
                                       'Yes',
                                       'No',
                                       'No',
                                       'No',
                                       'Yes',
                                       'Yes',
                                       'Yes'
                         ),
                         vectorize=FALSE)

#meja
hplc21_meta3$meja = 
  stri_replace_all_regex(hplc21_meta3$treatment,
                         pattern=c('NENMX',
                                   'YENMX',
                                   'NEYM2',
                                   'NEYM1',
                                   'NEYM0',
                                   'YEYM2',
                                   'YEYM1',
                                   'YEYM0'
                         ),
                         replacement=c('No',
                                       'No',
                                       'Yes',
                                       'Yes',
                                       'Yes',
                                       'Yes',
                                       'Yes',
                                       'Yes'
                         ),
                         vectorize=FALSE)

#meja
hplc21_meta3$meja_time = 
  stri_replace_all_regex(hplc21_meta3$treatment,
                         pattern=c('NENMX',
                                   'YENMX',
                                   'NEYM2',
                                   'NEYM1',
                                   'NEYM0',
                                   'YEYM2',
                                   'YEYM1',
                                   'YEYM0'
                         ),
                         replacement=c('',
                                       '',
                                       '2 weeks',
                                       '1 week',
                                       '0 weeks',
                                       '2 weeks',
                                       '1 week',
                                       '0 weeks'
                         ),
                         vectorize=FALSE)

#extract tree number by taking _,-, and numeric characters from our sample names subset
hplc21_meta3$tree_id <- gsub("[^_\\0-9]", "", 
                             substr(hplc21_meta3$sample_id_corrected, 
                                    start=10, #because we had a better naming scheme, we can just pull the characters according to their location within the string
                                    stop=13))

#lets replace "D", "0", and "_" with ""
hplc21_meta3$tree_id_new = 
  stri_replace_all_regex(hplc21_meta3$tree_id,
                         pattern=c('0' 
                         ),
                         replacement=c(''
                         ),
                         vectorize=FALSE)

#insert a hyphen between first and second character
hplc21_meta3$tree_id_new2 <- paste0(substr(hplc21_meta3$tree_id_new, 1, 1),
                                    "-", 
                                    substr(hplc21_meta3$tree_id_new, 2, 
                                           nchar(hplc21_meta3$tree_id_new)))

# create a new column "dbh_cat" based on "tree_id_new2"
#set tree size category according to the last two digits of the string we extracted and cleaned:
#NOTE: 9May23: NAs are sample ids with missing data to be addressed at a later date
#tree number 1-8 = 5-10 cm DBH
#            9-16 = 10.01-15 cm DBH


hplc21_meta3$dbh_cat <- ifelse(substr(gsub("[^0-9]", "", hplc21_meta3$tree_id_new2), 2, 3) %in% c("1", "2", "3","4","5","6","7","8"), "5-10",
                                                  ifelse(substr(gsub("[^0-9]", "", hplc21_meta3$tree_id_new2), 2, 3) %in% c("9", "10", "11","12","13","14","15"), "10.01-15",
                                                         NA))

#assign year
hplc21_meta3$sample_year <- 2021
  
#assign site
hplc21_meta3$site <- ifelse(grepl("TS", hplc21_meta3$sample_id_corrected), "Tuttle Swamp",
                                               ifelse(grepl("LR", hplc21_meta3$sample_id_corrected), "Little River",
                                                      ifelse(grepl("LF", hplc21_meta3$sample_id_corrected), "Lee Forest",
                                                             ifelse(grepl("PM", hplc21_meta3$sample_id_corrected), "Powder Major Forest", 
                                                                    NA))))
# print the updated columns
print(names(hplc21_meta3))

#final dataframe
hplc_2021_df = hplc21_meta3[,c(3,5:12,15:18)]
length(hplc_2021_df$hplc21_names) #474

#---alignment---#

#-subset by site-#

hplc_2021_df$row = seq(1:length(hplc_2021_df$sample_id_corrected))

#create subsets
tut = subset(hplc_2021_df,site =="Tuttle Swamp") #125 
pow = subset(hplc_2021_df,site =="Powder Major Forest") #119
lee =subset(hplc_2021_df,site =="Lee Forest") #112
lit = subset(hplc_2021_df,site =="Little River") #118

#hplc data length = 474
print(length(hplc21_list2))

rt_2021_tut = hplc21_list2
rt_2021_pow = hplc21_list2
rt_2021_lee = hplc21_list2
rt_2021_lit = hplc21_list2

length(tut$row)#125
length(pow$row)#119
length(lee$row)#112
length(lit$row)#118

#remove non-tut rows
rt_2021_tut[c(pow$row,lee$row,lit$row)] = NULL #list of elements to be removed from our list
length(rt_2021_tut) #length is now 125

#check ids to confirm we only have tut samples
rt_2021_tut_sample_ids = sapply(rt_2021_tut, function(x) attr(x,"sample_id")) #confirmed

#remove non-pow rows
rt_2021_pow[c(tut$row,lee$row,lit$row)] = NULL #list of elements to be removed from our list
length(rt_2021_pow) #length is now 119

#check names to confirm we only have pow samples
rt_2021_pow_sample_ids = sapply(rt_2021_pow, function(x) attr(x,"sample_id")) #confirmed

#remove non-lee rows
rt_2021_lee[c(tut$row,pow$row,lit$row)] = NULL #list of elements to be removed from our list
length(rt_2021_lee) #length is now 112

#check names to confirm we only have lee samples
rt_2021_lee_sample_ids = sapply(rt_2021_lee, function(x) attr(x,"sample_id")) #confirmed

#remove non-lit rows
rt_2021_lit[c(tut$row,pow$row,lee$row)] = NULL #list of elements to be removed from our list
length(rt_2021_lit) #length is now 118 (not sure why 2 more than length above)

#check names to confirm we only have lit samples
rt_2021_lit_sample_ids = sapply(rt_2021_lit, function(x) attr(x,"sample_id")) #confirmed

###---Align peaks---###

#---pre-process data---#

#2021 - Tuttle Swamp - Green Ash

#change column name to be consistent in list and dataframe
colnames(tut) <- c("hplc_name","sample_name","tree_species","sample_type","treatment","treatment_full","eab","meja","meja_time","tree_size","year","site","row")

#according to the ?find_peaks help
#It is recommended to do pre-processing using the preprocess function before peak detection. Overly high chromatographic resolution 
#can sometimes lead to peaks being split into multiple segments. In this case, it is recommended to increase the smooth_window or 
#reduce the resolution along the time axis by adjusting the dim2 argument before proceeding with further analyses.

#I originally thought I wanted a very fine time resolution but perhaps I do not, I am going to try and pre-process at 0.05 and then 0.1 to see how this
#alters the number of peaks found. I want to see ~70-90 as the HPLC output originally suggested
rt_2021_tut_end = sapply(rt_2021_tut, function(x) attr(x,"end_time"))


tut_pre = preprocess(rt_2021_tut, #for some reason the program ends at 50.016
                     dim1=seq(0,50,.05),
                     dim2=seq(210,340,1),
                     remove.time.baseline = TRUE,
                     cl = cl)#, #I don't know why I removed this previously
#mc.cores= 6,
#mc.cleanup=TRUE) #https://www.rdocumentation.org/packages/parallel/versions/3.4.0/topics/mclapply

#align data
tut_pre_05_corrected_75 = 
  correct_rt(
    tut_pre, 
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
    progress_bar = TRUE,
    cl = cl)

###---Peak finding---###
#adjusting the smooth window changes the number of peaks output/detected below
#0.75 = 55 peaks;0.05 = 60 peaks; 0.04 = 72 peaks (seems close to HPLC determination); 0.035 = 80 peaks (seems close to HPLC determination)  
#0.025 = 94 peaks; 0.02 = 103 peaks; 0.015 = 121 peaks

#use this subset for the report
tut_pre_corrected_75_pks_egh <- get_peaks(tut_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.025, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000, cl = cl) #works

#plot peak finding
plot(tut_pre_corrected_75_pks_egh, index=1, lambda='280')

###---Creating a peak table---###
tut_pre_corrected_75_egh_peak <- get_peaktable(tut_pre_corrected_75_pks_egh, response = 'height') #

#attach metadata
tut_pk_tab <- attach_metadata(peak_table = tut_pre_corrected_75_egh_peak, 
                              metadata = tut, 
                              column = "hplc_name") #this works! 

#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
tut_pk_tab_df = as.data.frame(do.call(cbind, tut_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
tut_pk_tab_rt_df = as.data.frame(do.call(cbind, tut_pk_tab[2]))

#2021 - Powder Major's Forest - Green

#---pre-process data---#
pow_pre = preprocess(rt_2021_pow, 
                     dim1=seq(0,50,.05),
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE,
                     cl = cl)#, #I don't know why I removed this previously

#align data
pow_pre_05_corrected_75 = 
  correct_rt(
    pow_pre, 
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
    progress_bar = TRUE,
    cl = cl)


###---Peak finding---###

#use this subset for the report
pow_pre_corrected_75_pks_egh <- get_peaks(pow_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.025, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000, cl = cl) #works


#plot peak finding
plot(pow_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035) looks like garbage but is it actually?

###---Creating a peak table---###
pow_pre_corrected_75_egh_peak <- get_peaktable(pow_pre_corrected_75_pks_egh, response = 'height') #

#change column name to be consistent in list and dataframe
colnames(pow)  <- c("hplc_name","sample_name","tree_species","sample_type","treatment","treatment_full","eab","meja","meja_time","tree_size","year","site","row")

#attach metadata - finally works
pow_pk_tab <- attach_metadata(peak_table = pow_pre_corrected_75_egh_peak, 
                              metadata = pow, 
                              column = "hplc_name")


#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
pow_pk_tab_df = as.data.frame(do.call(cbind, pow_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
pow_pk_tab_rt_df = as.data.frame(do.call(cbind, pow_pk_tab[2]))

#2021- Town of Lee Forest - White

#change column name to be consistent in list and dataframe
colnames(lee)  <- c("hplc_name","sample_name","tree_species","sample_type","treatment","treatment_full","eab","meja","meja_time","tree_size","year","site","row")

#rt_2021_lee = compact(rt_2021_lee) #NULL was causing an issue and still retained in the list for some reason

#pre-process
lee_pre = preprocess(rt_2021_lee,
                     dim1=seq(0,50,.05),
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE,
                     cl = cl) 

#align data
lee_pre_05_corrected_75 = 
  correct_rt(
    lee_pre, 
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
    progress_bar = TRUE, 
    cl = cl)

###---Peak finding---###
#adjusting the smooth window changes the number of peaks output/detected below
#using 0.035 here as it's what I used above = 73 peaks (seems close to HPLC determination)  

#use this subset for the report
lee_pre_corrected_75_pks_egh <- get_peaks(lee_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.025, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000, cl = cl) #works

#plot peak finding
plot(lee_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035) looks like garbage but is it actually?

###---Creating a peak table---###
lee_pre_corrected_75_egh_peak <- get_peaktable(lee_pre_corrected_75_pks_egh, response = 'height') #

#attach metadata
lee_pk_tab <- attach_metadata(peak_table = lee_pre_corrected_75_egh_peak, 
                              metadata = lee, 
                              column = "hplc_name") #this works! 

#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
lee_pk_tab_df = as.data.frame(do.call(cbind, lee_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
lee_pk_tab_rt_df = as.data.frame(do.call(cbind, lee_pk_tab[2]))

###---Little River Conservation Area - White Ash---###

#change column name to be consistent in list and dataframe
colnames(lit)  <- c("hplc_name","sample_name","tree_species","sample_type","treatment","treatment_full","eab","meja","meja_time","tree_size","year","site","row")

#pre-process
lit_pre = preprocess(rt_2021_lit,
                     dim1=seq(0,50,.05),
                     dim2=seq(274,284,2),
                     remove.time.baseline = TRUE,
                     cl = cl) 

#align data
lit_pre_05_corrected_75 = 
  correct_rt(
    lit_pre, 
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
    progress_bar = TRUE, 
    cl = cl)

###---Peak finding---###

#use this subset for the report
lit_pre_corrected_75_pks_egh <- get_peaks(lit_pre_05_corrected_75, 
                                          lambdas = c('280'), 
                                          sd.max=25, #adjusting this higher seems neglible.
                                          fit="egh", 
                                          smooth_window=0.025, #lower smooth finds more peaks but tightens up the peak approximation; 0.001
                                          max.iter=2000, cl = cl) #works

#plot peak finding
plot(lit_pre_corrected_75_pks_egh, index=1, lambda='280') #this plot (with 0.035)

###---Creating a peak table---###
lit_pre_corrected_75_egh_peak <- get_peaktable(lit_pre_corrected_75_pks_egh, response = 'height') #

#attach metadata
lit_pk_tab <- attach_metadata(peak_table = lit_pre_corrected_75_egh_peak, 
                              metadata = lit, 
                              column = "hplc_name") #this works! 

#convert our nested list into a df so we can do stats. According to: https://sparkbyexamples.com/r-programming/convert-list-to-r-dataframe/#:~:text=frame()%20converts%20the%20nested,will%20be%20nested%20list%20names.
lit_pk_tab_df = as.data.frame(do.call(cbind, lit_pk_tab[c(1,3)])) #1 is the absorbance information, 3 is metadata
lit_pk_tab_rt_df = as.data.frame(do.call(cbind, lit_pk_tab[2]))
 

#stop cluster for parallel processing
stopCluster(cl)

