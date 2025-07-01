###---Storing standards data in a matrix for easy iterative access---###

###---2020---###
#storing sd of the Gallic Acid peak in a variable
#gallic_acid_sd = as.numeric(standards_gallic_acid_peak_table[6]) 
#summing and averaging the difference between the start rt, peak rt, and end rt, and adding the sd to that to use as a standard deviation in which to search
gallic_acid_sd = (((standards_gallic_acid_peak_table$pk_meta[3,1])-(standards_gallic_acid_peak_table$pk_meta[4,1])
                   +(standards_gallic_acid_peak_table$pk_meta[5,1])-(standards_gallic_acid_peak_table$pk_meta[3,1]))/2)+standards_gallic_acid_peak_table$pk_meta[6,1]
#storing rt of Gallic Acid peak
gallic_acid_peak_rt = standards_gallic_acid_peak_table$pk_meta[3,1]


#storing sd of the Chlorogenic Acid peak in a variable
#chlorogenic_acid_sd = as.numeric(standards_chlorogenic_acid_peak_table[6]) 
#summing and averaging the difference between the start rt, peak rt, and end rt, and adding the sd to that to use as a standard deviation in which to search
chlorogenic_acid_sd = (((standards_chlorogenic_acid_peak_table$pk_meta[3,1])-(standards_chlorogenic_acid_peak_table$pk_meta[4,1])
                        +(standards_chlorogenic_acid_peak_table$pk_meta[5,1])-(standards_chlorogenic_acid_peak_table$pk_meta[3,1]))/2)+standards_chlorogenic_acid_peak_table$pk_meta[6,1]
#storing rt of Chlorogenic Acid peak
chlorogenic_acid_peak_rt = standards_chlorogenic_acid_peak_table$pk_meta[3,1]


#storing sd of the Vanillic Acid peak in a variable
#vanillic_acid_sd = as.numeric(standards_vanillic_acid_peak_table[6]) 
vanillic_acid_sd = (((standards_vanillic_acid_peak_table$pk_meta[3,1])-(standards_vanillic_acid_peak_table$pk_meta[4,1])
                    +(standards_vanillic_acid_peak_table$pk_meta[5,1])-(standards_vanillic_acid_peak_table$pk_meta[3,1]))/2)+standards_vanillic_acid_peak_table$pk_meta[6,1]
#storing rt of Vanillic Acid peak
vanillic_acid_peak_rt = standards_vanillic_acid_peak_table$pk_meta[3,1]


#storing sd of the eleutheroside b peak in a variable
#eleutheroside_b_sd = standards_eleutheroside_b_peak_table$pk_meta[6,1]
eleutheroside_b_sd = (((standards_eleutheroside_b_peak_table$pk_meta[3,1])-(standards_eleutheroside_b_peak_table$pk_meta[4,1])
                       +(standards_eleutheroside_b_peak_table$pk_meta[5,1])-(standards_eleutheroside_b_peak_table$pk_meta[3,1]))/2)+standards_eleutheroside_b_peak_table$pk_meta[6,1]
#storing rt of eleutheroside b peak
eleutheroside_b_peak_rt = standards_eleutheroside_b_peak_table$pk_meta[3,1]


#storing sd of the Coumaric Acid peak in a variable
#coumaric_acid_sd = as.numeric(standards_coumaric_acid_peak_table[6])
coumaric_acid_sd = (((standards_coumaric_acid_peak_table$pk_meta[3,1])-(standards_coumaric_acid_peak_table$pk_meta[4,1])
                    +(standards_coumaric_acid_peak_table$pk_meta[5,1])-(standards_coumaric_acid_peak_table$pk_meta[3,1]))/2)+standards_coumaric_acid_peak_table$pk_meta[6,1]
#storing rt of Coumaric Acid peak
coumaric_acid_peak_rt = standards_coumaric_acid_peak_table$pk_meta[3,1]


#storing sd of the Ferulic Acid peak in a variable
#ferulic_acid_sd = as.numeric(standards_ferulic_acid_peak_table[6])
ferulic_acid_sd = (((standards_ferulic_acid_peak_table$pk_meta[3,1])-(standards_ferulic_acid_peak_table$pk_meta[4,1])
                    +(standards_ferulic_acid_peak_table$pk_meta[5,1])-(standards_ferulic_acid_peak_table$pk_meta[3,1]))/2)+standards_ferulic_acid_peak_table$pk_meta[6,1]
#storing rt of Ferulic Acid peak
ferulic_acid_peak_rt = standards_ferulic_acid_peak_table$pk_meta[3,1]


#storing sd of the Verbascoside peak in a variable
#verbascoside_sd = as.numeric(standards_verbascoside_peak_table[6]) 
verbascoside_sd = (((standards_verbascoside_peak_table$pk_meta[3,1])-(standards_verbascoside_peak_table$pk_meta[4,1])
                    +(standards_verbascoside_peak_table$pk_meta[5,1])-(standards_verbascoside_peak_table$pk_meta[3,1]))/2)+standards_verbascoside_peak_table$pk_meta[6,1]
#storing rt of Verbascoside peak
verbascoside_peak_rt = standards_verbascoside_peak_table$pk_meta[3,1]


#storing sd of the Oleuropein peak in a variable
#oleuropein_sd = as.numeric(standards_oleuropein_peak_table[6]) 
oleuropein_sd = (((standards_oleuropein_peak_table$pk_meta[3,1])-(standards_oleuropein_peak_table$pk_meta[4,1])
                  +(standards_oleuropein_peak_table$pk_meta[5,1])-(standards_oleuropein_peak_table$pk_meta[3,1]))/2)+standards_oleuropein_peak_table$pk_meta[6,1]
#storing rt of Oleuropein peak
oleuropein_peak_rt = standards_oleuropein_peak_table$pk_meta[3,1]


#storing sd of the Pinoresinol peak in a variable
#pinoresinol_sd = as.numeric(standards_pinoresinol_peak_table[6]) 
pinoresinol_sd = (((standards_pinoresinol_peak_table$pk_meta[3,1])-(standards_pinoresinol_peak_table$pk_meta[4,1])
                   +(standards_pinoresinol_peak_table$pk_meta[5,1])-(standards_pinoresinol_peak_table$pk_meta[3,1]))/2)+standards_pinoresinol_peak_table$pk_meta[6,1]
#storing rt of Pinoresinol peak
pinoresinol_peak_rt = standards_pinoresinol_peak_table$pk_meta[3,1]


#storing sd of the Ligustrocide peak in a variable
#ligustrocide_sd = standards_ligustrocide_peak_table$pk_meta[6,1]
ligustrocide_sd = (((standards_ligustrocide_peak_table$pk_meta[3,1])-(standards_ligustrocide_peak_table$pk_meta[4,1])
                    +(standards_ligustrocide_peak_table$pk_meta[5,1])-(standards_ligustrocide_peak_table$pk_meta[3,1]))/2)+standards_ligustrocide_peak_table$pk_meta[6,1]
#storing rt of the Ligustrocide peak in a variable
ligustrocide_peak_rt = standards_ligustrocide_peak_table$pk_meta[3,1]


#storing sd of the Luteolin peak in a variable
#luteolin_sd = as.numeric(standards_luteolin_peak_table[6])
luteolin_sd = (((standards_luteolin_peak_table$pk_meta[3,1])-(standards_luteolin_peak_table$pk_meta[4,1])
                +(standards_luteolin_peak_table$pk_meta[5,1])-(standards_luteolin_peak_table$pk_meta[3,1]))/2)+standards_luteolin_peak_table$pk_meta[6,1]
#storing rt of Luteolin peak in a variable 
luteolin_peak_rt = standards_luteolin_peak_table$pk_meta[3,1]


#storing sd of the Apagenin peak in a variable
#apagenin_sd = standards_apagenin_peak_table$pk_meta[6,1]
apagenin_sd = (((standards_apagenin_peak_table$pk_meta[3,1])-(standards_apagenin_peak_table$pk_meta[4,1])
                +(standards_apagenin_peak_table$pk_meta[5,1])-(standards_apagenin_peak_table$pk_meta[3,1]))/2)+standards_apagenin_peak_table$pk_meta[6,1]
#storing rt of the Apagenin peak in a variable 
apagenin_peak_rt = standards_apagenin_peak_table$pk_meta[3,1]


#Creating the standard data table for iterative access

#storing like data in lists of the identical ordering
standard_name_list = list("Gallic Acid", "Chlorogenic Acid", "Vanillic Acid", "Eleutheroside B", "Coumaric Acid", 
             "Ferulic Acid", "Verbascoside", "Oleuropein", "Pinoresinol", "Ligustrocide", 
             "Luteolin", "Apagenin")
#creating a data frame from the list
standard_name = data.frame(standard_name_list)
#setting the column names to empty strings to prevent errors with rbind
colnames(standard_name) = rep("", 12)
#print to see if it looks right
print(standard_name)


#storing like data in lists of identical ordering
rt_list = list(gallic_acid_peak_rt, chlorogenic_acid_peak_rt, vanillic_acid_peak_rt, 
       eleutheroside_b_peak_rt, coumaric_acid_peak_rt, ferulic_acid_peak_rt, 
       verbascoside_peak_rt, oleuropein_peak_rt, pinoresinol_peak_rt, ligustrocide_peak_rt, 
       luteolin_peak_rt, apagenin_peak_rt)
#creating a data frame from the list
standard_peak_rt = data.frame(rt_list)
#setting the column names to empty strings to prevent errors with rbind
colnames(standard_peak_rt) = rep("", 12)
#print to see if it looks right
print(standard_peak_rt)

#storing like data in lists of identical ordering
sd_list = list(gallic_acid_sd, chlorogenic_acid_sd, vanillic_acid_sd, eleutheroside_b_sd, 
       coumaric_acid_sd, ferulic_acid_sd, verbascoside_sd, oleuropein_sd, pinoresinol_sd, 
       ligustrocide_sd, luteolin_sd, apagenin_sd)
#creating a data frame from the list
standard_peak_sd = data.frame(sd_list)
#using lapply to round all data to 3 sig figs
standard_peak_sd = as.data.frame(lapply(standard_peak_sd, function(x) round(x, 3)))
#setting column names to empty strings to prevent errors with rbind
colnames(standard_peak_sd) = rep("", 12)
#print to see if it looks right
print(standard_peak_sd)


#Binding all individual data frames into a table with rbind
standards_data_table = rbind(standard_name, standard_peak_rt, standard_peak_sd)
#printing to show data structure
print(standards_data_table)








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




###---2021---###
#storing sd of the Gallic Acid peak in a variable
#gallic_acid_21_sd = as.numeric(standards21_gallic_acid_peak_table[6]) 
gallic_acid_21_sd = (((standards21_gallic_acid_peak_table$pk_meta[3,1])-(standards21_gallic_acid_peak_table$pk_meta[4,1])
                      +(standards21_gallic_acid_peak_table$pk_meta[5,1])-(standards21_gallic_acid_peak_table$pk_meta[3,1]))/2)+standards21_gallic_acid_peak_table$pk_meta[6,1] #summing and averaging the difference between the start rt, peak rt, and end rt, and adding the sd to that to use as a standard deviation in which to search
#gallic_acid_21_sd = as.numeric(standards21_gallic_acid_peak_table[6]) * 2
#storing rt of Gallic Acid peak
gallic_acid_21_peak_rt = standards21_gallic_acid_peak_table$pk_meta[3,1]


#storing sd of the Chlorogenic Acid peak in a variable
#chlorogenic_acid_21_sd = as.numeric(standards21_chlorogenic_acid_peak_table[6]) 
chlorogenic_acid_21_sd = (((standards21_chlorogenic_acid_peak_table$pk_meta[3,1])-(standards21_chlorogenic_acid_peak_table$pk_meta[4,1])
                           +(standards21_chlorogenic_acid_peak_table$pk_meta[5,1])-(standards21_chlorogenic_acid_peak_table$pk_meta[3,1]))/2)+standards21_chlorogenic_acid_peak_table$pk_meta[6,1] #summing and averaging the difference between the start rt, peak rt, and end rt to use as a standard deviation in which to search
#chlorogenic_acid_21_sd = as.numeric(standards21_chlorogenic_acid_peak_table[6]) * 2
#storing rt of Chlorogenic Acid peak
chlorogenic_acid_21_peak_rt = standards21_chlorogenic_acid_peak_table$pk_meta[3,1]


#storing sd of the Vanillic Acid peak in a variable
#vanillic_acid_21_sd = as.numeric(standards21_vanillic_acid_peak_table[6]) 
vanillic_acid_21_sd = (((standards21_vanillic_acid_peak_table$pk_meta[3,1])-(standards21_vanillic_acid_peak_table$pk_meta[4,1])
                        +(standards21_vanillic_acid_peak_table$pk_meta[5,1])-(standards21_vanillic_acid_peak_table$pk_meta[3,1]))/2)+standards21_vanillic_acid_peak_table$pk_meta[6,1]
#vanillic_acid_21_sd = as.numeric(standards21_vanillic_acid_peak_table[6]) * 2
#storing rt of Vanillic Acid peak
vanillic_acid_21_peak_rt = standards21_vanillic_acid_peak_table$pk_meta[3,1]


#storing sd of the Syringin peak in a variable
#syringin_21_sd = standards21_syringin_peak_table$pl_meta[6,1]
syringin_21_sd = (((standards21_syringin_peak_table$pk_meta[3,1])-(standards21_syringin_peak_table$pk_meta[4,1])
                   +(standards21_syringin_peak_table$pk_meta[5,1])-(standards21_syringin_peak_table$pk_meta[3,1]))/2)+standards21_syringin_peak_table$pk_meta[6,1]
#syringin_21_sd = standards21_syringin_peak_table$pk_meta[6,1] * 2
#storing rt of Syringin peak
syringin_21_peak_rt = standards21_syringin_peak_table$pk_meta[3,1]


#storing sd of the Coumaric Acid peak in a variable
#coumaric_acid_21_sd = as.numeric(standards21_coumaric_acid_peak_table[6])
coumaric_acid_21_sd = (((standards21_coumaric_acid_peak_table$pk_meta[3,1])-(standards21_coumaric_acid_peak_table$pk_meta[4,1])
                        +(standards21_coumaric_acid_peak_table$pk_meta[5,1])-(standards21_coumaric_acid_peak_table$pk_meta[3,1]))/2)+standards21_coumaric_acid_peak_table$pk_meta[6,1]
#storing rt of Coumaric Acid peak
coumaric_acid_21_peak_rt = standards21_coumaric_acid_peak_table$pk_meta[3,1]


#storing sd of the Ferulic Acid peak in a variable
#ferulic_acid_21_sd = as.numeric(standards21_ferulic_acid_peak_table[6])
ferulic_acid_21_sd = (((standards21_ferulic_acid_peak_table$pk_meta[3,1])-(standards21_ferulic_acid_peak_table$pk_meta[4,1])
                       +(standards21_ferulic_acid_peak_table$pk_meta[5,1])-(standards21_ferulic_acid_peak_table$pk_meta[3,1]))/2)+standards21_ferulic_acid_peak_table$pk_meta[6,1]
#storing rt of Ferulic Acid peak
ferulic_acid_21_peak_rt = standards21_ferulic_acid_peak_table$pk_meta[3,1]


#storing sd of the Verbascoside peak in a variable
#verbascoside_21_sd = as.numeric(standards21_verbascoside_peak_table[6]) 
verbascoside_21_sd = (((standards21_verbascoside_peak_table$pk_meta[3,1])-(standards21_verbascoside_peak_table$pk_meta[4,1])
                       +(standards21_verbascoside_peak_table$pk_meta[5,1])-(standards21_verbascoside_peak_table$pk_meta[3,1]))/2)+standards21_verbascoside_peak_table$pk_meta[6,1]
#storing rt of Verbascoside peak
verbascoside_21_peak_rt = standards21_verbascoside_peak_table$pk_meta[3,1]


#storing sd of the Oleuropein peak in a variable
#oleuropein_21_sd = as.numeric(standards21_oleuropein_peak_table[6]) 
oleuropein_21_sd = (((standards21_oleuropein_peak_table$pk_meta[3,1])-(standards21_oleuropein_peak_table$pk_meta[4,1])
                     +(standards21_oleuropein_peak_table$pk_meta[5,1])-(standards21_oleuropein_peak_table$pk_meta[3,1]))/2)+standards21_oleuropein_peak_table$pk_meta[6,1]
#storing rt of Oleuropein peak
oleuropein_21_peak_rt = standards21_oleuropein_peak_table$pk_meta[3,1]


#storing sd of the Pinoresinol peak in a variable
#pinoresinol_21_sd = as.numeric(standards21_pinoresinol_peak_table[6]) 
pinoresinol_21_sd = (((standards21_pinoresinol_peak_table$pk_meta[3,1])-(standards21_pinoresinol_peak_table$pk_meta[4,1])
                      +(standards21_pinoresinol_peak_table$pk_meta[5,1])-(standards21_pinoresinol_peak_table$pk_meta[3,1]))/2)+standards21_pinoresinol_peak_table$pk_meta[6,1]
#storing rt of Pinoresinol peak
pinoresinol_21_peak_rt = standards21_pinoresinol_peak_table$pk_meta[3,1]


#storing sd of the Ligustrocide peak in a variable
#ligustrocide_21_sd = standards21_ligustrocide_peak_table$pk_meta[6,1]
ligustrocide_21_sd = (((standards21_ligustrocide_peak_table$pk_meta[3,1])-(standards21_ligustrocide_peak_table$pk_meta[4,1])
                       +(standards21_ligustrocide_peak_table$pk_meta[5,1])-(standards21_ligustrocide_peak_table$pk_meta[3,1]))/2)+standards21_ligustrocide_peak_table$pk_meta[6,1]
#storing rt of the Ligustrocide peak in a variable
ligustrocide_21_peak_rt = standards21_ligustrocide_peak_table$pk_meta[3,1]


#storing sd of the Luteolin peak in a variable
#luteolin_21_sd = as.numeric(standards21_luteolin_peak_table[6])
luteolin_21_sd = (((standards21_luteolin_peak_table$pk_meta[3,1])-(standards21_luteolin_peak_table$pk_meta[4,1])
                   +(standards21_luteolin_peak_table$pk_meta[5,1])-(standards21_luteolin_peak_table$pk_meta[3,1]))/2)+standards21_luteolin_peak_table$pk_meta[6,1]
#storing rt of Luteolin peak in a variable 
luteolin_21_peak_rt = standards21_luteolin_peak_table$pk_meta[3,1]


#storing sd of the Apagenin peak in a variable
#apagenin_21_sd = standards21_apagenin_peak_table$pk_meta[6,1]
apagenin_21_sd = (((standards21_apagenin_peak_table$pk_meta[3,1])-(standards21_apagenin_peak_table$pk_meta[4,1])
                   +(standards21_apagenin_peak_table$pk_meta[5,1])-(standards21_apagenin_peak_table$pk_meta[3,1]))/2)+standards21_apagenin_peak_table$pk_meta[6,1]
#storing rt of the Apagenin peak in a variable 
apagenin_21_peak_rt = standards21_apagenin_peak_table$pk_meta[3,1]


#Creating the standard data table for iterative access

#storing like data in lists of the identical ordering
standard_name_list_21 = list("Gallic Acid", "Chlorogenic Acid", "Vanillic Acid", "Syringin", "Coumaric Acid", 
                          "Ferulic Acid", "Verbascoside", "Oleuropein", "Pinoresinol", "Ligustrocide", 
                          "Luteolin", "Apagenin")
#creating a data frame from the list
standard_name_21 = data.frame(standard_name_list_21)
#setting the column names to empty strings to prevent errors with rbind
colnames(standard_name_21) = rep("", 12)
#print to see if it looks right
print(standard_name_21)


#storing like data in lists of identical ordering
rt_list_21 = list(gallic_acid_21_peak_rt, chlorogenic_acid_21_peak_rt, vanillic_acid_21_peak_rt, 
               syringin_21_peak_rt, coumaric_acid_21_peak_rt, ferulic_acid_21_peak_rt, 
               verbascoside_21_peak_rt, oleuropein_21_peak_rt, pinoresinol_21_peak_rt, ligustrocide_21_peak_rt, 
               luteolin_21_peak_rt, apagenin_21_peak_rt)
#creating a data frame from the list
standard_peak_rt_21 = data.frame(rt_list_21)
#setting the column names to empty strings to prevent errors with rbind
colnames(standard_peak_rt_21) = rep("", 12)
#print to see if it looks right
print(standard_peak_rt_21)

#storing like data in lists of identical ordering
sd_list_21 = list(gallic_acid_21_sd, chlorogenic_acid_21_sd, vanillic_acid_21_sd, syringin_21_sd, 
               coumaric_acid_21_sd, ferulic_acid_21_sd, verbascoside_21_sd, oleuropein_21_sd, pinoresinol_21_sd, 
               ligustrocide_21_sd, luteolin_21_sd, apagenin_21_sd)
#creating a data frame from the list
standard_peak_sd_21 = data.frame(sd_list_21)
#using lapply to round all data to 3 sig figs
standard_peak_sd_21 = as.data.frame(lapply(standard_peak_sd_21, function(x) round(x, 3)))
#setting column names to empty strings to prevent errors with rbind
colnames(standard_peak_sd_21) = rep("", 12)
#print to see if it looks right
print(standard_peak_sd_21)


#Binding all individual data frames into a table with rbind
standards_data_table_21 = rbind(standard_name_21, standard_peak_rt_21, standard_peak_sd_21)
#printing to show data structure
print(standards_data_table_21)


