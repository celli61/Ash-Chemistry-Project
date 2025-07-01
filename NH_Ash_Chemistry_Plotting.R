########################################################################################
# Using compiled and structured data for graphical analysis of data in different forms #
########################################################################################

okabe_ito_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)


###---Plotting Standard Presence by Treatment in Subset---###
# plot_standard_presence(peak_tab, chrom_list, treatment_index, data_table, subset_meta)
#
# This method takes in a specified treatment from a chromatographR chrom_list and plots a histogram of it with vertical
# lines at any peak in the histogram identified as a standard, with an accompanying diamond surrounding part of the 
# histogram identified as a standard, and a horizontal error bar representing the detection range around the identified
# peak.
# 
# @param peak_tab - A chromatographR peak_table object containing identified peaks, absorbance data, and peak metadata.
# @param Chrom_list - A chromatographR chrom_list object containing all treatments of a subset.
# @param treatment_index - The numeric index of the treatment to be analyzed and plotted.
# @param data_table - A data frame containing information about the standards, including RTs and SDs.
# @param subset_meta - A data frame containing metadata about the subset.
#
# @return - A combined plot displaying the standard presence in the specified treatment in the given subset. 
# 
plot_standard_presence <- function(peak_tab, run_index, data_table){
  #Extract and store subset meta data and chrom list
  chrom_list = get(peak_tab$args$chrom_list)
  subset_meta = peak_tab$sample_meta
  
  #Extract absorbance data for the specified treatment and convert to numeric 
  absorbance = get(peak_tab$args$chrom_list)[[run_index]][,'280']
  
  #Extract retention time data for the specified treatment and convert to numeric 
  time = as.numeric(attr(get(peak_tab$args$chrom_list)[[run_index]], "dimnames")[[1]])
  
  #Combine time and absorbance data into a data frame for plotting
  df = data.frame(time, absorbance) 
  
  #Initialize a list to store plots for each standard
  plot_list = vector(mode = "list", length = 12)
  
  
  #Iterate over the 12 external standards
  for(i in 1:length(data_table[1,])){
    #Plot the absorbance over time data for the current standard, label, and theme properly 
    plot = ggplot(df, aes(x = time, y = absorbance)) + geom_line(size = 0.185) +
      labs(x = "Retention Time", y = "Absorbance", title = data_table[1,i]) + 
      theme(plot.title = element_text(size = 7, color = "#E69F00"), 
            axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5), 
            axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))
    
    #Call the linear search method to find a peak corresponding to the current standard
    standard_index = linear_search(peak_tab$pk_meta[3,], as.numeric(data_table[2,i]), as.numeric(data_table[3,i])) 
    
    #Check if a matching peak was found and the intensity at the peak RT is > 0
    if(standard_index != -1 && peak_tab$tab[run_index, standard_index] != 0){ 
      #Extract and store the retention time of the peak
      pk_time = peak_tab$pk_meta[3,standard_index] 
      pk_time_char = as.character(round(pk_time / 0.05) * 0.05)
      
      #Use rounding and numeric to character conversions to extract absorbance at the peak RT
      pk_absorbance = absorbance[pk_time_char] 
      
      #Calculate the standard deviation range for the peak
      pk_time_min = (pk_time - as.double(data_table[3,i])) 
      pk_time_max = (pk_time + as.double(data_table[3,i]))
      
      #Create a data frame for plotting the peak data 
      pk_df = data.frame(pk_time, pk_absorbance, pk_time_min, pk_time_max) 
      
      #Add a point at the RT and intensity of the identified peak, with a value label,
      #and horizontal error bars to represent SD
      plot = plot + geom_point(data = pk_df, aes(x = pk_time, y = pk_absorbance),
                   size = 1.5, shape = 23, stroke = 0.375, fill = NA, color = "#0072b2") + 
        geom_vline(xintercept = pk_time, linetype = "dashed", color = "#D55E00", size = 0.25) +
        annotate("text", x = pk_time, y = max(unlist(absorbance)), 
                 label = pk_time_char, size = 1.5, color = "#D55E00", hjust = 1.4) +
        geom_errorbarh(data = pk_df, aes(x = pk_time, xmin = pk_time_min, xmax = pk_time_max, 
                                         y = pk_absorbance), size = 0.33, height = 0.04, alpha = 0.75) 
    }
    #Add the plot to the list of plots
    plot_list[[i]] = plot 
  }
  #Combine all the plots in plot_list and add a title to the combined plot
  combined_plot = wrap_plots(plot_list) + 
    plot_annotation(title = paste0(subset_meta[run_index, 8], " ", subset_meta[run_index, 7]," - ", 
                                   names(chrom_list)[run_index], " - ", subset_meta[run_index, 3], " Ash"))
  
  return(combined_plot)
}

#Plotting and saving standard presence plots for all subsets
plot_standard_presence(doe_pk_tab, 43, standards_data_table)
ggsave("doe_plot_standard_presence.pdf", plot = plot_standard_presence(doe_pk_tab, 43, standards_data_table))
plot_standard_presence(fre_pk_tab, 12, standards_data_table)
ggsave("fre_plot_standard_presence.pdf", plot = plot_standard_presence(fre_pk_tab, 12, standards_data_table))
plot_standard_presence(eas_pk_tab, 27, standards_data_table)
ggsave("eas_plot_standard_presence.pdf", plot = plot_standard_presence(eas_pk_tab, 27, standards_data_table))
plot_standard_presence(jen_pk_tab, 24, standards_data_table)
ggsave("jen_plot_standard_presence.pdf", plot = plot_standard_presence(jen_pk_tab, 24, standards_data_table))

plot_standard_presence(tut_pk_tab, 12, standards21_data_table)
ggsave("tut_plot_standard_presence.pdf", plot = plot_standard_presence(tut_pk_tab, 12, standards21_data_table))
plot_standard_presence(pow_pk_tab, 67, standards21_data_table)
ggsave("pow_plot_standard_presence.pdf", plot = plot_standard_presence(pow_pk_tab, 67, standards21_data_table))
plot_standard_presence(lee_pk_tab, 4, standards21_data_table)
ggsave("lee_plot_standard_presence.pdf", plot = plot_standard_presence(lee_pk_tab, 4, standards21_data_table))
plot_standard_presence(lit_pk_tab, 4, standards21_data_table)
ggsave("lit_plot_standard_presence.pdf", plot = plot_standard_presence(lit_pk_tab, 4, standards21_data_table))


###---Plotting Standard Presence by Tree Size---###
# plot_standard_presence_tree_size(std_pres_tab, subset_meta)
#
# This method uses ggplot2 to plot the detection of standards in the treatments of a subset by treatment tree size.
#
# @param std_pres_tab - A Boolean matrix indicating the presence of each of the 12 external standards in the treatments
#                       of a subset.
# @param susbset_meta - A data frame containing metadata about the subset, including tree size, site name, year, and 
#                      tree species.
#
# @return - A combined plot displaying the algorithmic detection of the 12 standards by tree sizes 
#
plot_standard_presence_tree_size = function(std_pres_tab, subset_meta) {
  #storing all unique tree sizes of the current subset in a sorted list
  tree_sizes = unlist(as.list(na.omit(mixedsort(unique(subset_meta[, 6]))))) 
  
  #Creating an empty list to store iteratively created ggplot2 plot objects
  plot_list = list() 
  
  #Iterate over the 12 standards (columns)
  for(i in 1:12) { 
    #creating data frame to store data for plotting
    standard_presence_rate_df = data.frame(Tree_Size = tree_sizes, 
                                           Algorithmic_Detection = vector("numeric", length = length(tree_sizes)), 
                                           Size_Total = vector("numeric", length = length(tree_sizes)))
    
    #Iterate over treatments (rows)
    for(j in 1:length(std_pres_tab[,1])) { 
      #Find which size our current treatment is 
      index = which(unlist(tree_sizes) == subset_meta[j, 6]) 
      
      #Increment the total number of trees with the current tree size
      standard_presence_rate_df[index, 3] = standard_presence_rate_df[index, 3] + 1 
      
      #Check if the current standard is present in the current treatment 
      if(std_pres_tab[j, i] == TRUE) { 
        #Increment the algorithmic detection count for the current tree size 
        standard_presence_rate_df[index, 2] = standard_presence_rate_df[index, 2] + 1 
      }
    }
    #Create and store a bar plot with the number of treatments on the Y axis and the tree sizes on the X axis
    plot_list[[i]] = ggplot(standard_presence_rate_df, aes(factor(Tree_Size, levels = Tree_Size), 
                                                           Algorithmic_Detection)) + 
      geom_col(fill = "#0072B2", color = "#000000") + 
      geom_text(aes(label = paste0("N = ", Size_Total, "\n", Algorithmic_Detection)), vjust = -0.25, size = 1.5) + 
      scale_y_continuous(limits = c(0, max(standard_presence_rate_df$Size_Total) + 10)) + 
      labs(title = colnames(std_pres_tab)[i], x = "Tree Size", y = "Algorithmic Detection") +
      theme(plot.title = element_text(size = 7, color = "#E69F00"), axis.title.x = element_text(size = 5), 
            axis.title.y = element_text(size = 5), axis.text.x = element_text(size = 5), 
            axis.text.y = element_text(size = 5))
  }
  #Combine all plots and add annotations showing the subset site, year, and tree species 
  combined_plot = wrap_plots(plot_list) + 
    plot_annotation(title = paste0(subset_meta[j, 8], " ", subset_meta[j, 7], " - ", 
                                   subset_meta[j, 3], " Ash | N = # of trees of respective size"))
  
  return(combined_plot)
}

#plotting and saving all standard presence detection plots for all subsets
plot_standard_presence_tree_size(doe_standard_presence_table, doe_meta)
ggsave("doe_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(doe_standard_presence_table, doe_meta))
plot_standard_presence_tree_size(fre_standard_presence_table, fre_meta)
ggsave("fre_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(fre_standard_presence_table, fre_meta))
plot_standard_presence_tree_size(eas_standard_presence_table, eas_meta)
ggsave("eas_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(eas_standard_presence_table, eas_meta))
plot_standard_presence_tree_size(jen_standard_presence_table, jen_meta)
ggsave("jen_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(jen_standard_presence_table, jen_meta))

plot_standard_presence_detection(tut_standard_presence_table, tut_meta)
ggsave("tut_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(tut_standard_presence_table, tut_meta))
plot_standard_presence_detection(pow_standard_presence_table, pow_meta)
ggsave("pow_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(pow_standard_presence_table, pow_meta))
plot_standard_presence_detection(lee_standard_presence_table, lee_meta)
ggsave("lee_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(lee_standard_presence_table, lee_meta))
plot_standard_presence_detection(lit_standard_presence_table, lit_meta)
ggsave("lit_plot_standard_presence_detection.png", plot = plot_standard_presence_detection(lit_standard_presence_table, lit_meta))


###---Plotting Algorithmic Detection Rate by Tree Size---###
# plot_standard_presence_rate(std_pres_tab, subset_meta)
#
# @param std_pres_tab - A Boolean matrix indicating the presence of the 12 standards in the treatments of the subset.
# @param subset_meta - A data frame containing the metadata for the susbet, including tree size, site, year, and species.
#
# @return - A combined plot displaying the detection rate of the standards by the different tree sizes.
#
plot_standard_presence_rate = function(std_pres_tab, subset_meta) {
  #Storing all unique tree sizes of the subset as a sorted list 
  tree_sizes = unlist(as.list(na.omit(mixedsort(unique(subset_meta[, 6]))))) 
  
  #Creating an empty list to store ggplot2 plot objects created iteratively 
  plot_list = list() 
  
  #Iterate over the 12 standards (columns)
  for(i in 1:12) {
    #Creating data frame to store data for plotting
    standard_presence_rate_df = data.frame(Tree_Size = tree_sizes, 
                                           Algorithmic_Detection = vector("numeric", length = length(tree_sizes)), 
                                           Size_Total = vector("numeric", length = length(tree_sizes)))
    
    #Iterate over the treatments (rows)
    for(j in 1:length(std_pres_tab[,1])) {
      #Find the size of the current treatment
      index = which(unlist(tree_sizes) == subset_meta[j, 6]) 
      
      #Increment the total number of trees with the current tree size 
      standard_presence_rate_df[index, 3] = standard_presence_rate_df[index, 3] + 1 
      
      #Check if the crrent standard is present in the current treatment 
      if(std_pres_tab[j, i] == TRUE) { 
        #Increment the algorithmic detection count for the treatment's tree size 
        standard_presence_rate_df[index, 2] = standard_presence_rate_df[index, 2] + 1 
      }
    }
    
    #Calculating the detection rate for each tree size 
    standard_presence_rate_df$Algorithmic_Detection_Rate = standard_presence_rate_df$Algorithmic_Detection / standard_presence_rate_df$Size_Total
    
    #Creating and storing a bar plot with the algorithmic detection rate on the Y axis and tree sizes on the X axis
    plot_list[[i]] = ggplot(standard_presence_rate_df, 
                            aes(factor(Tree_Size, levels = Tree_Size), Algorithmic_Detection_Rate)) + 
      geom_col(fill = "#0072B2", color = "#000000") +
      geom_text(aes(label = paste0("N = ", Size_Total, "\n", 
                                   round(Algorithmic_Detection_Rate, 2))), vjust = -0.2, size = 1.2) + 
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + 
      labs(title = colnames(std_pres_tab)[i], x = "Tree Size", y = "Algorithmic Detection Rate") + 
      theme(plot.title = element_text(size = 7, color = "#E69F00"), axis.title.x = element_text(size = 5), 
            axis.title.y = element_text(size = 5), axis.text.x = element_text(size = 5), 
            axis.text.y = element_text(size = 5)) + 
      coord_cartesian(clip = "off")
  }
  
  #Combine all plots and add annotations showing the subset site, year, and tree species
  combined_plot = wrap_plots(plot_list) + 
    plot_annotation(title = paste0(subset_meta[j, 8], " ", subset_meta[j, 7], " - ", 
                                   subset_meta[j, 3], " Ash | N = # of trees of respective size"))
  
  return(combined_plot)
}

#plotting and saving all standard presence rate table plots for all susbets
plot_standard_presence_rate(doe_standard_presence_table, doe_meta)
ggsave("doe_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(doe_standard_presence_table, doe_meta))
plot_standard_presence_rate(fre_standard_presence_table, fre_meta)
ggsave("fre_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(fre_standard_presence_table, fre_meta))
plot_standard_presence_rate(eas_standard_presence_table, eas_meta)
ggsave("eas_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(eas_standard_presence_table, eas_meta))
plot_standard_presence_rate(jen_standard_presence_table, jen_meta)
ggsave("jen_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(jen_standard_presence_table, jen_meta))

plot_standard_presence_rate(tut_standard_presence_table, tut_meta)
ggsave("tut_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(tut_standard_presence_table, tut_meta))
plot_standard_presence_rate(pow_standard_presence_table, pow_meta)
ggsave("pow_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(pow_standard_presence_table, pow_meta))
plot_standard_presence_rate(lee_standard_presence_table, lee_meta)
ggsave("lee_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(lee_standard_presence_table, lee_meta))
plot_standard_presence_rate(lit_standard_presence_table, lit_meta)
ggsave("lit_plot_standard_presence_rate.pdf", plot = plot_standard_presence_rate(lit_standard_presence_table, lit_meta))


###---Plotting Algorithmic Detetction by Treatment Type---###
# plot_standard_presence_rate_treatment(std_pres_tab, subset_meta)
# 
# This method uses ggplot2 to plot the detection rate of the 12 standards by treatment type.
#
# @param std_pres_tab - A Boolean matrix indicating the presence of the 12 standards in the treatments of the subset.
# @param subset_meta - A data frame containing the metadata for the susbet, including tree size, site, year, and species.
#
# @return - A combined plot displaying the detection rate of standards over the treatment types of the subset 
#
plot_standard_presence_rate_treatment = function(std_pres_tab, subset_meta) {
  #storing all unique treatments in the subset as a sorted list, excluding "Unkno" entries
  treatments = unlist(as.list(setdiff((mixedsort(unique(str_to_title(subset_meta[, 5])))), "Unkno")))
  
  #Creating an empty list to store ggplot2 plot objects created iteratively 
  plot_list = list() 
  
  #Iterate over the 12 standards (columns)
  for(i in 1:12) {
    #Creating data frame to store data for plotting
    standard_presence_rate_df = data.frame(Treatment = treatments, 
                                           Algorithmic_Detection = vector("numeric", length = length(treatments)),
                                           Treatment_Total = vector("numeric", length = length(treatments)))
    
    #Iterate over the treatments (rows)
    for(j in 1:length(std_pres_tab[,1])) { 
      #Find the treatment type of the current treatment
      index = which(unlist(treatments) == str_to_title(subset_meta[j, 5])) 
      
      #Increment the total number of trees of the cuurent treatment's treatment type
      standard_presence_rate_df[index, 3] = standard_presence_rate_df[index, 3] + 1
      
      #Check if the current standard is present in the current treatment 
      if(std_pres_tab[j, i] == TRUE) { 
        #Increment the algorithmic detection count for the current treatment's treatment type
        standard_presence_rate_df[index, 2] = standard_presence_rate_df[index, 2] + 1 
      }
    }
    #Calculate the detection rate for each treatment type 
    standard_presence_rate_df$Algorithmic_Detection_Rate = standard_presence_rate_df$Algorithmic_Detection / standard_presence_rate_df$Treatment_Total
    
    #Replace NA values with 0 in the detection rate column of the data frame 
    standard_presence_rate_df[is.na(standard_presence_rate_df)] = 0
    
    #Create and store a bar plot with detection rate on the Y axis and treatment type on the X axis 
    plot_list[[i]] = ggplot(standard_presence_rate_df, 
                            aes(factor(Treatment, levels = Treatment), Algorithmic_Detection_Rate)) + 
      geom_col(fill = "#0072B2", color = "#000000") +
      geom_text(aes(label = paste0("N = ", Treatment_Total, "\n", round(Algorithmic_Detection_Rate, 2))), 
                vjust = -0.2, size = 1.2) + 
      scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) + 
      labs(title = colnames(std_pres_tab)[i], x = "Tree Treatment", y = "Algorithmic Detection Rate") + 
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      theme(plot.title = element_text(size = 7, color = "#E69F00"),
            axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5), 
            axis.text.x = element_text(size = 3), axis.text.y = element_text(size = 5)) + 
      coord_cartesian(clip = "off")
  }
  #Combine all plots and add annotations with subset site, year, and tree species
  combined_plot = wrap_plots(plot_list) + 
    plot_annotation(title = paste0(subset_meta[j, 8], " ", subset_meta[j, 7],
                                   " - ", subset_meta[j, 3], " Ash | N = # of trees of respective treatment"))
  
  return(combined_plot)
}

#plotting and saving att standard presence rate by treatment plots for all susbets
plot_standard_presence_rate_treatment(doe_standard_presence_table, doe_meta)
ggsave("doe_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(doe_standard_presence_table, doe_meta))
plot_standard_presence_rate_treatment(fre_standard_presence_table, fre_meta)
ggsave("fre_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(fre_standard_presence_table, fre_meta))
plot_standard_presence_rate_treatment(eas_standard_presence_table, eas_meta)
ggsave("eas_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(eas_standard_presence_table, eas_meta))
plot_standard_presence_rate_treatment(jen_standard_presence_table, jen_meta)
ggsave("jen_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(jen_standard_presence_table, jen_meta))

plot_standard_presence_rate_treatment(tut_standard_presence_table, tut_meta)
ggsave("tut_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(tut_standard_presence_table, tut_meta))
plot_standard_presence_rate_treatment(pow_standard_presence_table, pow_meta)
ggsave("pow_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(pow_standard_presence_table, pow_meta))
plot_standard_presence_rate_treatment(lee_standard_presence_table, lee_meta)
ggsave("lee_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(lee_standard_presence_table, lee_meta))
plot_standard_presence_rate_treatment(lit_standard_presence_table, lit_meta)
ggsave("lit_plot_standard_presence_rate_treatment.pdf", plot = plot_standard_presence_rate_treatment(lit_standard_presence_table, lit_meta))




#Fixing 2021 metadata column names
colnames(tut_meta) = c("hplc_name", "sample_name", "tree_species", "sample_type", "treatment", "tree_size", "year", "site")
colnames(pow_meta) = c("hplc_name", "sample_name", "tree_species", "sample_type", "treatment", "tree_size", "year", "site")
colnames(lee_meta) = c("hplc_name", "sample_name", "tree_species", "sample_type", "treatment", "tree_size", "year", "site")
colnames(lit_meta) = c("hplc_name", "sample_name", "tree_species", "sample_type", "treatment", "tree_size", "year", "site")



plot_subset_peaks_quantity = function(pk_quant_tab, meta, group_var = c("sample_type", "treatment", "tree_size")) {
  #ensure x variable parameter is from meta data
  group_var = match.arg(group_var)
  #print(group_var)
  
  #store unique values for group var
  group_types = setdiff(mixedsort(unique(na.omit(meta[, group_var]))), "Unkno")
  #group_types = mixedsort(unique(meta[, group_var]))
  #print(group_types)
  
  #make plot list
  plot_list = list()
  
  #iterate over peaks of peak quantity table
  for(pk_index in seq_along(pk_quant_tab)) {
    #create plotting data df
    quant_by_group_var_df = data.frame(group_var = group_types, quant_sum = vector("numeric", length = length(group_types)), 
                                       group_total = vector("numeric",length(group_types)))
    #print(quant_by_group_var_df)
    
    #iterate over runs
    for(run_index in seq_along(pk_quant_tab[,1])) {
      #store current runs group var type
      run_group_var = meta[meta$hplc_name == row.names(pk_quant_tab)[run_index], group_var]
      #print(run_group_var)
      
      #increment the group total for the current runs group var
      row_index = which(quant_by_group_var_df$group_var == run_group_var)
      quant_by_group_var_df[row_index, "group_total"] = quant_by_group_var_df[row_index, "group_total"] + 1
      
      #add quantity value to quantity sum
      quant_by_group_var_df[row_index, "quant_sum"] = quant_by_group_var_df[row_index, "quant_sum"] + pk_quant_tab[run_index, pk_index]
    }
    #print(quant_by_group_var_df)
    
    #store group variables, qauntity sums, and group totals as vectors for plotting
    group_var_vect = as.character(quant_by_group_var_df$group_var)
    #print(group_var_vect)
    quant_sum_vect = unlist(quant_by_group_var_df$quant_sum)
    group_total_vect = unlist(quant_by_group_var_df$group_total)
    #print(group_total_vect)
    
    plot_data = data.frame(group_var = factor(group_var_vect, levels = group_types), average_quant = quant_sum_vect / group_total_vect, group_total = group_total_vect)
    #print(plot_data)
    
    quant_barplot = ggplot(plot_data, aes(x = group_var, y = average_quant, fill = group_var)) + 
      geom_bar(stat = "identity") +
      scale_fill_manual(values = okabe_ito_colors[1:length(levels(plot_data$group_var))]) +
      geom_text(aes(label = paste0(group_total, "/", length(pk_quant_tab[, 1]))), vjust = -0.5, size = 1.5) +
      scale_x_discrete(labels = abbreviate) + 
      theme_minimal(base_size = 6) + 
      theme(legend.position = "none", 
            axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
            axis.text.y = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            plot.title = element_text(size = 7.5)) +
      labs(
        title = paste("RT:", names(pk_quant_tab)[pk_index]),
        x = tools::toTitleCase(gsub("_", " ", as.character(group_var))),
        y = "Average Proportion \nof Sample"
      ) +
      expand_limits(y = max(plot_data$average_quant) * 1.4)
    
    plot_list[[pk_index]] = quant_barplot
  }
  combined_plot = wrap_plots(plot_list, ncols = 8) + 
    plot_annotation(title = paste0(meta$site[1], " Average Proportion of Sample for All Peaks"), 
                    theme = theme(plot.title = element_text(size = 9)))
  return(combined_plot)
}



#Doe
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/doe_sample_type.png",
       plot = plot_subset_peaks_quantity(doe_quant_tab, doe_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/doe_treatment.png",
       plot = plot_subset_peaks_quantity(doe_quant_tab, doe_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/doe_tree_size.png",
       plot = plot_subset_peaks_quantity(doe_quant_tab, doe_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Fre
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/fre_sample_type.png",
       plot = plot_subset_peaks_quantity(fre_quant_tab, fre_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/fre_treatment.png",
       plot = plot_subset_peaks_quantity(fre_quant_tab, fre_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/fre_tree_size.png",
       plot = plot_subset_peaks_quantity(fre_quant_tab, fre_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Eas
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/eas_sample_type.png",
       plot = plot_subset_peaks_quantity(eas_quant_tab, eas_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/eas_treatment.png",
       plot = plot_subset_peaks_quantity(eas_quant_tab, eas_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/eas_tree_size.png",
       plot = plot_subset_peaks_quantity(eas_quant_tab, eas_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Jen
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/jen_sample_type.png",
       plot = plot_subset_peaks_quantity(jen_quant_tab, jen_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/jen_treatment.png",
       plot = plot_subset_peaks_quantity(jen_quant_tab, jen_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2020/jen_tree_size.png",
       plot = plot_subset_peaks_quantity(jen_quant_tab, jen_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Tut
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/tut_sample_type.png",
       plot = plot_subset_peaks_quantity(tut_quant_tab, tut_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/tut_treatment.png",
       plot = plot_subset_peaks_quantity(tut_quant_tab, tut_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/tut_tree_size.png",
       plot = plot_subset_peaks_quantity(tut_quant_tab, tut_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Pow
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/pow_sample_type.png",
       plot = plot_subset_peaks_quantity(pow_quant_tab, pow_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/pow_treatment.png",
       plot = plot_subset_peaks_quantity(pow_quant_tab, pow_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/pow_tree_size.png",
       plot = plot_subset_peaks_quantity(pow_quant_tab, pow_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Lee
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/lee_sample_type.png",
       plot = plot_subset_peaks_quantity(lee_quant_tab, lee_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/lee_treatment.png",
       plot = plot_subset_peaks_quantity(lee_quant_tab, lee_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/lee_tree_size.png",
       plot = plot_subset_peaks_quantity(lee_quant_tab, lee_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Lit
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/lit_sample_type.png",
       plot = plot_subset_peaks_quantity(lit_quant_tab, lit_meta, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/lit_treatment.png",
       plot = plot_subset_peaks_quantity(lit_quant_tab, lit_meta, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_peaks/2021/lit_tree_size.png",
       plot = plot_subset_peaks_quantity(lit_quant_tab, lit_meta, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")






plot_subset_peaks_quantity_external_standards = function(pk_quant_tab, meta, standards_data_tab, group_var = c("sample_type", "treatment", "tree_size")) {
  #ensure x variable is a group variable from meta data
  group_var = match.arg(group_var)
  
  #store external standard peak RT and SD values as numeric values
  standard_rts = as.numeric(standards_data_tab[2, ])
  standard_sds = as.numeric(standards_data_tab[3, ])
  
  #store subset peak RTs
  peak_rts = as.numeric(names(pk_quant_tab))
  
  #store unique values for x variable
  group_types = setdiff(mixedsort(unique(na.omit(meta[, group_var]))), "Unkno")
  
  #make a list to store multiple plots
  plot_list = list()
  
  #iterate over 12 external standards
  for (standard_index in seq_along(standard_rts)) {
    #store data for current standard
    standard_name = standards_data_tab[1, standard_index]
    standard_rt = standard_rts[standard_index]
    standard_sd = standard_sds[standard_index]
    
    #get all peak RTs that fall within SD of standard RT
    matching_rts = get_rt_in_sd_list(peak_rts, standard_rt, standard_sd)
    
    #if any peaks were identified within SD
    if (length(matching_rts) != 0) {
      #select the peak with the highest sum of proportions of sample for current standard
      peak_rt = matching_rts[which.max(colSums(pk_quant_tab[, as.character(matching_rts), drop = FALSE]))]
      pk_index = which(names(pk_quant_tab) == as.character(peak_rt))
      
      #put data in DF to add to iteratively
      quant_by_group_var_df = data.frame(group_var = group_types,
                                         quant_sum = numeric(length(group_types)),
                                         group_total = numeric(length(group_types)))
      
      #iterate over samples in subset
      for (run_index in seq_along(pk_quant_tab[,1])) {
        run_group_var = meta[meta$hplc_name == row.names(pk_quant_tab)[run_index], group_var]
        row_index = which(quant_by_group_var_df$group_var == run_group_var)
        quant_by_group_var_df[row_index, "group_total"] = quant_by_group_var_df[row_index, "group_total"] + 1
        quant_by_group_var_df[row_index, "quant_sum"] = quant_by_group_var_df[row_index, "quant_sum"] + pk_quant_tab[run_index, pk_index]
      }
      
      #format data for plotting
      plot_data = data.frame(
        group_var = factor(as.character(quant_by_group_var_df$group_var), levels = group_types),
        average_quant = quant_by_group_var_df$quant_sum / quant_by_group_var_df$group_total,
        group_total = quant_by_group_var_df$group_total
      )
      
      #plot data in barplot
      quant_barplot = ggplot(plot_data, aes(x = group_var, y = average_quant, fill = group_var)) + 
        geom_bar(stat = "identity") +
        scale_fill_manual(values = okabe_ito_colors[1:length(levels(plot_data$group_var))]) +
        geom_text(aes(label = paste0(group_total, " / ", length(pk_quant_tab[, 1]))), vjust = -0.5, size = 2.5) +
        theme_minimal(base_size = 8) + 
        theme(legend.position = "none", 
              axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.title = element_text(size = 10)) +
        labs(
          title = paste0("Standard: ", standard_name, " (",peak_rt, ")"),
          x = as.character(group_var),
          y = "Average Proportion \nof Sample"
        ) 
      
      #add plot to list
      plot_list[[standard_index]] = quant_barplot
      
      #if no peaks found in SD for current standard
    } else {
      #make a dummy plot and put big N/A label on it
      dummy_data = data.frame(group_var = factor(group_types, levels = group_types), average_quant = 0)
      empty_barplot = ggplot(dummy_data, aes(x = group_var, y = average_quant)) +
        geom_blank() +
        annotate("text", x = median(1:length(group_types)), y = 0.5, label = "N/A", size = 6) +
        theme_minimal(base_size = 8) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 10)
        ) +
        labs(
          title = paste0("Standard: ", standard_name, " (", standard_rt, ")"),
          x = tools::toTitleCase(gsub("_", " ", as.character(group_var))),
          y = "Average Proportion \nof Sample"
        ) +
        expand_limits(y = 1)
      
      #add dummy plot to plot list
      plot_list[[standard_index]] = empty_barplot
    }
  }
  #combine and return plots 
  combined_plot = wrap_plots(plot_list, ncol = 4) + 
    plot_annotation(title = paste0(meta$site[1], " Average Proportion of Sample for External Standard Peaks"))
  return(combined_plot)
}

#Doe
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/doe_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(doe_quant_tab, doe_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/doe_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(doe_quant_tab, doe_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/doe_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(doe_quant_tab, doe_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Fre
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/fre_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(fre_quant_tab, fre_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/fre_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(fre_quant_tab, fre_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/fre_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(fre_quant_tab, fre_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Eas
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/eas_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(eas_quant_tab, eas_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/eas_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(eas_quant_tab, eas_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/eas_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(eas_quant_tab, eas_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Jen
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/jen_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(jen_quant_tab, jen_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/jen_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(jen_quant_tab, jen_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2020/jen_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(jen_quant_tab, jen_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")



#Tut
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/tut_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(tut_quant_tab, tut_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/tut_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(tut_quant_tab, tut_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/tut_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(tut_quant_tab, tut_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Pow
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/pow_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(pow_quant_tab, pow_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/pow_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(pow_quant_tab, pow_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/pow_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(pow_quant_tab, pow_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Lee
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/lee_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(lee_quant_tab, lee_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/lee_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(lee_quant_tab, lee_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/lee_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(lee_quant_tab, lee_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")

#Lit
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/lit_sample_type.png",
       plot = plot_subset_peaks_quantity_external_standards(lit_quant_tab, lit_meta, standards_data_table, "sample_type"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/lit_treatment.png",
       plot = plot_subset_peaks_quantity_external_standards(lit_quant_tab, lit_meta, standards_data_table, "treatment"),
       width = 15, height = 10, units = "in", bg = "white")
ggsave(filename = "G:/Work/Ash Chem Plots/6.24/quantity_by_standards/2021/lit_tree_size.png",
       plot = plot_subset_peaks_quantity_external_standards(lit_quant_tab, lit_meta, standards_data_table, "tree_size"),
       width = 15, height = 10, units = "in", bg = "white")





plot_labelled_chrom = function(peak_tab) {
  #find good sample in subset to plot
  sample_index = which.max(apply(peak_tab$tab, 1, function(row) sum(row != 0)))
  sample_name = rownames(peak_tab$tab)[sample_index]
  
  #store subset peak data
  rts = as.numeric(peak_tab$pk_meta["rt", ])
  ints = peak_tab$tab[sample_index, ]
  labels_df = data.frame(rt = rts[ints != 0], intensity = ints[ints != 0])
  
  #get samples chrom data
  chrom_list = get(attr(get(peak_tab$args$peak_list), "chrom_list"))
  chrom = chrom_list[[sample_index]]
  
  chrom_data_df = data.frame(
    time = as.numeric(rownames(chrom)),
    intensity = unname(chrom[, "280"])
  )
  
  #plot chrom with peak labels
  ggplot(chrom_data_df, aes(x = time, y = intensity)) +
    geom_line() +
    geom_text(data = labels_df,
              aes(x = rt, y = intensity * 1.05, label = round(rt, 2)),
              size = 3, vjust = 0, inherit.aes = FALSE) +
    labs(title = paste0("Chromatogram with Peaks: ", sample_name),
         x = "Retention Time", y = "Intensity") +
    theme_minimal()
}

plot_labelled_chrom(doe_pk_tab)




#fed the single chrom plotting method to ChatGPT and asked it to pick a 
#representative subset to cover all identified peaks and plot them overlayed 
#with peaks labelled so that the labels do not overlap
plot_labelled_chrom_overlay = function(peak_tab) {
  # Step 1: Extract info
  all_peak_rts = as.numeric(peak_tab$pk_meta["rt", ])
  tab = peak_tab$tab
  chrom_list = get(attr(get(peak_tab$args$peak_list), "chrom_list"))
  
  # Step 2: Create binary peak matrix
  coverage_matrix = tab != 0
  rownames(coverage_matrix) = rownames(tab)
  colnames(coverage_matrix) = all_peak_rts
  
  # Step 3: Select subset of chromatograms to cover all peaks
  peaks_remaining = colnames(coverage_matrix)
  selected_samples = c()
  while (length(peaks_remaining) > 0) {
    cover_counts = rowSums(coverage_matrix[, peaks_remaining, drop = FALSE])
    best_sample = names(which.max(cover_counts))
    selected_samples = c(selected_samples, best_sample)
    covered_peaks = names(which(coverage_matrix[best_sample, ]))
    peaks_remaining = setdiff(peaks_remaining, covered_peaks)
  }
  
  # Step 4: Gather all chrom data into a single data frame
  chrom_overlay_df = do.call(rbind, lapply(selected_samples, function(sample) {
    chrom = chrom_list[[sample]]
    data.frame(
      time = as.numeric(rownames(chrom)),
      intensity = unname(chrom[, "280"]),
      sample = sample
    )
  }))
  
  # Step 5: Identify best label positions for each peak
  label_positions = lapply(all_peak_rts, function(rt) {
    nearby_points = chrom_overlay_df[
      chrom_overlay_df$time >= rt - 0.05 &
        chrom_overlay_df$time <= rt + 0.05,
    ]
    if (nrow(nearby_points) > 0) {
      max_idx = which.max(nearby_points$intensity)
      peak = nearby_points[max_idx, ]
      return(data.frame(rt = peak$time, intensity = peak$intensity))
    } else {
      return(data.frame(rt = rt, intensity = 0))
    }
  })
  labels_df = do.call(rbind, label_positions)
  
  # Step 6: Plot all chromatograms overlaid, label each peak once
  ggplot(chrom_overlay_df, aes(x = time, y = intensity, color = sample)) +
    geom_line(alpha = 0.6) +
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(x = rt, y = intensity, label = round(rt, 2)),
      inherit.aes = FALSE,
      size = 2.7,                     # Smaller label text
      max.overlaps = Inf,            # Allow more labels
      box.padding = 0.4,             # More space around label
      point.padding = 0.25,
      segment.color = "gray60",
      segment.size = 0.3,
      min.segment.length = 0,
      force = 1.5,                   # Stronger repelling force
      force_pull = 0.1
    ) +
    labs(
      title = paste0(peak_tab$sample_meta$site[1], " All Identified Peaks Labelled"),
      x = "Retention Time",
      y = "Intensity"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2020/doe.png",
       plot = plot_labelled_chrom_overlay(doe_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2020/fre.png",
       plot = plot_labelled_chrom_overlay(fre_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2020/eas.png",
       plot = plot_labelled_chrom_overlay(eas_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2020/jen.png",
       plot = plot_labelled_chrom_overlay(jen_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")


ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2021/tut.png",
       plot = plot_labelled_chrom_overlay(tut_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2021/pow.png",
       plot = plot_labelled_chrom_overlay(pow_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2021/lee.png",
       plot = plot_labelled_chrom_overlay(lee_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")

ggsave(filename = "G:/Work/Ash Chem Plots/6.24/labelled_peaks/2021/lit.png",
       plot = plot_labelled_chrom_overlay(lit_pk_tab),
       width = 15, height = 10, units = "in", bg = "white")



###---Plotting Algorithmic Detetction by Treatment Type---###
# plot_standard_presence_rate_treatment(std_pres_tab, subset_meta)
# 
# This method uses ggplot2 to plot the detection rate of the 12 standards by treatment type.
#
# @param std_pres_tab - A Boolean matrix indicating the presence of the 12 standards in the treatments of the subset.
# @param subset_meta - A data frame containing the metadata for the susbet, including tree size, site, year, and species.
#
# @return - A combined plot displaying the detection rate of standards over the treatment types of the subset 
#
plot_standard_presence_rate_treatment_multi_alg = function(std_pres_tab_list, subset_meta) {
  #storing all unique treatments in the subset as a sorted list, excluding "Unkno" entries
  treatments = unlist(as.list(setdiff((mixedsort(unique(str_to_title(subset_meta[, 5])))), "Unkno")))
  
  #Creating an empty list to store ggplot2 plot objects created iteratively 
  plot_list = list() 
  
  #Iterate over the 12 standards (columns)
  for(i in 1:12) {
    #Creating data frame to store data for plotting
    standard_presence_rate_df = data.frame(Treatment = rep(treatments, each = 3), 
                                           Algorithm = rep(names(std_pres_tab_list), times = length(treatments)),
                                           Detection_Count = 0,
                                           Treatment_Total = 0)
    
    #Iterate over the treatments (rows)
    for(j in 1:length(std_pres_tab_list[[1]][,1])) { 
      #Find the treatment type of the current treatment
      treatment_index = which(unlist(treatments) == str_to_title(subset_meta[j, 5])) 
      
      #Iterate over 3 algorithms
      for(k in 1:3) {
        #Check if the current standard is present in the current treatment 
        if(std_pres_tab_list[[k]][j, i] == TRUE) { 
          #Increment the algorithmic detection count for the current treatment's treatment type
          standard_presence_rate_df[(treatment_index - 1) * 3 + k, "Detection_Count"] =
            standard_presence_rate_df[(treatment_index - 1) * 3 + k, "Detection_Count"] + 1
        }
        #Increment the total number of runs of the current treatment type by 1
        standard_presence_rate_df[(treatment_index - 1) * 3 + k, "Treatment_Total"] =
          standard_presence_rate_df[(treatment_index - 1) * 3 + k, "Treatment_Total"] + 1
      }
    }
    #Calculate the detection rates 
    standard_presence_rate_df$Detection_Rate = standard_presence_rate_df$Detection_Count / standard_presence_rate_df$Treatment_Total
    
    #Replace NA values with 0 in the detection rate column of the data frame 
    standard_presence_rate_df[is.na(standard_presence_rate_df)] = 0
    
    #Create and store a bar plot with detection rate on the Y axis and treatment type on the X axis 
    plot_list[[i]] = ggplot(standard_presence_rate_df, 
                            aes(x = factor(Treatment, levels = treatments), y = Detection_Rate, fill = Algorithm)) + 
      geom_col(position = "dodge", color = "black") +
      geom_text(aes(label = paste0("N = ", Treatment_Total, "\n", round(Detection_Rate, 2))), 
                position = position_dodge(width = 0.9), vjust = -0.2, size = 1.2) + 
      scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) + 
      labs(title = colnames(std_pres_tab_list[[1]])[i], x = "Tree Treatment", y = "Detection Rate", fill = "Algorithm") + 
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      theme(plot.title = element_text(size = 7, color = "#E69F00"),
            axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5), 
            axis.text.x = element_text(size = 3), axis.text.y = element_text(size = 5), 
            legend.position = "bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 5)) + 
      coord_cartesian(clip = "off")
  }
  #Combine all plots and add annotations with subset site, year, and tree species
  combined_plot = wrap_plots(plot_list) + 
    plot_annotation(title = paste0(subset_meta[j, 8], " ", subset_meta[j, 7],
                                   " - ", subset_meta[j, 3], " Ash | N = # of trees of respective treatment"))
  
  return(combined_plot)
}

#Get and store the three different standard presence tables for doe
doe_alg_det_lin = get_standard_presence_table(doe_pk_tab, standards_data_table)
doe_alg_det_ref_spec = get_standard_presence_table_spec_check(doe_pk_tab, standards_data_table)
doe_alg_det_opt = get_standard_presence_table_multiple(doe_pk_tab, standards_data_table)

#Store all three in a vector
doe_alg_det_list = list(Linear_RT_Search = doe_alg_det_lin, Reference_Spectra_Check = doe_alg_det_ref_spec, Optimized_RT_Search = doe_alg_det_opt)

plot_standard_presence_rate_treatment_multi_alg(doe_alg_det_list, doe_meta)







compare_identified_peaks = function(peak_tab, standards_data, standard_index) {
  #Call and store the algorithmic detection data for the given peak table
  peak_detection_rt_tab = std_detection_rt_table(peak_tab, standards_data)
  
  #Make a subset with all peaks identified as standard of interest
  std_subset_pk_det_rt_tab = peak_detection_rt_tab[peak_detection_rt_tab$Standard == standards_data[1, standard_index], ]
  #print to make sure we only have what we wanted
  #print(std_subset_pk_det_rt_tab)
  
  detected_sample_id = std_subset_pk_det_rt_tab$Sample_ID[which.max(std_subset_pk_det_rt_tab$Intensity)]
  empty_sample_id = peak_detection_rt_tab$Sample_ID[match(FALSE, peak_detection_rt_tab$Standard == standards_data[1, standard_index])]
  #print statements to make sure the Sample_IDs match criteria 
  print(detected_sample_id)
  print(empty_sample_id)
  
  #Use a get call on the args data from the peak list to store the chrom list
  pk_tab_chrom_list = get(peak_tab$args$chrom_list)
  
  detected_chrom = NULL
  empty_chrom = NULL
  
  #Loop through chrom list to find detected and empty chroms
  for(chrom in pk_tab_chrom_list) {
    if(attr(chrom, "sample_name") == detected_sample_id) {
      detected_chrom = chrom
    } else if(attr(chrom, "sample_name") == empty_sample_id) {
      empty_chrom = chrom
    } 
    #print(attr(chrom, "sample_name"))
    if(!is.null(detected_chrom) && !is.null(empty_chrom)) {
      break
    }
  }
  
  #print(detected_chrom)
  #print(empty_chrom)
  
  #combine intensity and time data in a data frame together for plotting the chroms for detected and empty
  detected_chrom_data = data.frame(time = as.numeric(attr(detected_chrom, "dimnames")[[1]]), intensity = as.numeric(detected_chrom[,'280']))
  empty_chrom_data = data.frame(time = as.numeric(attr(empty_chrom, "dimnames")[[1]]), intensity = as.numeric(empty_chrom[,'280']))
  
  #print(detected_chrom_data)
  #print(empty_chrom_data)
  
  #store the rt of the detected peak and it's intensity value
  detected_rt = std_subset_pk_det_rt_tab$Detected_RT[std_subset_pk_det_rt_tab$Sample_ID == detected_sample_id]
  detected_rt_intensity = std_subset_pk_det_rt_tab$Intensity[std_subset_pk_det_rt_tab$Sample_ID == detected_sample_id]
  #print(detected_rt)
  #print(detected_rt_intensity)
  
  #Define the window of the RT to plot over
  window_min = detected_rt - 2.5
  window_max = detected_rt + 2.5
  
  #Plot the intensity over time data for the detected and empty runs, label and theme properly 
  detected_chrom_plot = ggplot(detected_chrom_data, aes(x = time, y = intensity)) + geom_line(size = 0.185) +
    labs(x = "Time", y = "Intensity", title = paste0(detected_sample_id, " | ", standards_data[1, standard_index])) + 
    theme(plot.title = element_text(size = 7, color = "#E69F00"), 
          axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5), 
          axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5)) + xlim(window_min, window_max)
  
  empty_chrom_plot = ggplot(empty_chrom_data, aes(x = time, y = intensity)) + geom_line(size = 0.185) +
    labs(x = "Time", y = "Intensity", title = paste0(empty_sample_id, " | ", standards_data[1, standard_index])) + 
    theme(plot.title = element_text(size = 7, color = "#E69F00"), 
          axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5), 
          axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5)) + xlim(window_min, window_max)
  
  #Calculate the standard deviation range for the peak
  detected_rt_min = (detected_rt - as.double(standards_data[3, standard_index]))
  detected_rt_max = (detected_rt + as.double(standards_data[3, standard_index]))
  
  print(paste0(detected_rt, " ", detected_rt_intensity, " ", detected_rt_min, " ", detected_rt_max))
  #Create a data frame for plotting the peak data 
  detected_rt_df = data.frame(detected_rt, detected_rt_intensity, detected_rt_min, detected_rt_max) 
  
  #Add a point at the RT and intensity of the identified peak, with a value label,
  #and horizontal error bars to represent SD
  detected_chrom_plot = detected_chrom_plot + geom_point(data = detected_rt_df, aes(x = detected_rt, y = detected_rt_intensity),
                                               size = 1.5, shape = 23, stroke = 0.375, fill = NA, color = "#0072b2") + 
    geom_vline(xintercept = detected_rt, linetype = "dashed", color = "#D55E00", size = 0.25) +
    annotate("text", x = detected_rt, y = max(detected_chrom[,2]), 
             label = as.character(detected_rt), size = 2.5, color = "#D55E00", hjust = 1.4)
  chroms = list(detected_chrom_plot, empty_chrom_plot)
  
  chroms = wrap_plots(chroms, nrow = 2)
  
  return(chroms)
  
}

compare_identified_peaks(doe_pk_tab, standards_data_table, 1)
compare_identified_peaks(doe_pk_tab, standards_data_table, 2)
compare_identified_peaks(doe_pk_tab, standards_data_table, 3)
compare_identified_peaks(doe_pk_tab, standards_data_table, 4)
compare_identified_peaks(doe_pk_tab, standards_data_table, 5)
compare_identified_peaks(doe_pk_tab, standards_data_table, 6)
compare_identified_peaks(doe_pk_tab, standards_data_table, 7)
compare_identified_peaks(doe_pk_tab, standards_data_table, 8)
compare_identified_peaks(doe_pk_tab, standards_data_table, 9)
compare_identified_peaks(doe_pk_tab, standards_data_table, 10)
compare_identified_peaks(doe_pk_tab, standards_data_table, 11)
compare_identified_peaks(doe_pk_tab, standards_data_table, 12)



