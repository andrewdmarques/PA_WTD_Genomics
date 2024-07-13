################################################################################
# Functions for WTD Report Generator
################################################################################

get_old_data <- function(met1,old1){
  # Reformat the old data into the same format to be appended to the end of the new data.
  # Make a blank data frame
  col <- colnames(met1)
  old2 <- data.frame(matrix(NA, nrow = nrow(old1), ncol = length(col)))
  colnames(old2) <- col
  
  # Clean up the old data frame first.
  old1$result <- ifelse(old1$N1.Ct != 'Undetermined' & old1$N2.Ct != 'Undetermined', 'pos', 'neg')
  old1$AnimalStatus <- old1$Cause.of.death
  old1$AnimalStatus <- gsub('Hunter Harvested','Hunter / Trapper',old1$AnimalStatus)
  old1$AnimalStatus <- gsub('Road Kill','Roadkill',old1$AnimalStatus)
  old1$AnimalStatus <- gsub('Unkown','Undetermiend',old1$AnimalStatus)
  old1$Gender <- gsub('M','Male',old1$Gender)
  old1$Gender <- gsub('F','Female',old1$Gender)
  old1$County.Location <- toupper(old1$County.Location)
  
  # Populate the data frame.
  old2$VSP <- old1$VSP
  old2$IncidentDate <- old1$Date.of.collection
  old2$sars_cov_2 <- old1$result
  old2$n1_average_ct <- as.numeric(gsub('Undetermined',NA,old1$N1.Ct))
  old2$n2_average_ct <- as.numeric(gsub('Undetermined',NA,old1$N2.Ct))
  old2$seq_variant <- old1$Variant
  old2$Age <- gsub('Unknown','Undetermined',old1$Age)
  old2$CountyValue <- old1$County.Location
  old2$CWD <- 'neg'
  old2$Gender <- old1$Gender
  old2$Animalstatus <- old1$AnimalStatus
  
  
  # Combine together
  met2 <- rbind(met1,old2)
  
}

linearize <- function(data_frame){
  r <- rownames(data_frame)
  c <- colnames(data_frame)
  
  # Make a blank data frame.
  col <- c('dim1','dim2', 'value')
  linear <- data.frame(matrix(NA, nrow = length(r)*length(c), ncol = length(col)))
  colnames(linear ) <- col
  
  x <- 1
  y <- 1
  for(i in 1:length(linear$dim1)){
    linear$dim1[i] <- r[x]
    linear$dim2[i] <- c[y]
    linear$value[i] <- data_frame[x,y]
    y <- y + 1
    if(y > length(c)){
      y <- 1
      x <- x + 1
    }
  }
  return(linear)
}

################################################################################
# Functions for prevalence over time
################################################################################

# Function that determines the available dates -- prevalence over time 
get_date <- function(met2){
  # Determine the dates that should be included.
  date_incl1 <- unique(met2$date)
  date_incl2 <- seq.Date(from = min(date_incl1), to = max(date_incl1), by = "day")
  date_incl3 <- date_incl2
  # Get the date rounded to the first day of the week (Sunday).
  weekday_numbers = as.integer(format(date_incl3, "%w"))
  # Subtract the weekday number from the date to get the previous Sunday
  date_incl4 <- date_incl3 - weekday_numbers
  date_incl5 <- unique(date_incl4)
  
  return(date_incl5)
}

# Function that groups the data by age. -- prevalence over time
get_all <- function(prev2, group, dir_data_main,prefix,suffix){
  
  # Create a vector of numbers for the 4 groups that will be made. Not all will be needed for the analysis.
  numbers <- 1:4
  
  # Generate the strings
  file_group <- sapply(numbers, function(n) {
    paste0(dir_data_main, prefix, '_', suffix, '_',tolower(group),'-', n, '.csv')
  })
  
  # Make sure that all columns are populated
  prev2_temp <- subset(prev2,prev2$sars_cov_2 != '')
  
  
  # Make group 1 (non-CWD samples)
  gg1 <- prev2_temp
  group1 <- int_format(gg1,group,date_incl5)
  write.csv(group1,file_group[1],row.names = F)
  
  # Make group 2 (CWD samples)
  gg2 <- subset(prev2_temp,prev2_temp$Age == 'none')
  group2 <- int_format(gg2,group,date_incl5)
  write.csv(group2,file_group[2],row.names = F)
  
  # Make group 3 (all samples)
  gg3 <- subset(prev2_temp,prev2_temp$Age == 'nonenone')
  group3 <- int_format(gg3,group,date_incl5)
  write.csv(group3,file_group[3],row.names = F)
  
  # Make group 4 (all samples)
  gg4 <- subset(prev2_temp,prev2_temp$Age == 'nonenonenone')
  group4 <- int_format(gg4,group,date_incl5)
  write.csv(group4,file_group[4],row.names = F)
  
  return(file_group)
}

# Function that groups the data by age. -- prevalence over time
get_cwd <- function(prev2, group, dir_data_main,prefix,suffix){
  
  # Create a vector of numbers for the 4 groups that will be made. Not all will be needed for the analysis.
  numbers <- 1:4
  
  # Generate the strings
  file_group <- sapply(numbers, function(n) {
    paste0(dir_data_main, prefix, '_', suffix, '_',tolower(group),'-', n, '.csv')
  })
  
  # Make sure that all columns are populated
  prev2_temp <- subset(prev2,prev2$sars_cov_2 != '')
  
  # Make group 1 (non-CWD samples)
  # gg1 <- subset(prev2_temp,prev2_temp$CWD == 'neg')
  gg1 <- prev2_temp
  group1 <- int_format(gg1,group,date_incl5)
  write.csv(group1,file_group[1],row.names = F)
  
  # Make group 2 (CWD samples)
  gg2 <- prev2_temp
  group2 <- int_format(gg2,group,date_incl5)
  write.csv(group2,file_group[2],row.names = F)
  
  # Make group 3 (all samples)
  gg3 <- prev2_temp
  # gg3 <- subset(prev2_temp,prev2_temp$Age == 'Fawn')
  group3 <- int_format(gg3,group,date_incl5)
  write.csv(group3,file_group[3],row.names = F)
  
  # Make group 4 (all samples)
  gg4 <- prev2_temp
  # gg4 <- subset(prev2_temp,prev2_temp$Age == 'Yearling')
  group4 <- int_format(gg4,group,date_incl5)
  write.csv(group4,file_group[4],row.names = F)
  
  return(file_group)
}

# Function that formats the intermediate files that are saved as csv. -- prevalence over time
int_format <- function(gg1,group,date_incl5){
  tt1 <- gg1
  
  # Make a data frame of the correct dimensions.
  col <- c('date_repeat', 'mutant_label','counts','counts_total','percent')
  tt2 <- data.frame(matrix(NA, nrow = length(date_incl5), ncol = length(col)))
  colnames(tt2) <- col
  tt2$date_repeat <- date_incl5
  tt2$mutant_label <- group
  
  # Iterate through the data frame and populate it.
  for(ii in 1:nrow(tt2)){
    tt <- tt2$date_repeat[ii]
    aa <- subset(gg1,gg1$date == tt)
    bb <- table(aa$positive) # Make a table that has the number of positive and negative samples (note that bb[1] should be pos and bb[2] should be neg). This requires that$positive is a factor object with values of pos or neg. 
    
    # Populate the table.
    tt2$counts[ii] <- bb[1]
    tt2$counts_total[ii] <- bb[1] + bb[2]
    tt2$percent[ii] <- tt2$counts[ii] / tt2$counts_total[ii] * 100
  }
  
  # Remove the rows that had no data in them.
  tt3 <- tt2[!is.na(tt2$percent),]
  
  return(tt3)
}

# Function that linearizes data frames so that it is fit for plotting using ggplot.
# Linearize the data
linearize <- function(data_frame){
  r <- rownames(data_frame)
  c <- colnames(data_frame)
  # Make a blank data frame.
  col <- c('dim1','dim2', 'value')
  linear <- data.frame(matrix(NA, nrow = length(r)*length(c), ncol = length(col)))
  colnames(linear ) <- col
  
  x <- 1
  y <- 1
  for(i in 1:length(linear$dim1)){
    linear$dim1[i] <- r[x]
    linear$dim2[i] <- c[y]
    linear$value[i] <- data_frame[x,y]
    y <- y + 1
    if(y > length(c)){
      y <- 1
      x <- x + 1
    }
  }
  return(linear)
}

# Define function to prepare data for bayesian input. -- prevalence over time
bayesian_prep <- function(file_group, date_incl5,dir_data_main,prefix,suffix){
  a1 <- read.csv(file_group[1])
  b1 <- read.csv(file_group[2])
  c1 <- read.csv(file_group[3])
  d1 <- read.csv(file_group[4])
  
  # Get all of the dates in the date range. 
  date4 <- date_incl5
  date_length <- length(date4)
  date4 <- rep(date4, each = 2)
  date_real_repeat <- rep(date4, each = length(unique(a1$mutant_label)))
  date <- 1:date_length
  date <- rep(date, each = 2)
  date_repeat <- rep(date, each = length(unique(a1$mutant_label)))
  
  # Adapt the data frame to have all the dates.
  bay <- data.frame(matrix(NA, nrow = length(date_repeat), ncol = 6))
  colnames(bay) <- c("date_real","date", "rationale", "mutation", "category", "count")
  bay$date_real <- date_real_repeat
  bay$date <- date_repeat
  mutations <- rep(unique(a1$mutant_label), date_length)
  mutations <- rep(mutations, each = 2)
  bay$mutation <- mutations
  
  # Populate the data frame with group1.
  bay1 <- bay
  for(i in 1:length(bay1$date)){
    sub <- subset(a1, as.character(a1$date_repeat) == (bay1$date_real[i]))
    sub <- subset(sub, sub$mutant_label == bay1$mutation[i])
    if((i %% 2) == 0) {
      bay1$category[i] <- "FALSE"
      bay1$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay1$rationale[i] <- "group1"
    } else {
      bay1$category[i] <- "TRUE"
      bay1$count[i] <- sub$counts[1]
      bay1$rationale[i] <- "group1"
    }
  }
  bay1[is.na(bay1)] <- 0
  
  # Populate the data frame with group2.
  bay2 <- bay
  for(i in 1:length(bay2$date)){
    #i <- 1
    sub <- subset(b1, b1$date_repeat == bay2$date_real[i])
    sub <- subset(sub, sub$mutant_label == bay2$mutation[i])
    if((i %% 2) == 0) {
      bay2$category[i] <- "FALSE"
      bay2$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay2$rationale[i] <- "group2"
    } else {
      bay2$category[i] <- "TRUE"
      bay2$count[i] <- sub$counts[1]
      bay2$rationale[i] <- "group2"
    }
  }
  bay2[is.na(bay2)] <- 0
  
  # Populate the data frame with group3.
  bay3 <- bay
  for(i in 1:length(bay3$date)){
    #i <- 1
    sub <- subset(c1, c1$date_repeat == bay3$date_real[i])
    sub <- subset(sub, sub$mutant_label == bay3$mutation[i])
    if((i %% 2) == 0) {
      bay3$category[i] <- "FALSE"
      bay3$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay3$rationale[i] <- "group3"
    } else {
      bay3$category[i] <- "TRUE"
      bay3$count[i] <- sub$counts[1]
      bay3$rationale[i] <- "group3"
    }
  }
  bay3[is.na(bay3)] <- 0
  
  # Populate the data frame with group4.
  bay4 <- bay
  for(i in 1:length(bay4$date)){
    #i <- 1
    sub <- subset(d1, d1$date_repeat == bay4$date_real[i])
    sub <- subset(sub, sub$mutant_label == bay4$mutation[i])
    if((i %% 2) == 0) {
      bay4$category[i] <- "FALSE"
      bay4$count[i] <- sub$counts_total[1] - sub$counts[1]
      bay4$rationale[i] <- "group4"
    } else {
      bay4$category[i] <- "TRUE"
      bay4$count[i] <- sub$counts[1]
      bay4$rationale[i] <- "group4"
    }
  }
  bay4[is.na(bay4)] <- 0
  
  # Combine the data frames to one data frame.
  bay5 <- rbind(bay1, bay2, bay3, bay4)
  
  return(bay5)
}


# Define function for bayesian calculations. -- prevalence over time
bayesian_calculations <- function(bay5,group,date_incl5,dir_data_main,prefix,suffix){
  
  moi <- unique(bay5$mutation)
  # If this is to be done iterating through each of the groups specified in the 'mutation' category. For WTD SARS-CoV-2 this is ignored and run across all data.
  for(ll in 1:length(moi)){
    mutation_target <- moi[ll]
    
    # Reformat into the dat structure.
    bay6 <- subset(bay5, bay5$mutation == mutation_target)
    bay7 <- structure(bay6$count,
                      .Dim = c(2L, length(bay6$date)/8, 4L), # 8 comes from 4 different groups (group1, group2, group3, group4) and 2 categories (True and False); 2x4=8
                      .Dimnames = list(c("TRUE", "FALSE"), unique(bay6$date), 
                                       c("group1", "group2", "group3", "group4")))
    dat <- bay7
    #source('functions_prevalence_time.R')
    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
    
    # Uncomment these lines to run the model again.
    mod <- rstan::stan_model("model_prevalence_time.stan")
    fit<-runMutationStan(dat[,,'group1'],dat[,,'group2'],dat[,,'group3'],mod, nIter=200,nChain = 20)
    # Save the stan object.
    saveRDS(fit$stan, paste(dir_data_main, "/",prefix,"_",suffix,"_",tolower(group),"_bayesian.rds", sep = ""))
    fit_stan <- fit$stan
    
    # Plot the graph of the specific mutation with the confidence interval
    # Prepare the stan data object.
    mut_stan <- as.matrix(fit_stan)
    mut_stan <- invLogit(mut_stan)
    # Make the credible interval
    cri<-apply(mut_stan[,grep('means',colnames(mut_stan))],2,quantile,c(.025,.975))
    mut_stan <- apply(mut_stan[,grep('means',colnames(mut_stan))],2,mean)
    
    # Convert the week counts back to the week names.
    mut_stan2 <- data.frame(matrix(NA, nrow = length(mut_stan), ncol = 4))
    colnames(mut_stan2) <- c("Date", "Proportion", "2.5%", "97.5%")
    mut_stan2$Date <- as.POSIXct(date_incl5) 
    # Add the proportions
    for(i in 1:length(mut_stan)){
      mut_stan2$Proportion[i] <- mut_stan[i]
    }
    # Add the credible interval
    for(i in 1:(length(cri)/2)){
      mut_stan2$`2.5%`[i] <- cri[1,i]
      mut_stan2$`97.5%`[i] <- cri[2,i]
    }
    
    # Save the estimated proportion as .csv file.
    write.csv(mut_stan2, paste(dir_data_main, "/",prefix,"_",tolower(group), "_proportion-over-time",suffix,"_", ".csv", sep =""), quote=FALSE, row.names = F)
    write.csv(mut_stan2, paste0(dir_data_main,prefix,'_temp_prevalence_mut_stan2_',suffix,'.csv'), quote=FALSE, row.names = F)
    
    #XXX Added 2024-04-26 showing raw data.
    # Get the main data.
    temp_raw1 <- data.frame(t(dat[,,'group1']))
    rownames(temp_raw1) <- as.character(as.Date(mut_stan2$Date))
    colnames(temp_raw1) <- c('Positive','Negative')
    temp_raw2 <- linearize(temp_raw1)
    colnames(temp_raw2) <- c('YearMonth2','Positive','Count')
    temp_raw2$YearMonth2 <-  as.POSIXct(temp_raw2$YearMonth2)
    
    # Continue with plotting as before...
    plot_prevalence_stan <- ggplot(temp_raw2, aes(x = YearMonth2, y = Count, fill = Positive)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Negative" = "lightblue", "Positive" = "darkred")) +
      labs(x = "", y = "Counts", fill = "SARS-CoV-2 Positivity") +
      scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months", date_minor_breaks = "1 month", expand = c(0, 0)) +
      theme_classic()  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    plot_prevalence_stan
    
    # Plot the data for a smoothing curve.
    mut_stan2_perc <- mut_stan2
    mut_stan2_perc$Proportion <- mut_stan2_perc$Proportion * 100
    mut_stan2_perc$`2.5%` <- mut_stan2_perc$`2.5%` * 100
    mut_stan2_perc$`97.5%` <- mut_stan2_perc$`97.5%`  * 100
    plot_occurrence <- ggplot() +
      geom_line(data = mut_stan2_perc, aes(x=Date, y=Proportion, color="Proportion"), lwd = 1.01) +
      geom_ribbon(data = mut_stan2_perc, aes(x=Date, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1, fill = "black")+
      scale_color_manual(name = "", values = c("Prortion" = "black"))+
      scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months", date_minor_breaks = "1 month", expand = c(0, 0)) +
      theme_classic()  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")  +
      ylim(0, 50) +
      ylab("Percent")
    plot_occurrence
    
    # Break up the data so that there are two groups grouped by date. date1 is samples Jan to Jun, date2 is Jul to Dec.
    temp_raw3 <- temp_raw2
    temp_raw3$month  <- format(temp_raw3$YearMonth2, "%m")
    temp_raw3$month <- as.numeric(temp_raw3$month)
    # Assign the date group for jan to jun or jul to dec
    date_group <- c('Jan_to_Jun','Jul_to_Dec')
    temp_raw3$date_group <- ifelse(temp_raw3$month >= 1 & temp_raw3$month <= 6, date_group[1], date_group[2])
    
    # Determine the total number of positive and negative results for each of the date groups.
    unique_dates <- unique(temp_raw2$YearMonth2)
    col <- c('date_group', 'Positive','Negative')
    temp_raw4 <- data.frame(matrix(NA, nrow = length(date_group), ncol = length(col)))
    colnames(temp_raw4) <- col
    temp_raw4$date_group <- date_group
    # Populate the data frame.
    for(qq in 1:nrow(temp_raw4)){
      qq1 <- temp_raw4$date_group[qq]
      qq2 <- subset(temp_raw3,temp_raw3$date_group == qq1)
      qq_pos <- subset(qq2,qq2$Positive=='Positive')
      qq_neg <- subset(qq2,qq2$Positive=='Negative')
      temp_raw4$Positive[qq] <- sum(qq_pos$Count)
      temp_raw4$Negative[qq] <- sum(qq_neg$Count)
    }
    temp_raw4$Total <- temp_raw4$Positive + temp_raw4$Negative
    prop_test_date_group <- prop.test(x = temp_raw4$Positive, 
                                      n = temp_raw4$Total)
    
    
    plot_occurrence   
  }
  
  # Save the plotting files to be plotted later.
  file_plot_prevalence_stan <- paste(dir_data_main,prefix,"_", tolower(group), "_plot_prevalence_raw","_",suffix,".rds", sep ="")
  saveRDS(plot_prevalence_stan, file = file_plot_prevalence_stan) 
  file_plot_occurrence <- paste(dir_data_main,prefix,"_",suffix,"_", tolower(group), "_plot_prevalence_smooth","_",suffix,".rds", sep ="")
  saveRDS(plot_occurrence, file = file_plot_occurrence) 
  
  file_plot_bayes <- c(file_plot_prevalence_stan,file_plot_occurrence)
  return(file_plot_bayes)
  # Combine all of the bayesian plots into one pdf
  # combine_pdf(paste("_Bayesian", sep = ""))
}

get_cwd_positivity <- function(met1,dir_data_main,prefix,suffix){
  sta1 <- subset(met1, met1$VSP != 'VSP30247')
  sta1 <- subset(sta1,sta1$sars_cov_2 != '')
  
  # Curate the data for the positive values.
  sta1_p <- subset(sta1,sta1$CWD == 'pos')
  # Get the date format.
  sta1_p$date_exact <- as.Date(sta1_p$IncidentDate)
  # Remove any rows that are missing metadata.
  sta1_p <- sta1_p[!is.na(sta1_p$date_exact),]
  # Get the date rounded to the first day of the week (Sunday).
  weekday_numbers = as.integer(format(sta1_p$date_exact, "%w"))
  # Subtract the weekday number from the date to get the previous Sunday
  sta1_p$date <- sta1_p$date_exact - weekday_numbers
  # Get the table for positivity by dates.
  sta2_p <- data.frame(table(sta1_p$date,sta1_p$sars_cov_2))
  colnames(sta2_p) <- c('YearMonth2','Positive','Count')
  sta2_p$YearMonth2 <- as.POSIXct(sta2_p$YearMonth2)
  sta2_p$Positive <- as.character(sta2_p$Positive)
  sta2_p$Positive <- gsub('pos','Positive',sta2_p$Positive)
  sta2_p$Positive <- gsub('neg','Negative',sta2_p$Positive)
  
  # Curate the data for the negative values.
  sta1_n <- subset(sta1,sta1$CWD == 'neg')
  # Get the date format.
  sta1_n$date_exact <- as.Date(sta1_n$IncidentDate)
  # Remove any rows that are missing metadata.
  sta1_n <- sta1_n[!is.na(sta1_n$date_exact),]
  # Get the date rounded to the first day of the week (Sunday).
  weekday_numbers = as.integer(format(sta1_n$date_exact, "%w"))
  # Subtract the weekday number from the date to get the previous Sunday
  sta1_n$date <- sta1_n$date_exact - weekday_numbers
  # Get the table for positivity by dates.
  sta2_n <- data.frame(table(sta1_n$date,sta1_n$sars_cov_2))
  colnames(sta2_n) <- c('YearMonth2','Positive','Count')
  sta2_n$YearMonth2 <- as.POSIXct(sta2_n$YearMonth2)
  sta2_n$Positive <- as.character(sta2_n$Positive)
  sta2_n$Positive <- gsub('pos','Positive',sta2_n$Positive)
  sta2_n$Positive <- gsub('neg','Negative',sta2_n$Positive)
  
  # Subset the negative values so that only those weeks that have positive values are plotted.
  sta3_p <- sta2_p
  sta_date_keep <- unique(sta3_p$YearMonth2)
  sta_max_date <- max(sta2_n$YearMonth2)
  sta3_n <- subset(sta2_n, sta2_n$YearMonth2 %in% c(sta_date_keep))
  # If the min date is missing then add it in. (this is adding the date from the negative to the date to the positive only!!!)
  sta_min_date <- min(sta3_p$YearMonth2)
  if(min(sta3_n$YearMonth2) > sta_min_date){
    # Data frame to be appended to make a minimum date for the surveillance samples.
    col <- c('YearMonth2','Positive','Count')
    sta_temp <- data.frame(matrix(NA, nrow = 2, ncol = length(col)))
    colnames(sta_temp) <- col
    sta_temp$YearMonth2 <- sta_min_date
    sta_temp$Positive <- c('Negative','Positive')
    sta_temp$Count <- 0
    
    sta3_n <- rbind(sta3_n,sta_temp)
  }
  # If the max date is missing then add it in. (this is adding the date from the positive to the date to the negative only!!!)
  # Data frame to be appended to make a minimum date for the surveillance samples.
  col <- c('YearMonth2','Positive','Count')
  sta_temp <- data.frame(matrix(NA, nrow = 2, ncol = length(col)))
  colnames(sta_temp) <- col
  sta_temp$YearMonth2 <- sta_max_date
  sta_temp$Positive <- c('Negative','Positive')
  sta_temp$Count <- 0
  sta3_n <- rbind(sta3_n,sta_temp)
  sta3_p <- rbind(sta3_p,sta_temp)
  
  
  # Plot prevalence
  plot_prevalence_p <- ggplot(sta3_p, aes(x = YearMonth2, y = Count, fill = Positive)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Negative" = "lightblue", "Positive" = "darkred")) +
    labs(x = "", y = "CWD Positive", fill = "SARS-CoV-2 Positivity") +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months", date_minor_breaks = "1 month", expand = c(0, 0)) +
    theme_classic()  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  # plot_prevalence_p
  # Plot prevalence
  plot_prevalence_n <- ggplot(sta3_n, aes(x = YearMonth2, y = Count, fill = Positive)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Negative" = "lightblue", "Positive" = "darkred")) +
    labs(x = "", y = "CWD Negative", fill = "SARS-CoV-2 Positivity") +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months", date_minor_breaks = "1 month", expand = c(0, 0)) +
    theme_classic()  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  # plot_prevalence_n
  
  # Save the plots for later.
  file_plot_prevalence_n <- paste(dir_data_main,prefix,"_",suffix,"_", tolower(group), "_plot_prevalence_n.rds", sep ="")
  saveRDS(plot_prevalence_n, file = file_plot_prevalence_n) 
  file_plot_prevalence_p <- paste(dir_data_main,prefix,"_",suffix,"_plot_prevalence_p.rds", sep ="")
  saveRDS(plot_prevalence_p, file = file_plot_prevalence_p) 
  
  # Run stats test on the data to determine if there is a significant differnece across the times where collection was done.
  temp1 <- subset(sta3_p,sta3_p$Positive == 'Positive')
  sta_cwd_pos <- sum(temp1$Count)
  sta_cwd_tot <- sum(sta3_p$Count)
  temp2 <- subset(sta3_n,sta3_n$Positive == 'Positive')
  sta_noncwd_pos <- sum(temp2$Count)
  sta_noncwd_tot <- sum(sta3_n$Count)
  sta_prop_test <- prop.test(x = c(sta_cwd_pos, sta_noncwd_pos), 
                             n = c(sta_cwd_tot,  sta_noncwd_tot))
  file_sta_prop_test <- paste(dir_data_main,prefix,"_",suffix,"_sta_prop_test.rds", sep ="")
  saveRDS(sta_prop_test, file = file_sta_prop_test) 
  
  # Get the table for the CWD status.
  col <- c('CWD', 'Counts','Counts_Positive','Percent_Positive')
  sta4 <- data.frame(matrix(NA, nrow = 2, ncol = length(col)))
  colnames(sta4) <- col
  sta4$CWD <- c('Negative','Positive')
  sta4$Counts <- c(sta_noncwd_tot,sta_cwd_tot)
  sta4$Counts_Positive <- c(sta_noncwd_pos,sta_cwd_pos)
  sta4$Percent_Positive <- c(as.character(round(sta_noncwd_pos/sta_noncwd_tot*100,2)),as.character(round(sta_cwd_pos/sta_cwd_tot*100,2)))
  file_sta_table <- paste(dir_data_main,prefix,"_",suffix,"_sta_table.rds", sep ="")
  saveRDS(sta4, file = file_sta_table) 
  
  file_sta <- c(file_plot_prevalence_n,file_plot_prevalence_p,file_sta_prop_test,file_sta_table)
  return(file_sta)
}

# Define function for bayesian calculations. -- prevalence over time
bayesian_calculations <- function(bay5,group,date_incl5,dir_data_main,prefix,suffix){
  
  
  # Function that linearizes data frames so that it is fit for plotting using ggplot.
  # Linearize the data
  linearize <- function(data_frame){
    r <- rownames(data_frame)
    c <- colnames(data_frame)
    # Make a blank data frame.
    col <- c('dim1','dim2', 'value')
    linear <- data.frame(matrix(NA, nrow = length(r)*length(c), ncol = length(col)))
    colnames(linear ) <- col
    
    x <- 1
    y <- 1
    for(i in 1:length(linear$dim1)){
      linear$dim1[i] <- r[x]
      linear$dim2[i] <- c[y]
      linear$value[i] <- data_frame[x,y]
      y <- y + 1
      if(y > length(c)){
        y <- 1
        x <- x + 1
      }
    }
    return(linear)
  }
  
  moi <- unique(bay5$mutation)
  for(ll in 1:length(moi)){
    mutation_target <- moi[ll]
    
    # Reformat into the dat structure.
    bay6 <- subset(bay5, bay5$mutation == mutation_target)
    bay7 <- structure(bay6$count,
                      .Dim = c(2L, length(bay6$date)/8, 4L), # 8 comes from 4 different groups (group1, group2, group3, group4) and 2 categories (True and False); 2x4=8
                      .Dimnames = list(c("TRUE", "FALSE"), unique(bay6$date), 
                                       c("group1", "group2", "group3", "group4")))
    dat <- bay7
    #source('functions_prevalence_time.R')
    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
    
    # Uncomment these lines to run the model again.
    mod <- rstan::stan_model("model_prevalence_time.stan")
    fit<-runMutationStan(dat[,,'group1'],dat[,,'group2'],dat[,,'group3'],mod, nIter=200,nChain = 20)
    # Save the stan object.
    saveRDS(fit$stan, paste(dir_data_main, "/",prefix,"_",suffix,"_",tolower(group),"_bayesian.rds", sep = ""))
    fit_stan <- fit$stan
    
    # Plot the graph of the specific mutation with the confidence interval
    # Prepare the stan data object.
    mut_stan <- as.matrix(fit_stan)
    mut_stan <- invLogit(mut_stan)
    # Make the credible interval
    cri<-apply(mut_stan[,grep('means',colnames(mut_stan))],2,quantile,c(.025,.975))
    mut_stan <- apply(mut_stan[,grep('means',colnames(mut_stan))],2,mean)
    
    # Convert the week counts back to the week names.
    mut_stan2 <- data.frame(matrix(NA, nrow = length(mut_stan), ncol = 4))
    colnames(mut_stan2) <- c("Date", "Proportion", "2.5%", "97.5%")
    mut_stan2$Date <- as.POSIXct(date_incl5) 
    # Add the proportions
    for(i in 1:length(mut_stan)){
      mut_stan2$Proportion[i] <- mut_stan[i]
    }
    # Add the credible interval
    for(i in 1:(length(cri)/2)){
      mut_stan2$`2.5%`[i] <- cri[1,i]
      mut_stan2$`97.5%`[i] <- cri[2,i]
    }
    
    # Save the estimated proportion as .csv file.
    write.csv(mut_stan2, paste(dir_data_main, "/",prefix,"_",suffix,"_", tolower(group), "_proportion-over-time.csv", sep =""), quote=FALSE, row.names = F)
    write.csv(mut_stan2, paste0(dir_data_main,prefix,'_temp_prevalence_mut_stan2_',suffix,'.csv'), quote=FALSE, row.names = F)
    
    #XXX Added 2024-04-26 showing raw data.
    # Get the main data.
    temp_raw1 <- data.frame(t(dat[,,'group1']))
    rownames(temp_raw1) <- as.character(as.Date(mut_stan2$Date))
    colnames(temp_raw1) <- c('Positive','Negative')
    temp_raw2 <- linearize(temp_raw1)
    colnames(temp_raw2) <- c('YearMonth2','Positive','Count')
    temp_raw2$YearMonth2 <-  as.POSIXct(temp_raw2$YearMonth2)
    
    # Continue with plotting as before...
    plot_prevalence_stan <- ggplot(temp_raw2, aes(x = YearMonth2, y = Count, fill = Positive)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Negative" = "lightblue", "Positive" = "darkred")) +
      labs(x = "", y = "All WTD Counts", fill = "SARS-CoV-2 Positivity") +
      scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months", date_minor_breaks = "1 month", expand = c(0, 0)) +
      theme_classic()  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    plot_prevalence_stan
    
    # Plot the data for a smoothing curve.
    mut_stan2_perc <- mut_stan2
    mut_stan2_perc$Proportion <- mut_stan2_perc$Proportion * 100
    mut_stan2_perc$`2.5%` <- mut_stan2_perc$`2.5%` * 100
    mut_stan2_perc$`97.5%` <- mut_stan2_perc$`97.5%`  * 100
    plot_occurrence <- ggplot() +
      geom_line(data = mut_stan2_perc, aes(x=Date, y=Proportion, color="Proportion"), lwd = 1.01) +
      geom_ribbon(data = mut_stan2_perc, aes(x=Date, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.1, fill = "black")+
      scale_color_manual(name = "", values = c("Prortion" = "black"))+
      scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months", date_minor_breaks = "1 month", expand = c(0, 0)) +
      theme_classic()  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")  +
      ylim(0, 50) +
      ylab("Percent")
    plot_occurrence
    
    # Break up the data so that there are two groups grouped by date. date1 is samples Jan to Jun, date2 is Jul to Dec.
    temp_raw3 <- temp_raw2
    temp_raw3$month  <- format(temp_raw3$YearMonth2, "%m")
    temp_raw3$month <- as.numeric(temp_raw3$month)
    # Assign the date group for jan to jun or jul to dec
    date_group <- c('Jan_to_Jun','Jul_to_Dec')
    temp_raw3$date_group <- ifelse(temp_raw3$month >= 1 & temp_raw3$month <= 6, date_group[1], date_group[2])
    
    # Determine the total number of positive and negative results for each of the date groups.
    unique_dates <- unique(temp_raw2$YearMonth2)
    col <- c('date_group', 'Positive','Negative')
    temp_raw4 <- data.frame(matrix(NA, nrow = length(date_group), ncol = length(col)))
    colnames(temp_raw4) <- col
    temp_raw4$date_group <- date_group
    # Populate the data frame.
    for(qq in 1:nrow(temp_raw4)){
      qq1 <- temp_raw4$date_group[qq]
      qq2 <- subset(temp_raw3,temp_raw3$date_group == qq1)
      qq_pos <- subset(qq2,qq2$Positive=='Positive')
      qq_neg <- subset(qq2,qq2$Positive=='Negative')
      temp_raw4$Positive[qq] <- sum(qq_pos$Count)
      temp_raw4$Negative[qq] <- sum(qq_neg$Count)
    }
    temp_raw4$Total <- temp_raw4$Positive + temp_raw4$Negative
    prop_test_date_group <- prop.test(x = temp_raw4$Positive, 
                                      n = temp_raw4$Total)
    text <- paste0(
      'Two sided t test performed on the groups ',date_group[1],' and ', date_group[2],'.\n',
      'p=',as.character(prop_test_date_group$p.value),'\n',
      'Proportion of ',date_group[1],' = ',as.character(prop_test_date_group$estimate[1]),' (n = ',as.character(temp_raw4$Total[1]),')\n',
      'Proportion of ',date_group[2],' = ',as.character(prop_test_date_group$estimate[2]),' (n = ',as.character(temp_raw4$Total[2]),')\n'
    )
    write.csv(text, paste(dir_data_main,prefix,"_",suffix,"_", tolower(group), "_t-test-for-time-groups.txt", sep =""), quote=FALSE, row.names = F)
    
    
    
    plot_occurrence   
    
    # Save the plotting files to be plotted later.
    file_plot_prevalence_stan <- paste(dir_data_main,prefix,"_",suffix,"_", tolower(group), "_plot_prevalence_stan.rds", sep ="")
    saveRDS(plot_prevalence_stan, file = file_plot_prevalence_stan) 
    file_plot_occurrence <- paste(dir_data_main,prefix,"_",suffix,"_", tolower(group), "_plot_occurrence.rds", sep ="")
    saveRDS(plot_occurrence, file = file_plot_occurrence) 
    
    file_plot_bayes <- c(file_plot_prevalence_stan,file_plot_occurrence)
    return(file_plot_bayes)
  }
}

# Function that formats the data into the correct format. -- prevalence over time
format_deer_prev <- function(met1){
  temp1 <- met1
  
  # Get the date format.
  temp1$date_exact <- as.Date(temp1$IncidentDate)
  # Remove any rows that are missing metadata.
  temp1 <- temp1[!is.na(temp1$date_exact),]
  # Get the date rounded to the first day of the week (Sunday).
  weekday_numbers = as.integer(format(temp1$date_exact, "%w"))
  # Subtract the weekday number from the date to get the previous Sunday
  temp1$date <- temp1$date_exact - weekday_numbers
  
  # Determine if it is a positive sample or not.
  temp1$positive <- temp1$sars_cov_2
  # temp1$positive <- 'neg'
  # # Set positive to TRUE where nq_average_ct is not NA
  # temp1$positive[!is.na(temp1$n1_average_ct)] <- 'pos'
  temp1$positive <- factor(temp1$positive,levels = c('pos','neg'))
  # temp1 <- temp1[!grepl("extraction error", temp1$notes), ] # Remove the samples with extraction errors.
  temp1 <- temp1[!grepl("VSP35", temp1$VSP), ] # Remove the experimental samples.
  
  met2 <- temp1
  return(met2)
}


# Original functions below.

# Earlier paper. -- prevalence over time
runStan<-function(countTab,mod,iter=2000){
  stan_sample <- rstan::sampling(
    mod,
    data=list(
      counts=countTab,
      nTime=ncol(countTab),
      nLineage=nrow(countTab)
    ),
    iter=iter,
    chains=50,
    control=list(max_treedepth=15),
  )
  return(stan_sample)
}

# Earlier paper. -- prevalence over time
calcCredInt<-function(stan_sample,quants=c(.025,.975),names=NULL){
  mat<-as.matrix(stan_sample)
  props<-mat[,grepl('^props\\[',colnames(mat))]
  meanCrI<-function(xx)c(mean(xx),quantile(xx,quants))
  meanLowerUpper<-apply(props,2,meanCrI)
  convertToMat<-function(xx){
    row<-as.numeric(sub('[^0-9]+([0-9]+),[0-9]+.*','\\1',names(xx)))
    col<-as.numeric(sub('[^0-9]+[0-9]+,([0-9]+).*','\\1',names(xx)))
    tapply(xx,list(row,col),c)
  }
  means<-convertToMat(meanLowerUpper[1,])
  lower<-convertToMat(meanLowerUpper[2,])
  upper<-convertToMat(meanLowerUpper[3,])
  rownames(means)<-rownames(lower)<-rownames(upper)<-names
  return(list('mean'=means,'upper'=upper,'lower'=lower))
}

calcDensity<-function(stan_sample,names=NULL){
  mat<-as.matrix(stan_sample)
  group3s<-mat[,grepl('^group3Change\\[',colnames(mat))]
  group2s<-mat[,grepl('^group2Change\\[',colnames(mat))]
  minMax<-range(cbind(group3s,group2s))
  group2Density<-apply(group2s,2,density,from=minMax[1],to=minMax[2])
  group3Density<-apply(group3s,2,density,from=minMax[1],to=minMax[2])
  names(group2Density)<-names(group3Density)<-names
  return(list('group2'=group2Density,'group3'=group3Density))
}

# Earlier paper. -- prevalence over time
plotMeanUpperLower<-function(means,upper,lower,countTab,baseDate,cols=NULL){
  #https://sashamaps.net/docs/resources/20-colors/
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i')
  for(ii in 1:nrow(means)){
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],rownames(countTab),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
}

# Earlier paper. -- prevalence over time
plotBars<-function(countTab,baseDate,cols=NULL,showLegend=TRUE,showSum=TRUE){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[rownames(propTab)])
  uniqDate<-data.frame('week'=which(counts>0),'rdate'=baseDate+(which(counts>0)-1)*7)
  dnar::slantAxis(1,uniqDate$week,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),location=.5,cex=.6)
  if(showSum)dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  if(showLegend)legend('bottom',names(cols),fill=cols,inset=-.33,ncol=ceiling(nrow(countTab)/2),xpd=NA,cex=.8)
}

# Earlier paper. -- prevalence over time
runStan2<-function(countTab,group3Tab,group2Tab,mod,iter=2000){#,hospitalTab
  #should do more double checking
  stan_sample <- rstan::sampling(
    mod,
    data=list(
      counts=countTab,
      #hospital=hospitalTab,
      group3=group3Tab,
      group2=group2Tab,
      nTime=ncol(countTab),
      nLineage=nrow(countTab)
    ),
    iter=iter,
    chains=50,
    thin=2,
    control=list(max_treedepth=15),
    pars=c('means','propsGroup3','propsGroup2'),
    include=FALSE
  )
  return(stan_sample)
}

# Earlier paper. -- prevalence over time
plotIndivStan<-function(means,upper,lower,countTab,baseDate,cols=NULL,nCol=5){
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  counts<-apply(countTab,2,sum)
  barCol<-rev(grey.colors(max(counts)+1))
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  #https://sashamaps.net/docs/resources/20-colors/
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21','7/1/21'))
  for(ii in 1:nrow(means)){
    plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-4,4),ylim=c(0,1.2),las=1,xlab='',bty='l',xaxt='n',yaxs='i',xaxs='i',yaxt='n',xaxt='n')
    prettyY<-c(0,.5,1)
    if(ii %% nCol==1)axis(2,las=1,prettyY)
    else axis(2,prettyY,rep('',length(prettyY)))
    if(ii > nrow(countTab)-nCol) dnar::slantAxis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.2,location=.5)
    else axis(1,(prettyDates-baseDate)/7+1,rep('',length(prettyDates)),tcl=-.2)
    if(ii==(2*nCol+1))mtext('Estimated proportion',2,line=2.2,at=1)
    title(main=rownames(countTab)[ii],line=-2,cex=.9)
    rect(1:ncol(propTab)-.5,0,1:ncol(propTab)+.5,propTab[ii,],col=barCol[counts+1],border=NA)
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s77',cols[rownames(countTab)[ii]]))
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
    box()
  }
  abline(h=0)
}

# Earlier paper. -- prevalence over time
plotIndivDense<-function(dense,cols=NULL,xlim=range(exp(dense[[1]]$x)),ylab='Fold change'){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  #https://sashamaps.net/docs/resources/20-colors/
  yMax<-max(sapply(dense,function(xx)max(xx$y)))
  for(ii in 1:length(dense)){
    print(names(dense)[ii])
    plot(1,1,type='n',xlim=xlim,ylim=c(0,yMax),las=1,xlab='',bty='n',xaxt='n',yaxs='i',xaxs='i',yaxt='n',xaxt='n',log='x')
    polygon(exp(dense[[ii]]$x),dense[[ii]]$y,col=cols[names(dense)[ii]])
    title(main=names(dense)[ii],line=0,cex=.9)
    dnar::logAxis(1,axisVals=c(-2,0,2,4))
    abline(v=1,lty=2)
    #if(ii==(2*nCol+1))mtext(2,line=2.2,at=par('usr')[4])
  }
  text(grconvertX(.015,'ndc','user'),grconvertY(.5,'ndc','user'),'Estimated posterior probability',xpd=NA,cex=1.2,srt=90)
  text(grconvertX(.5,'ndc','user'),grconvertY(.015,'ndc','user'),ylab,xpd=NA,cex=1.2)
}

#deprecated
plotStan<-function(stan_sample,countTab,baseDate,cols=NULL){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  mat<-as.matrix(stan_sample)
  props<-mat[,grepl('^props\\[',colnames(mat))]
  meanCrI<-function(xx)c(mean(xx),quantile(xx,c(.025,.975)))
  meanLowerUpper<-apply(props,2,meanCrI)
  convertToMat<-function(xx){
    row<-as.numeric(sub('[^0-9]+([0-9]+),[0-9]+.*','\\1',names(xx)))
    col<-as.numeric(sub('[^0-9]+[0-9]+,([0-9]+).*','\\1',names(xx)))
    tapply(xx,list(row,col),c)
  }
  means<-convertToMat(meanLowerUpper[1,])
  lower<-convertToMat(meanLowerUpper[2,])
  upper<-convertToMat(meanLowerUpper[3,])
  #https://sashamaps.net/docs/resources/20-colors/
  par(mar=c(3,3,.5,4),tcl=-.2,mgp=c(2,.3,0))
  layout(matrix(c(1,0,2,0),ncol=1),height=c(1,.01,1,.15))
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i')
  for(ii in 1:nrow(means)){
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],rownames(countTab),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
  #
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[rownames(propTab)])
  uniqDate<-data.frame('week'=which(counts>0),'rdate'=baseDate+(which(counts>0)-1)*7)
  dnar::slantAxis(1,uniqDate$week,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),location=.5,cex=.6)
  dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  legend('bottom',names(cols),fill=cols,inset=-.33,ncol=ceiling(nrow(countTab)/2),xpd=NA,cex=.8)
  invisible(return(list(means,lower,upper)))
}

# Earlier paper. -- prevalence over time
runStan3<-function(countTab,group3Tab,group2Tab,greek,mod,iter=2000){
  #should do more double checking
  if(length(greek)!=nrow(countTab))stop('Sublineage groupings not same length as count table')
  abundantGreek<-names(table(greek)[table(greek)>1])
  greekIds<-structure(1:(length(abundantGreek)+1),.Names=c('__BASE__',abundantGreek))
  dat<-list(
    counts=countTab,
    group3=group3Tab,
    group2=group2Tab,
    nTime=ncol(countTab),
    nLineage=nrow(countTab),
    nSublineageGroup=max(greekIds),
    sublineageGroup=greekIds[ifelse(greek %in% names(greekIds),greek,'__BASE__')]
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=iter,
    chains=50,
    thin=2,
    control=list(max_treedepth=15),
    pars=c('means'),
    include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'greek'=greekIds,'lineage'=rownames(countTab)))
}

# Earlier paper. -- prevalence over time
calcDensity2<-function(stan_sample,names=NULL,greekNames=NULL){
  mat<-as.matrix(stan_sample)
  group3s<-mat[,grepl('^group3Change\\[',colnames(mat))]
  group2s<-mat[,grepl('^group2ChangeWithSub\\[',colnames(mat))]
  greeks<-mat[,grepl('^group2Change\\[',colnames(mat))]
  minMax<-range(cbind(group3s,group2s,greeks))
  group2Density<-apply(group2s,2,density,from=minMax[1],to=minMax[2])
  group3Density<-apply(group3s,2,density,from=minMax[1],to=minMax[2])
  greekDensity<-apply(greeks,2,density,from=minMax[1],to=minMax[2])
  names(group2Density)<-names(group3Density)<-names
  names(greekDensity)<-greekNames
  return(list('group2'=group2Density,'group3'=group3Density,'greek'=greekDensity))
}


logit<-function(xx)log(xx)-log(1-xx)
invLogit<-function(xx)1/(1+exp(-xx))

# Earlier paper. -- prevalence over time
runMutationStan<-function(countTab,group3Tab,group2Tab,mod,nChain=50,nIter=2000){
  dat<-list(
    counts=countTab[1,],
    nCounts=apply(countTab,2,sum),
    group3=group3Tab[1,],
    nGroup3=apply(group3Tab,2,sum),
    group2=group2Tab[1,],
    nGroup2=apply(group2Tab,2,sum),
    nTime=ncol(countTab)
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=nIter,
    chains=nChain,
    thin=2,
    control=list(max_treedepth=15)
    #pars=c('means'),
    #include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'mod'=mod))
}

################################################################################
# Functions for anatomy of the lymph node
################################################################################

# Function that cleans the data. -- anatomy of the lymph node
##############################################################################

# Function that cleans the data.
get_clean_data <- function(comp1){
  # Extract only the lines from the data that we need (VSP35111 to VSP35173, only the odd ones)
  # Function to check if the number part of the string is odd and within the specified range
  valid_vsp <- paste0('VSP',as.character(seq(from = 35111, to = 35173, by = 1)))
  comp2 <- subset(comp1, comp1$VSP %in% valid_vsp)
  
  # Remove any samples that were not prepped (there were some blank VSPs to skip over).
  comp2 <- comp2[!is.na(comp2$n1_batch),]
  comp2 <- subset(comp2,comp2$n1_batch != '')
  
  # When uncommented, this section counts half positives as positive.
  for(ii in 1:nrow(comp2)){
    if(is.na(comp2$n1_average_ct[ii])){
      if(!is.na(comp2$n1_half_positive[ii])){
        comp2$n1_average_ct[ii] <- comp2$n1_half_positive[ii]
      }
    }
  }
  comp2$ct <- comp2$n1_average_ct
  comp4 <- comp2
  
  # Split up the notes section from the data.
  comp4$sample <- NA
  comp4$node_num <- NA
  comp4$buffer <- NA
  comp4$location <- NA
  for(ii in 1:nrow(comp4)){
    split <- unlist(strsplit(comp4$notes[ii], "_")) 
    comp4$sample[ii] <- split[1]
    comp4$node_num[ii] <- split[2]
    comp4$buffer[ii] <- split[3]
    comp4$location[ii] <- split[4]
  }
  comp5 <- comp4 
  
  return(comp5)
}

# Function that determines the statistics for each location. -- anatomy of the lymph node
get_stat <- function(dir_data,prefix,suffix,nod5,name){
  
  # Plot the data by location.
  colfunc <- colorRampPalette(c('blue', 'grey90','green','brown'))
  picked_colors <- colfunc(length(unique(nod5$sample)))
  p1 <- ggplot(nod5, aes(x = location, y = ct)) +
    geom_boxplot() +
    geom_jitter(aes(color = sample), width = 0.2, size = 1.5) +  # Color the dots based on 'sample'
    scale_color_manual(values = picked_colors) +  # Use custom colors
    theme_minimal() +
    labs(title = paste0("Ct by Location: ",name),
         x = "Location",
         y = "Ct",
         color = "Sample")  # Add label for the color legend
  p1
  ggsave(paste0(dir_data,'/',prefix,'_',suffix,'_ct-by-location_',tolower(name),'.pdf'), plot = p1, width = 170 / 25.4, height = 100 / 25.4)
  
  
  # Plot the data by sample.
  colfunc <- colorRampPalette(c('red','black','orange','grey90','deeppink3'))
  picked_colors <- colfunc(length(unique(nod5$location)))
  p2 <- ggplot(nod5, aes(x = sample, y = ct)) +
    geom_boxplot() +
    geom_jitter(aes(color = location), width = 0.2, size = 1.5) +  # Color the dots based on 'location'
    scale_color_manual(values = picked_colors) +  # Use custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
    labs(title = paste0("Ct by Sample: ", name),
         x = "Sample",
         y = "Ct",
         color = "Location")  # Add label for the color legend
  p2
  ggsave(paste0(dir_data,'/',prefix,'_',suffix,'_ct-by-sample_',tolower(name),'.pdf'),plot = p2,width = 170 / 25.4, height = 100 / 25.4)
  
  
  # Determine the averages and the standard deviation.
  nod_summary <- nod5 %>%
    group_by(location) %>%
    summarise(
      average = mean(ct, na.rm = TRUE),
      std_dev = sd(ct, na.rm = TRUE)
    )
  
  # Determine if any of the averages are outliers
  find_outliers <- function(numbers) {
    # Calculate mean and standard deviation
    mean_value <- mean(numbers)
    sd_value <- sd(numbers)
    
    # Calculate Z-scores
    z_scores <- (numbers - mean_value) / sd_value
    
    # Identify outliers
    outliers <- numbers[abs(z_scores) > 3]
    
    return(outliers)
  }
  outliers_average <- find_outliers(nod_summary$average)
  paste0('Number of significantly different averages: ',as.character(length(outliers_average)))
  outliers_std_dev <- find_outliers(nod_summary$std_dev)
  paste0('Number of significantly different standard deviations: ',as.character(length(outliers_std_dev)))
  
  # Run an ANOVA test 
  nod5$location <- as.factor(nod5$location)
  
  # Run the ANOVA test
  result <- aov(ct ~ location, data = nod5)
  summary(result)
  
  return(p2)
}

# Function that gets the normalized values for the Ct, where they are normalized against the average for that sample. -- anatomy of the lymph node
get_normalized <- function(nod5){
  nod7 <- nod5
  nod7$ct_original <- nod7$ct
  
  for(ii in 1:nrow(nod7)){
    tt <- nod7$sample[ii]
    temp1 <- subset(nod7,nod7$sample == tt)
    temp1 <- temp1[!is.na(temp1$ct_original),]
    nod7$ct[ii] <- nod7$ct_original[ii] - mean(temp1$ct_original)
  }
  
  return(nod7)
}

# Define the directory and file name pattern -- anatomy of the lymph node
get_summary <- function(dir_data,prefix,suffix){
  # Create the pattern to match files
  pattern <- paste0("^", prefix, "_", suffix, "_")
  
  # List all files in the directory
  all_files <- list.files(path = dir_data, pattern = pattern, full.names = TRUE)
  
  # Filter out files that contain 'summary'
  filtered_files <- all_files[!grepl("summary", all_files)]
  filtered_files <- filtered_files[grepl(".pdf", filtered_files)]
  
  # Determine the output file.
  output_file <- paste0(dir_data,'/', prefix, "_", suffix, "_summary.pdf")
  
  # Find PDF files containing "subject" in their name
  pdf_files <- filtered_files
  if (length(pdf_files) == 0) {
    print("No PDF files containing 'subject' in their name found in the directory.")
    return(NULL)
  }
  # Combine them with the other one
  pdf_combine(pdf_files, output = output_file)
  # Combine the pdf files.
  pdf_length(output_file)
  
  print("PDF files stitched successfully!")
}

# Function that cleans the data.
get_clean_data_noda <- function(noda1){
  # Extract only the lines from the data that we need (VSP35001 to VSP35110, only the odd ones)
  # Function to check if the number part of the string is odd and within the specified range
  is_valid_vsp <- function(vsp) {
    num <- as.numeric(sub("VSP", "", vsp))
    num >= 35001 && num <= 35110 && num %% 2 == 1
  }
  # Filter dat1 to create dat2
  noda2 <- noda1 %>%
    filter(sapply(VSP, is_valid_vsp))
  # Remove any samples that were not prepped (there were some blank VSPs to skip over).
  noda2 <- noda2[!is.na(noda2$n1_batch),]
  noda2 <- subset(noda2,noda2$n1_batch != '')
  
  # When uncommented, this section counts half positives as positive.
  for(ii in 1:nrow(noda2)){
    if(is.na(noda2$n1_average_ct[ii])){
      if(!is.na(noda2$n1_half_positive[ii])){
        noda2$n1_average_ct[ii] <- noda2$n1_half_positive[ii]
      }
    }
  }
  
  # Format the data into a data frame so that it has columns sample, location, ct, relative ct from center.
  noda3 <- noda2
  noda3 <- noda3[ , which(names(noda3) %in% c('VSP','notes','n1_average_ct'))]
  
  noda4 <- noda3
  
  # Split the 'notes' column
  split_notes <- strsplit(noda4$notes, "_")
  
  # Determine the maximum number of parts any string is split into
  max_parts <- max(sapply(split_notes, length))
  
  # Create a new data frame from the split strings
  # Pad shorter rows with NA
  split_df <- do.call(rbind, lapply(split_notes, function(x) {
    length(x) <- max_parts
    return(as.data.frame(t(x)))
  }))
  
  # Rename the columns (optional, but recommended for clarity)
  names(split_df) <- paste0("notes_part_", 1:max_parts)
  
  # Combine with original data frame
  noda5 <- cbind(noda4, split_df)
  
  # Rename the first columns as sample, location, buffer.
  names(noda5)[names(noda5) == "notes_part_1"] <- "sample"
  names(noda5)[names(noda5) == "notes_part_2"] <- "location_code"
  names(noda5)[names(noda5) == "notes_part_3"] <- "buffer"
  names(noda5)[names(noda5) == "n1_average_ct"] <- "ct"
  
  # Reclassify the location to the correct position
  locations_decode <- data.frame(
    location_code = c('A', 'B', 'C', 'D', 'E', 'F', 'G'),
    location = c('blood', 'adipose', 'efferent', 'distal_medial','cranial_medial', 'cranial_lateral', 'center')
    # location = c('blood', 'adipose', 'efferent', 'bottom_left', 'top_left', 'top_right', 'center')
  )
  
  # Assign the data to be found in the table.
  lookup_tab_locations <- locations_decode$location
  # Give the data names.
  names(lookup_tab_locations) <- locations_decode$location_code
  # Search a vector of names and outputs the associated data.
  noda5$location <- lookup_tab_locations[noda5$location_code]
  
  # Assign the location as a factored object.
  noda5$location <- factor(noda5$location, levels = c( 'blood', 'cranial_medial', 'center', 'distal_medial', 'cranial_lateral', 'efferent', 'adipose'))
  return(noda5)
}


################################################################################
# Functions for early brms model
################################################################################

# Make a dataframe for all the comparisons to be made. -- brms model
multiple_comparisons <- function(model,title){
  dat_sex <- model %>% spread_draws(seasonspring) 
  dat_season <- model %>% spread_draws(b_season[season])
  
  
  
  model %>%
    spread_draws(b_seasonwinter[condition,Intercept]) -> dat1
  model %>% spread_draws(b_seasonwinter) -> dat_t1
  dat_t1$condition <- 'Alpha'
  model %>% spread_draws(b_variantDelta) -> dat_t2
  dat_t2$condition <- 'Delta'
  model %>% spread_draws(b_variantOmicron) -> dat_t3
  dat_t3$condition <- 'Omicron'
  dat1 <- rbind(dat_t1,dat_t2,dat_t3)
  comps <- c('Delta-Alpha','Omicron-Alpha','Delta-Omicron')
  col <- c('comparison', 'difference')
  dat2 <- data.frame(matrix(NA, nrow = 0, ncol = length(col)))
  colnames(dat2) <- col
  # Record the comparisons for all of the groups
  for(ii in 1:length(comps)){
    v1 <- sub("-.*", "", comps[ii])
    v2 <- sub(".*-", "", comps[ii])
    temp1 <- subset(dat1,dat1$condition == v1)
    temp2 <- subset(dat1,dat1$condition == v2)
    for(jj in 1:length(temp1$condition)){
      dat2_temp <- data.frame(matrix(NA, nrow = 1, ncol = length(col)))
      colnames(dat2_temp) <- col
      dat2_temp$comparison[1] <- comps[ii]
      dat2_temp$difference[1] <- temp1$r_variant[jj] - temp2$r_variant[jj]
      dat2 <- rbind(dat2,dat2_temp)
    }
  }
  # Plot
  cols <- c("#998ad1", "grey80",  "#bf3434")
  comparison <- ggplot(dat2, aes(x = difference, fill = comparison)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = cols)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    labs(y = 'Density', x = 'Expected Difference in Ct Value')+
    guides(fill=guide_legend(title="Comparison"))
  comparison
  
  pdf(title,width = 4, height = 6)
  par(mfrow = c(3,1), mai = c(0.4, 1, 0.4, 1))
  cols <- c("#998ad1","#bf3434","grey80")
  for(ii in 1:length(comps)){
    t1 <- comps[ii]
    temp <- subset(dat2,dat2$comparison == t1)
    lower <- as.numeric(quantile(temp$difference,probs=c(.025,.975))[1])
    upper <- as.numeric(quantile(temp$difference,probs=c(.025,.975))[2])
    dens <- density(temp$difference)
    xlabel <- ''
    if(ii == length(comps)){xlabel <- 'Difference in Ct Effect Size'}
    plot(dens, frame = FALSE, col = cols[ii], main = t1, xaxs = 'i',yaxs = 'i',xlim = c(-4,1),ylim = c(0,2.6),xlab = xlabel) 
    polygon(dens, col = cols[ii])
    abline(v = 0, col="black", lwd=2, lty=2)
    abline(v = lower, col=cols[ii], lwd=1, lty=1)
    abline(v = upper, col=cols[ii], lwd=1, lty=1)
    # abline(v = mean(temp$difference), col=cols[ii], lwd=4, lty=1)
  }
  dev.off()
  
}

# Plots the spline over time by machine. -- brms model
plot_by_machine <- function(m1_2,title){
  temp <- m1_2$data %>%as_tibble() %>% expand(days_from_variant_appearance = 0:200, variant = unique(variant), machine = unique(dat_ADO$machine), specimen_upper = unique(dat_ADO$specimen_upper))
  preds <-  add_epred_draws(newdata = temp, object = m1_2)
  # preds <- preds[ , -which(names(preds) %in% c("specimen_upper"))]
  
  meanCrI<-function(xx,groupings,lower=.1,upper=.9)list('mean'=tapply(xx,groupings,mean),'lower'=tapply(xx,groupings,quantile,lower),'upper'=tapply(xx,groupings,quantile,upper))
  means<-meanCrI(preds$.epred,preds[,c('variant','days_from_variant_appearance','machine','specimen_upper')])
  
  machineCol<-structure(dnar::rainbow.lab(length(unique(dat$machine)),alpha=.6),.Names=unique(dat$machine))
  pdf(paste0(title),width=8.5,height=11)
  yy <- 1
  par(mfrow=c(6,3),las=1,mar=c(4,4,1,.5),mgp=c(2,.5,0),tcl=-.3)
  for(machine in c(names(machineCol))){ # Do not include the machines 5 and 6 because they do not span all variants
    xx <- 1
    for(ii in rownames(means[[1]])){
      if(machine=='All') thisDat<-as.data.frame(dat[dat$variant==ii,])
      else thisDat<-as.data.frame(dat[dat$variant==ii&dat$machine==machine,])
      plot(1,1,type='n',xlim=c(0,200),ylim=range(dat$ct),xlab='Time after wave start',ylab='Ct',main=paste(machine,ii))
      if(nrow(thisDat)==0)next()
      points(thisDat$days_from_variant_appearance,thisDat$ct,bg=machineCol[thisDat$machine],col=NA,pch=21,cex=.6)
      polygon(c(0:200,200:0),c(means[[2]][ii,,machine,2],rev(means[[3]][ii,,machine,2])),border=NA,col='#FF000033')
      lines(0:200,means[[1]][ii,,machine,2],col='#FF5555')
      xx <- xx + 1  
    }
    yy <- yy + 1
  }
  dev.off()
}


################################################################################
# Functions for bayesian county positivity map
################################################################################

# Function -- county positivity map
runStan<-function(countTab,mod,iter=2000){
  stan_sample <- rstan::sampling(
    mod,
    data=list(
      counts=countTab,
      nTime=ncol(countTab),
      nLineage=nrow(countTab)
    ),
    iter=iter,
    chains=50,
    control=list(max_treedepth=15),
  )
  return(stan_sample)
}

# Function -- county positivity map
calcCredInt<-function(stan_sample,quants=c(.025,.975),names=NULL){
  mat<-as.matrix(stan_sample)
  props<-mat[,grepl('^props\\[',colnames(mat))]
  meanCrI<-function(xx)c(mean(xx),quantile(xx,quants))
  meanLowerUpper<-apply(props,2,meanCrI)
  convertToMat<-function(xx){
    row<-as.numeric(sub('[^0-9]+([0-9]+),[0-9]+.*','\\1',names(xx)))
    col<-as.numeric(sub('[^0-9]+[0-9]+,([0-9]+).*','\\1',names(xx)))
    tapply(xx,list(row,col),c)
  }
  means<-convertToMat(meanLowerUpper[1,])
  lower<-convertToMat(meanLowerUpper[2,])
  upper<-convertToMat(meanLowerUpper[3,])
  rownames(means)<-rownames(lower)<-rownames(upper)<-names
  return(list('mean'=means,'upper'=upper,'lower'=lower))
}

# Function -- county positivity map
calcDensity<-function(stan_sample,names=NULL){
  mat<-as.matrix(stan_sample)
  group3s<-mat[,grepl('^group3Change\\[',colnames(mat))]
  group2s<-mat[,grepl('^group2Change\\[',colnames(mat))]
  minMax<-range(cbind(group3s,group2s))
  group2Density<-apply(group2s,2,density,from=minMax[1],to=minMax[2])
  group3Density<-apply(group3s,2,density,from=minMax[1],to=minMax[2])
  names(group2Density)<-names(group3Density)<-names
  return(list('group2'=group2Density,'group3'=group3Density))
}

# Function -- county positivity map
plotMeanUpperLower<-function(means,upper,lower,countTab,baseDate,cols=NULL){
  #https://sashamaps.net/docs/resources/20-colors/
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i')
  for(ii in 1:nrow(means)){
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],rownames(countTab),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
}

# Function -- county positivity map
plotBars<-function(countTab,baseDate,cols=NULL,showLegend=TRUE,showSum=TRUE){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[rownames(propTab)])
  uniqDate<-data.frame('week'=which(counts>0),'rdate'=baseDate+(which(counts>0)-1)*7)
  dnar::slantAxis(1,uniqDate$week,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),location=.5,cex=.6)
  if(showSum)dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  if(showLegend)legend('bottom',names(cols),fill=cols,inset=-.33,ncol=ceiling(nrow(countTab)/2),xpd=NA,cex=.8)
}

# Function -- county positivity map
runStan2<-function(countTab,group3Tab,group2Tab,mod,iter=2000){#,hospitalTab
  #should do more double checking
  stan_sample <- rstan::sampling(
    mod,
    data=list(
      counts=countTab,
      #hospital=hospitalTab,
      group3=group3Tab,
      group2=group2Tab,
      nTime=ncol(countTab),
      nLineage=nrow(countTab)
    ),
    iter=iter,
    chains=50,
    thin=2,
    control=list(max_treedepth=15),
    pars=c('means','propsGroup3','propsGroup2'),
    include=FALSE
  )
  return(stan_sample)
}

# Function -- county positivity map
plotIndivStan<-function(means,upper,lower,countTab,baseDate,cols=NULL,nCol=5){
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  counts<-apply(countTab,2,sum)
  barCol<-rev(grey.colors(max(counts)+1))
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  #https://sashamaps.net/docs/resources/20-colors/
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21','7/1/21'))
  for(ii in 1:nrow(means)){
    plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-4,4),ylim=c(0,1.2),las=1,xlab='',bty='l',xaxt='n',yaxs='i',xaxs='i',yaxt='n',xaxt='n')
    prettyY<-c(0,.5,1)
    if(ii %% nCol==1)axis(2,las=1,prettyY)
    else axis(2,prettyY,rep('',length(prettyY)))
    if(ii > nrow(countTab)-nCol) dnar::slantAxis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.2,location=.5)
    else axis(1,(prettyDates-baseDate)/7+1,rep('',length(prettyDates)),tcl=-.2)
    if(ii==(2*nCol+1))mtext('Estimated proportion',2,line=2.2,at=1)
    title(main=rownames(countTab)[ii],line=-2,cex=.9)
    rect(1:ncol(propTab)-.5,0,1:ncol(propTab)+.5,propTab[ii,],col=barCol[counts+1],border=NA)
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s77',cols[rownames(countTab)[ii]]))
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
    box()
  }
  abline(h=0)
}

# Function -- county positivity map
plotIndivDense<-function(dense,cols=NULL,xlim=range(exp(dense[[1]]$x)),ylab='Fold change'){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  #https://sashamaps.net/docs/resources/20-colors/
  yMax<-max(sapply(dense,function(xx)max(xx$y)))
  for(ii in 1:length(dense)){
    print(names(dense)[ii])
    plot(1,1,type='n',xlim=xlim,ylim=c(0,yMax),las=1,xlab='',bty='n',xaxt='n',yaxs='i',xaxs='i',yaxt='n',xaxt='n',log='x')
    polygon(exp(dense[[ii]]$x),dense[[ii]]$y,col=cols[names(dense)[ii]])
    title(main=names(dense)[ii],line=0,cex=.9)
    dnar::logAxis(1,axisVals=c(-2,0,2,4))
    abline(v=1,lty=2)
    #if(ii==(2*nCol+1))mtext(2,line=2.2,at=par('usr')[4])
  }
  text(grconvertX(.015,'ndc','user'),grconvertY(.5,'ndc','user'),'Estimated posterior probability',xpd=NA,cex=1.2,srt=90)
  text(grconvertX(.5,'ndc','user'),grconvertY(.015,'ndc','user'),ylab,xpd=NA,cex=1.2)
}

#deprecated -- county positivity map
plotStan<-function(stan_sample,countTab,baseDate,cols=NULL){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  mat<-as.matrix(stan_sample)
  props<-mat[,grepl('^props\\[',colnames(mat))]
  meanCrI<-function(xx)c(mean(xx),quantile(xx,c(.025,.975)))
  meanLowerUpper<-apply(props,2,meanCrI)
  convertToMat<-function(xx){
    row<-as.numeric(sub('[^0-9]+([0-9]+),[0-9]+.*','\\1',names(xx)))
    col<-as.numeric(sub('[^0-9]+[0-9]+,([0-9]+).*','\\1',names(xx)))
    tapply(xx,list(row,col),c)
  }
  means<-convertToMat(meanLowerUpper[1,])
  lower<-convertToMat(meanLowerUpper[2,])
  upper<-convertToMat(meanLowerUpper[3,])
  #https://sashamaps.net/docs/resources/20-colors/
  par(mar=c(3,3,.5,4),tcl=-.2,mgp=c(2,.3,0))
  layout(matrix(c(1,0,2,0),ncol=1),height=c(1,.01,1,.15))
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i')
  for(ii in 1:nrow(means)){
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],rownames(countTab),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
  #
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[rownames(propTab)])
  uniqDate<-data.frame('week'=which(counts>0),'rdate'=baseDate+(which(counts>0)-1)*7)
  dnar::slantAxis(1,uniqDate$week,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),location=.5,cex=.6)
  dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  legend('bottom',names(cols),fill=cols,inset=-.33,ncol=ceiling(nrow(countTab)/2),xpd=NA,cex=.8)
  invisible(return(list(means,lower,upper)))
}

# Function -- county positivity map
runStan3<-function(countTab,group3Tab,group2Tab,greek,mod,iter=2000){
  #should do more double checking
  if(length(greek)!=nrow(countTab))stop('Sublineage groupings not same length as count table')
  abundantGreek<-names(table(greek)[table(greek)>1])
  greekIds<-structure(1:(length(abundantGreek)+1),.Names=c('__BASE__',abundantGreek))
  dat<-list(
    counts=countTab,
    group3=group3Tab,
    group2=group2Tab,
    nTime=ncol(countTab),
    nLineage=nrow(countTab),
    nSublineageGroup=max(greekIds),
    sublineageGroup=greekIds[ifelse(greek %in% names(greekIds),greek,'__BASE__')]
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=iter,
    chains=50,
    thin=2,
    control=list(max_treedepth=15),
    pars=c('means'),
    include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'greek'=greekIds,'lineage'=rownames(countTab)))
}

# Function -- county positivity map
calcDensity2<-function(stan_sample,names=NULL,greekNames=NULL){
  mat<-as.matrix(stan_sample)
  group3s<-mat[,grepl('^group3Change\\[',colnames(mat))]
  group2s<-mat[,grepl('^group2ChangeWithSub\\[',colnames(mat))]
  greeks<-mat[,grepl('^group2Change\\[',colnames(mat))]
  minMax<-range(cbind(group3s,group2s,greeks))
  group2Density<-apply(group2s,2,density,from=minMax[1],to=minMax[2])
  group3Density<-apply(group3s,2,density,from=minMax[1],to=minMax[2])
  greekDensity<-apply(greeks,2,density,from=minMax[1],to=minMax[2])
  names(group2Density)<-names(group3Density)<-names
  names(greekDensity)<-greekNames
  return(list('group2'=group2Density,'group3'=group3Density,'greek'=greekDensity))
}

logit<-function(xx)log(xx)-log(1-xx)
invLogit<-function(xx)1/(1+exp(-xx))

# Function -- county positivity map
runMutationStan<-function(countTab,group3Tab,group2Tab,mod,nChain=50,nIter=2000){
  dat<-list(
    counts=countTab[1,],
    nCounts=apply(countTab,2,sum),
    group3=group3Tab[1,],
    nGroup3=apply(group3Tab,2,sum),
    group2=group2Tab[1,],
    nGroup2=apply(group2Tab,2,sum),
    nTime=ncol(countTab)
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=nIter,
    chains=nChain,
    thin=2,
    control=list(max_treedepth=15)
    #pars=c('means'),
    #include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'mod'=mod))
}

logit<-function(xx)log(xx)-log(1-xx)
invLogit<-function(xx)1/(1+exp(-xx))
meanCrI<-function(xx,quants=c(.025,.975))c(mean(xx),quantile(xx,quants))

# Function -- county positivity map
runCountyStan<-function(countTab,adjacencyTab,mod,regions=rep(1,nrow(countTab)),nChain=50,nIter=10000){
  if(any(colnames(adjacencyTab)!=rownames(adjacencyTab))||any(rownames(countTab)!=rownames(adjacencyTab)))stop('Mismatch between counts and adjacency')
  regionId<-structure(1:length(unique(regions)),.Names=sort(unique(regions)))
  dat<-list(
    nCounty=nrow(countTab),
    nRegion=length(unique(regions)),
    regions=regionId[regions],
    counts=apply(countTab,1,sum),
    positives=countTab[,'TRUE'],
    adjacency=adjacencyTab
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=nIter,
    chains=nChain,
    thin=2
  )
  return(list('stan'=stan_sample,'dat'=dat,'mod'=mod,'tab'=countTab,'region'=regionId))
}

# Function -- county positivity map
insetScale<-function(breaks,col,insetPos=c(.025,.015,.04,.25),main='',offset=1e-3,at=NULL,labels=NULL,cex=1,labXOffset=0,labYOffset=0){
  if(length(breaks)!=length(col)+1)stop('Number of breaks must be one more than colors')
  insetPos<-c(graphics::grconvertY(insetPos[1],'nfc','user'),graphics::grconvertX(insetPos[2],'nfc','user'),graphics::grconvertY(insetPos[3],'nfc','user'),graphics::grconvertX(insetPos[4],'nfc','user'))
  breakPos<-((breaks)-(min(breaks)))/max((breaks)-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  #add a bit of offset to avoid pdf viewers displaying breaks between exact rectangle border meeting
  offsetPos<-breakPos[-1]+c(rep(offset*diff(range(breakPos)),length(breakPos)-2),0)
  graphics::rect(breakPos[-length(breakPos)],insetPos[1],offsetPos,insetPos[3],col=col,xpd=NA,border=NA)
  graphics::rect(insetPos[2],insetPos[1],insetPos[4],insetPos[3],xpd=NA)
  if(is.null(at)){
    at<-pretty(breaks)
    at<-at[at<=max(breaks)&at>=min(breaks)]
  }
  if(is.null(labels))labels<-at
  convertPos<-(at-(min(breaks)))/((max(breaks))-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  graphics::segments(convertPos,insetPos[1],convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.1,xpd=NA)
  graphics::text(convertPos+labXOffset*diff(insetPos[c(2,4)]),insetPos[1]-diff(insetPos[c(1,3)])*.175+labYOffset*diff(insetPos[c(1,3)]),labels,xpd=NA,adj=c(.5,1),cex=.85*cex)
  graphics::text(mean(insetPos[c(2,4)]),insetPos[3]+diff(insetPos[c(1,3)])*.45,main,xpd=NA,adj=c(.5,0),cex=cex)
  invisible(NULL)
}

# Function -- county positivity map
fillDown<-function(x,emptyStrings=c(NA,''),errorIfFirstEmpty=TRUE){
  #depending on %in% to catch NAs if necessary
  isEmpty<-x %in% emptyStrings
  if(isEmpty[1]&errorIfFirstEmpty)stop(simpleError('First value empty'))
  #if first is empty and we don't want errors then have to just fill down from it anyway
  isEmpty[1]<-FALSE
  ids<-1:length(x)
  ids[isEmpty]<-0
  ids<-cummax(ids)
  return(x[ids])
}

# Function -- county positivity map
get_geo_year <- function(date_start,date_stop,met4,deer,file_plot_geo_year,plot_title){
  props <- readRDS(paste0(dir_data_main,prefix,'_temp_county_props_',suffix,'.rds'))
  
  # Subset to have just the date of interest. 
  geo1 <- subset(met4,met4$IncidentDate >= date_start)
  geo1 <- subset(geo1,geo1$IncidentDate < date_stop)
  
  # Determine the proportion positive by county.
  geo2 <- data.frame(table(geo1$CountyValue,geo1$sars_cov_2))
  rr <- unique(as.character(geo2$Var1))
  cc <- c('name',unique(as.character(geo2$Var2)))
  
  # Make a blank data frame
  geo3 <- data.frame(matrix(NA, nrow = length(rr), ncol = length(cc)))
  colnames(geo3) <- cc
  geo3$name <- rr
  
  # Populate the data frame.
  for(kk in 1:nrow(geo3)){
    tt <- geo3$name[kk]
    temp1 <- subset(geo2,geo2$Var1 == tt)
    for(ll in 1:nrow(temp1)){
      geo3[kk,ll+1] <- temp1$Freq[ll]
    }
  }
  
  # Determine the total and proprtion positive
  geo3$total <- geo3$neg + geo3$pos
  geo3$prop <- geo3$pos / geo3$total
  
  geo3$name <- paste0('pennsylvania,',tolower(geo3$name))
  
  # Give the props data frame the correct proportions based off the selected date range.
  props1 <- props
  for(ii in 1:length(props[1,])){
    props1[1,ii] <- 0
    # Determine which name we are working with.
    name1 <- names(props1[1,])[ii]
    temp1 <- subset(geo3,geo3$name == name1)
    if(nrow(temp1)>0){
      props1[1,ii] <- temp1$prop[1]
      
      if(props1[1,ii] == 0){
        props1[1,ii] <- 0.001
      }
    }
    # print(props1[1,ii])
  }
  
  # Variables needed to plot: adjacency, counts, penn, labels, select
  # Get the adjacency data 
  #https://www.census.gov/geographies/reference-files/2010/geo/county-adjacency.html
  adj<-read.table('county_adjacency.txt',sep='\t')
  colnames(adj)<-c('county','id','neighbor','neighborId')
  adj$county<-fillDown(adj$county)
  adj<-adj[grepl(', PA$',adj$county)&grepl(', PA$',adj$neighbor),]
  adj$simple<-sub(' County, PA','',adj$county)
  adj$neighborSimple<-sub(' County, PA','',adj$neighbor)
  if(any(!deer$County.Location %in% adj$simple))stop('Unknown county')
  adjacencyAll<-table(adj$simple,adj$neighborSimple)
  diag(adjacencyAll)<-0
  adj<-adj[adj$simple %in% deer$County.Location &adj$neighborSimple %in% deer$County.Location,]
  adjacency<-table(adj$simple,adj$neighborSimple)
  diag(adjacency)<-0
  
  # Prepare the necessary files required for generating the figures.
  counts<-table(deer$County.Location,deer$pos)
  
  penn <- maps::map("county","Pennsylvania",plot=FALSE)
  penn$prettyName<-sub('pennsylvania,','',penn$name)
  substring(penn$prettyName,1,1)<-toupper(substring(penn$prettyName,1,1))
  
  colnames(props1)<-sprintf('pennsylvania,%s',tolower(rownames(adjacency)))
  if(any(!colnames(props1) %in% penn$names))stop('Unknown county name')
  
  cols<-1+penn$names %in% colnames(props1)
  # breaks<-seq(0,ceiling(max(props1[1,])*10)/10,.005) # This will select it so that it only goes as high as it needs to go
  breaks<-seq(0,1,.005)
  cuts<-cut(props1[1,],breaks)
  cuts2<-cut(props1[2,],breaks)
  nCut<-length(levels(cuts))
  nCutTop<-0;nCutBottom<-150
  propCol<-structure(tail(head(rev(colorRampPalette(c(viridis::rocket(30),'white'),space='Lab')(nCut+nCutTop+nCutBottom)),nCut+nCutTop),nCut),.Names=levels(cuts))
  cols<-structure(propCol[cuts],.Names=colnames(props1))[penn$names]
  labels<-rep('',length(penn$prettyName))
  select<-penn$prettyName %in% rownames(counts)
  
  # Determine the new names to be used as the labels.
  # lab0 <- sprintf('%s\n%d/%d',penn$prettyName[select],counts[penn$prettyName[select],'TRUE'],rowSums(counts[penn$prettyName[select],]))
  lab1 <- penn$prettyName[select]
  # Make placeholder labels that will be reassigned with their new values.
  lab2 <- lab1
  lab3 <- lab1
  # Repopulate those labels with the information of number of positives and total number of counts.
  for(ii in 1:length(lab1)){
    tt <- paste0('pennsylvania,',tolower(lab1[ii]))
    temp1 <- subset(geo3,geo3$name == tt)
    if(nrow(temp1) > 0 ){
      lab2[ii] <- temp1$pos[1]
      lab3[ii] <- temp1$total[1]
      # Add this to make it not grey since grey indicates that there is no data
    }else{
      lab2[ii] <- 0
      lab3[ii] <- 0
    }
  }
  lab2 <- as.numeric(lab2)
  lab3 <- as.numeric(lab3)
  lab4 <- sprintf('%s\n%d/%d',lab1,lab2,lab3)
  
  labels[select]<-lab4
  
  pdf(file_plot_geo_year)
  naCol <- '#EAEAEA'
  maps::map("county", "Pennsylvania", col=ifelse(is.na(cols), naCol, cols), fill=TRUE)
  maps::map.text("county", "Pennsylvania", add=TRUE, label=labels, cex=0.5)
  insetScale(breaks + rep(c(0, .0001), c(length(breaks) - 1, 1)), propCol, main='Proportion positive', insetPos = c(0.025, 0.1, 0.04, 0.3))
  mtext(plot_title, side = 3, line = 2)
  dev.off()
  
}

# Saves a pdf of the Bayesian plot -- county positivity map
get_geo_all <- function(cols,labels,breaks,propCol,file_geo_all){
  pdf(file_geo_all,width = 5.5, height = 5.5)
  # Plot
  naCol <- '#EAEAEA'
  maps::map("county", "Pennsylvania", col=ifelse(is.na(cols), naCol, cols), fill=TRUE)
  maps::map.text("county", "Pennsylvania", add=TRUE, label=labels, cex=0.5)
  insetScale(breaks + rep(c(0, .0001), c(length(breaks) - 1, 1)), propCol, main='Estimated proportion positive', insetPos = c(0.025, 0.1, 0.04, 0.3))
  dev.off()
}

# Function that helps determine which county gets which variant and for which date. -- county positivity map
get_county_date <- function(met4){
  cp1 <- subset(met4,met4$seq_variant != '')
  cp1$county_date <- paste0(cp1$CountyValue,'_',as.character(format(cp1$IncidentDate, "%Y-%m")))
  cp1$county_variant_date <- paste0(cp1$CountyValue,'_',cp1$seq_variant,'_',as.character(format(cp1$IncidentDate, "%Y-%m")))
  cp2 <- data.frame(table(cp1$county_variant_date))
  
  cp2
}
