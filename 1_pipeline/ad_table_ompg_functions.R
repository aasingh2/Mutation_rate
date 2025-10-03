
##############Updated and modified successfully on 29th january 2021####################################

######################################FUNCTIONS#######################################
######################################################################################
#################### To filter out de novo mutations##################################
#Inputs required for the function
#1.new_table.df=the table with the ad values for just a single generation where data starts from row one and a valid header
#2.alt_ratio_min= cut off alt_ratio to say that a sample is ref
#3.alt_ratio_max=cut off alt_ratio to say that a sample is alt
#4.start_sample_column= column from which the samples start
# alt_ratio is defined as no of alternate reads/(no of alternate reads+no of reference reads)
library(dplyr)

#obsolute function, its new modified version is adding_info_to_tables which also takes into account the coverage
ompg_ad <- function(new_table.df,alt_ratio_min,alt_ratio_max,start_sample_column) 
{
  testing.df=new_table.df
  ##sample size is total no of columns - the columns from where the sample starts, -1 because start sample column includes the column alraedy containing the sample
  sample_size=ncol(testing.df)-(start_sample_column-1)
  Reference <- c()
  Alternate <- c()
  Mixed <- c()
  No_info<- c()
  for (i in 1:nrow(testing.df)) {
    alt=0
    ref=0
    mixed=0
    no_avail=0
    for (j in start_sample_column:ncol(testing.df)) {
      if (isTRUE(testing.df[i,j] < alt_ratio_min)){
        ref=ref+1}
      else if (isTRUE(testing.df[i,j] > alt_ratio_min & testing.df[i,j]< alt_ratio_max)){
        mixed=mixed+1}
      else if (isTRUE(testing.df[i,j] > alt_ratio_max)){
        alt=alt+1} else 
        {no_avail=no_avail+1}
    }
    total=ref+alt+mixed+no_avail
    Reference <- c(Reference,ref/total)
    Alternate <- c(Alternate,alt/total)
    Mixed <- c(Mixed,mixed/total)
    No_info <- c(No_info,no_avail/total)
  }
  ##################adding the columns to the table ##########################################
  testing.df <- cbind(testing.df,Reference)
  testing.df <- cbind(testing.df,Alternate)
  testing.df <- cbind(testing.df,Mixed)
  testing.df <- cbind(testing.df,No_info)
  
  ########################Filtering rows and saving the csv file ##############################
  ###Only keeping the rows where alternate is less than or equal to 0.10
  ##Adjusting the ompg frequency to +- 5 percent
  ompg=1/sample_size
  ompg_adjust=ompg+(0.05*ompg)
  filtered_table <- testing.df %>% filter(testing.df$Alternate>0,testing.df$Alternate<=ompg_adjust,testing.df$Reference>0.60)
  return(filtered_table)
}

###########################################################################################################
##########################Function to build Ad ratio table from the gatk table#############################
#Requirements
#1.input_table.df The AD table received from extracting all the AD values should have been pre-processed(prefereable bash)
#to replace the comma with a space so that eache value of AD falls in a different column. From the pipeline I used, the first
#value is ref and the next value that appears is the alt so I calculate the alt/(ref+alt) for each sample.
#IMPORTANT NOTE- WORKS ONLY FOR BIALLELIC VARIANTS
#2.start_sample_column- the column number from which the actual data starts
#3.header_table- before preprocessing your data to replace all the commas with a string extract the header of the original file in a vcf
#preferable separated by a comma or a space, mention this, when you call the header file using the read table command and
# do not put the header=TRUE argument when calling the read table.
ad_table <- function(input_table.df,start_sample_column,header_table) 
{testing_ad.df<- input_table.df
new_table.df=data.frame()
for (i in 1:nrow(testing_ad.df)) 
{
  freq_colm=0
  for(j in seq(from = start_sample_column, to = ncol(testing_ad.df), by = 2))
  {
    freq=round(testing_ad.df[i,j+1]/(testing_ad.df[i,j]+testing_ad.df[i,j+1]), digits = 2)
    freq_colm=c(freq_colm,freq)
  }
  new_table.df=rbind(new_table.df,freq_colm)  
}

if (nrow(new_table.df)==nrow(input_table.df)) # & ncol(input_table.df) == ((ncol(new_table.df))/2)+2 )
{new_table.df=cbind(input_table.df[1:start_sample_column-1],new_table.df)}
new_table.df <- new_table.df[ -c(3) ] ### For some reason there is a empty column at 3rd positin so i remove it, not sure why or how it gets there
#####Applying header to the newly built table
data_char2 <- header_table                                             # Duplicate data
fac_cols <- sapply(data_char2, is.factor)                           # Identify all factor columns
data_char2[fac_cols] <- lapply(data_char2[fac_cols], as.character)   # Convert all factors to characters
  ## I insert a loop to add coverage after each each sample sample so that it works when filtering 
  ##(If i do not do that then, dply filter wont work in ompg filtering due to duplicate column names in the ad and coverage table)
for(k in 3:ncol(data_char2))
{titlename=as.character(data_char2[1,k])
new_titlename=paste(titlename,"alt_prop",sep="_")
data_char2[1,k]=new_titlename }
colnames(new_table.df)<-(data_char2[1,]) #put it as header of the table
return(new_table.df)
}

###########################################################################################################
##########################Function to build coverage table from the gatk table#############################
#Requirements
#1.input_table.df The AD table received from extracting all the AD values should have been pre-processed(prefereable bash)
#to replace the comma with a space so that eache value of AD falls in a different column. From the pipeline I used, the first
#value is ref and the next value that appears is the alt so I calculate the alt/(ref+alt) for each sample.
#IMPORTANT NOTE- WORKS ONLY FOR BIALLELIC VARIANTS
#2.start_sample_column- the column number from which the actual data starts
#3.header_table- before preprocessing your data to replace all the commas with a string extract the header of the original file in a vcf
#preferable separated by a comma or a space, mention this, when you call the header file using the read table command and
# do not put the header=TRUE argument when calling the read table.
coverage_table <- function(input_table.df,start_sample_column,header_table) 
{testing_ad.df<- input_table.df
new_table.df=data.frame()
for (i in 1:nrow(testing_ad.df)) 
{
  freq_colm=0
  for(j in seq(from = start_sample_column, to = ncol(testing_ad.df), by = 2))
  {
    freq=testing_ad.df[i,j]+testing_ad.df[i,j+1]
    freq_colm=c(freq_colm,freq)
  }
  new_table.df=rbind(new_table.df,freq_colm)  
}

if (nrow(new_table.df)==nrow(input_table.df)) # & ncol(input_table.df) == ((ncol(new_table.df))/2)+2 )
{new_table.df=cbind(input_table.df[1:start_sample_column-1],new_table.df)}
new_table.df <- new_table.df[ -c(3) ] ### For some reason there is a empty column at 3rd positin so i remove it, not sure why or how it gets there
#####Applying header to the newly built table
data_char2 <- header_table                                             # Duplicate data
fac_cols <- sapply(data_char2, is.factor)                           # Identify all factor columns
data_char2[fac_cols] <- lapply(data_char2[fac_cols], as.character)   # Convert all factors to characters
         ## I insert a loop to add coverage after each each sample sample so that it works when filtering 
         ##(If i do not do that then, dply filter wont work in ompg filtering due to duplicate column names in the ad and coverage table)
for(k in 3:ncol(data_char2))
{titlename=as.character(data_char2[1,k])
new_titlename=paste(titlename,"coverage",sep="_")
data_char2[1,k]=new_titlename }
colnames(new_table.df)<-(data_char2[1,]) #put it as header of the table
return(new_table.df)
}

###########################################FUNCTIONS END#####################################################
#############################################################################################################
#############################################################################################################

#########Example of how to call the function#########
######ad_table function#######
#Reading the raw Ad table
#ad_vcf <- read.table("/Users/aakankshasingh/Desktop/tmp_folder/3d7_gatk_indels/ready-indels.vcf", sep="", quote="")
#reading the header file
#header.df<-read.table("/Users/aakankshasingh/Desktop/tmp_folder/3d7_gatk_indels/header.txt",sep=",")
#Producing the AD table which calculate the alt/alt+ref ratio
#ADtable=ad_table(ad_vcf,3,header.df)
########ompg_ad function#########
#here i use the table generated by the previous table as an input
#ompg_in_sample_table<-ompg_ad(ADtable,0.2,0.3,3)

####Modifying the ompg table#######(Works as planned)
#What this function does
#1. Calculate the number of the samples that are that a reference, alternate or mixed status at a particular loci according to the cutoff range inputted by the user
#2. Check1- if only one sample has the alternate status, then it checks its coverage and outputs 1 if the coverage is equal to or greater than the cutoff value inputtted by 
#the user
#3. Check 2- if the value of the alteranate column and the coverage column is one at the particular loci then it gives the sample name and its alt ratio at that position
#4. Input required for the function
##  new_table.df=the table with the ad values for just a single generation where data starts from row one and a valid header
##  alt_ratio_min= cut off alt_ratio to say that a sample is ref
##  alt_ratio_max=cut off alt_ratio to say that a sample is alt
##  start_sample_column= column from which the samples start
##  total_no_of_samples = total number of the samples that you have put in the input table 
##  coverage_limit- The minimum cut off value of covearge to qualify as a true positive
#5.Additional Info
##  alt_ratio is defined as no of alternate reads/(no of alternate reads+no of reference reads)
##  total number of samples-( mind you it is not the total no of columns in most of the cases is ncol(table)-2 ie excluding the CHROM and POS columns)
##        
adding_info_to_tables <- function(new_table.df,alt_ratio_min,alt_ratio_max,start_sample_column,total_no_of_samples,coverage_limit) 
{
  testing.df=new_table.df
  sample_size=total_no_of_samples
  Reference <- c()
  Alternate <- c()
  Mixed <- c()
  No_info<- c()
  Check_coverage <-c()
  Sample_name <-c()
  Alt_proportation_of_pass_sample <- c()
  for (i in 1:nrow(testing.df)) {
    alt=0
    ref=0
    mixed=0
    no_avail=0
    pass_coverage=0
    samp_name=NA
    alt_prop_of_pass_sample=NA
    for (j in start_sample_column:(start_sample_column+total_no_of_samples-1)) { #-1 because the start sample column already inculdes the first column
      if (isTRUE(testing.df[i,j] <= alt_ratio_min)){
        ref=ref+1
        }
      else if (isTRUE(testing.df[i,j] > alt_ratio_min & testing.df[i,j] < alt_ratio_max)){
        mixed=mixed+1
        }
      else if (isTRUE(testing.df[i,j] >= alt_ratio_max)){
        alt=alt+1
        if (isTRUE(alt == 1))
        {if (isTRUE(testing.df[i,(j+total_no_of_samples)] >= coverage_limit))
          {pass_coverage = as.numeric(testing.df[i,(j+total_no_of_samples)])
           samp_name=colnames(testing.df[j]) 
           alt_prop_of_pass_sample=as.numeric(testing.df[i,j])}} ###The value of pass_coverage will be one only and only when the coverage is greater than the given value and only one one value of alt
         else  
           {pass_coverage = 0
           samp_name=NA
           alt_prop_of_pass_sample=NA}
        }
      else 
      {no_avail=no_avail+1
      }
    }
    Check_coverage <-c(Check_coverage,pass_coverage)
    Reference <- c(Reference,ref)
    Alternate <- c(Alternate,alt)
    Mixed <- c(Mixed,mixed)
    No_info <- c(No_info,no_avail)
    Sample_name <-c(Sample_name,samp_name)
    Alt_proportation_of_pass_sample <- c(Alt_proportation_of_pass_sample,alt_prop_of_pass_sample)
  }
  ##################adding the columns to the table ##########################################
  testing.df <- cbind(testing.df,Reference)
  testing.df <- cbind(testing.df,Alternate)
  testing.df <- cbind(testing.df,Mixed)
  testing.df <- cbind(testing.df,No_info)
  testing.df <- cbind(testing.df,Check_coverage)
  testing.df <- cbind(testing.df,Sample_name)
  testing.df <- cbind(testing.df,Alt_proportation_of_pass_sample)
  ########################Filtering rows and saving the csv file ##############################
  ###Only keeping the rows where alternate is less than or equal to 0.10
  ##Adjusting the ompg frequency to +- 5 percent
  #filtered_table <- testing.df %>% filter(testing.df$Alternate==1 |testing.df$Reference==1)
  return(testing.df)
}

ompg_filtering <- function(testing.df) 
{ filtered_table <- testing.df %>% filter(testing.df$Alternate==1,testing.df$Check_coverage>0)
  return(filtered_table)
  }

##Common Commands required are
#assembled_ad_table=ad_table(input_table.df,start_sample_column,header_table)
#assembled_coverage_table=coverage_table(input_table.df,start_sample_column,header_table)
#Joining the two tables together
#assembled_ad_and_coverage_table=cbind(assembled_ad_table,assembled_coverage_table[,]
#tmp_filter_table=adding_info_to_tables(assembled_ad_and_coverage_table,alt_ratio_min,alt_ratio_max,start_sample_column,number_of_samples,coverage_limt)
#final_filtered_table=ompg_filtering(tmp_filter_table)

ad_table_modified <- function(input_table.df,start_sample_column,header_table) 
{testing_ad.df<- input_table.df
new_table.df=data.frame()
for (i in start_sample_column:nrow(testing_ad.df)) 
{ new_table.df[i]=testing_ad.df[i]/(testing_ad.df[i]+testing_ad.df[i+1])
  }

if (nrow(new_table.df)==nrow(input_table.df)) # & ncol(input_table.df) == ((ncol(new_table.df))/2)+2 )
{new_table.df=cbind(input_table.df[1:start_sample_column-1],new_table.df)}
new_table.df <- new_table.df[ -c(3) ] ### For some reason there is a empty column at 3rd positin so i remove it, not sure why or how it gets there
#####Applying header to the newly built table
data_char2 <- header_table                                             # Duplicate data
fac_cols <- sapply(data_char2, is.factor)                           # Identify all factor columns
data_char2[fac_cols] <- lapply(data_char2[fac_cols], as.character)   # Convert all factors to characters
## I insert a loop to add coverage after each each sample sample so that it works when filtering 
##(If i do not do that then, dply filter wont work in ompg filtering due to duplicate column names in the ad and coverage table)
for(k in 3:ncol(data_char2))
{titlename=as.character(data_char2[1,k])
new_titlename=paste(titlename,"alt_prop",sep="_")
data_char2[1,k]=new_titlename }
colnames(new_table.df)<-(data_char2[1,]) #put it as header of the table
return(new_table.df)
}










