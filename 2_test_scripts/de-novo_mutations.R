## Author: Aakanksha Singh
## Date:January 2021
## Version: 1.0
## Description: The script contains functions for the identification of de-novo mutations after the GATK variant calling using the Allele Depth (AD) values


### Part 1 - Creating Functions ####
#######################################


library(dplyr)

ad_table <- function(input_table.df,start_sample_column,header_table)
# This functons builds ad table
# Args: input_table.df : input table with AD values
#       start_sample_column: columns where the Ad values start
#       header_table : The table containing the header to input_table.df
# Returns: AD table output
{
    testing_ad.df<- input_table.df
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
    if (nrow(new_table.df)==nrow(input_table.df))
        {
            new_table.df=cbind(input_table.df[1:start_sample_column-1],new_table.df)
        }
    new_table.df <- new_table.df[ -c(3) ] ### For some reason there is a empty column at 3rd positin so i remove it, not sure why or how it gets there

    #####Applying header to the newly built table
    data_char2 <- header_table                                             # Duplicate data
    fac_cols <- sapply(data_char2, is.factor)                              # Identify all factor columns
    data_char2[fac_cols] <- lapply(data_char2[fac_cols], as.character)     # Convert all factors to characters
    ## I insert a loop to add coverage after each each sample sample so that it works when filtering
    ##(If i do not do that then, dply filter wont work in ompg filtering due to duplicate column names in the ad and coverage table)
    for(k in 3:ncol(data_char2))
        {
            titlename=as.character(data_char2[1,k])
            new_titlename=paste(titlename,"alt_prop",sep="_")
            data_char2[1,k]=new_titlename
         }
    colnames(new_table.df)<-(data_char2[1,]) #put it as header of the table
    return(new_table.df)
}



coverage_table <- function(input_table.df,start_sample_column,header_table)
# This functons builds coverage table
# Args: input_table.df : input table with AD values
#       start_sample_column: columns where the Ad values start
#       header_table : The table containing the header to input_table.df
# Returns: coverage table output
{
    testing_ad.df<- input_table.df
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
        {
            new_table.df=cbind(input_table.df[1:start_sample_column-1],new_table.df)
        }
    new_table.df <- new_table.df[ -c(3) ] ### For some reason there is a empty column at 3rd positin so i remove it, not sure why or how it gets there
    #####Applying header to the newly built table
    data_char2 <- header_table                                             # Duplicate data
    fac_cols <- sapply(data_char2, is.factor)                           # Identify all factor columns
    data_char2[fac_cols] <- lapply(data_char2[fac_cols], as.character)   # Convert all factors to characters
    ## I insert a loop to add coverage after each each sample sample so that it works when filtering
    ##(If i do not do that then, dply filter wont work in ompg filtering due to duplicate column names in the ad and coverage table)
    for(k in 3:ncol(data_char2))
        {
            titlename=as.character(data_char2[1,k])
            new_titlename=paste(titlename,"coverage",sep="_")
            data_char2[1,k]=new_titlename
        }
    colnames(new_table.df)<-(data_char2[1,]) #put it as header of the table
    return(new_table.df)
}

combined_ad_and_coverage <- function(input_table.df,start_sample_column,header_table)
#This functons combines the ad and the coverage table
# Args: input_table.df : input table with AD values
#       start_sample_column: columns where the Ad values start
#       header_table : The table containing the header to input_table.df
# Returns: combined Ad and coverage table
{
    ad_ratio_table=ad_table(input_table.df,start_sample_column,header_table)
    coverage_table_from_gatk_values=coverage_table(input_table.df,start_sample_column,header_table)
    final_ad_and_coverage_table=cbind(ad_ratio_table[,c(1:ncol(ad_ratio_table))],coverage_table_from_gatk_values[,c(3:ncol(coverage_table_from_gatk_values))])
    return(final_ad_and_coverage_table)
}


adding_info_to_tables <- function(new_table.df,alt_ratio_min,alt_ratio_max,start_sample_column,total_no_of_samples,coverage_limit)
#This functons adds additional info for identifying de-novo mutations
    # Args: new_table.df : input table with AD values
    #       alt_ratio_min: minimum alt_ratio for 'Reference'
    #       alt_ratio_max: max alt_ratio for 'Alternate'
    #       start_sample_column: columns where the Ad values start
    #       total_no_of_samples: total no of samples analysed
    #       coverage_limit: The table containing the header to input_table.df
    # Returns: A dataframe with 7 additional columns
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
    for (i in 1:nrow(testing.df))
        {
            alt=0
            ref=0
            mixed=0
            no_avail=0
            pass_coverage=0
            samp_name=NA
            alt_prop_of_pass_sample=NA
            for (j in start_sample_column:(start_sample_column+total_no_of_samples-1))
                { # (-1) because the start sample column already inculdes the first column
                    if (isTRUE(testing.df[i,j] <= alt_ratio_min))
                        {
                            ref=ref+1
                        }
                    else if (isTRUE(testing.df[i,j] > alt_ratio_min & testing.df[i,j] < alt_ratio_max))
                        {
                            mixed=mixed+1
                        }
                    else if (isTRUE(testing.df[i,j] >= alt_ratio_max))
                        {
                            alt=alt+1
                            if (isTRUE(alt == 1))
                            {
                                if (isTRUE(testing.df[i,(j+total_no_of_samples)] >= coverage_limit))
                                    {
                                        pass_coverage = as.numeric(testing.df[i,(j+total_no_of_samples)])
                                        samp_name=colnames(testing.df[j])
                                        alt_prop_of_pass_sample=as.numeric(testing.df[i,j])
                                        
                                    }
                                    
                            } ###The value of pass_coverage will be one only and only when the coverage is greater than the given value and only one one value of alt
                    else
                    {
                        pass_coverage = 0
                        samp_name=NA
                        alt_prop_of_pass_sample=NA}
                    }
                    else
                    {
                        no_avail=no_avail+1
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
  #adding the columns to the table
  testing.df <- cbind(testing.df,Reference)
  testing.df <- cbind(testing.df,Alternate)
  testing.df <- cbind(testing.df,Mixed)
  testing.df <- cbind(testing.df,No_info)
  testing.df <- cbind(testing.df,Check_coverage)
  testing.df <- cbind(testing.df,Sample_name)
  testing.df <- cbind(testing.df,Alt_proportation_of_pass_sample)
  return(testing.df)
}



ompg_filtering <- function(testing.df)
#This functons filters de-novo mutations from  additional info table
# Args: testing.df : output from additional_info table
    {
        filtered_table <- testing.df %>% filter(testing.df$Alternate==1,testing.df$Check_coverage>0)
        return(filtered_table)
    }

#######  PART 2 - template for calling the Functions #####
##########################################################

snps<-read.table("snps.txt", sep="")
snps_header<-read.table("header.txt", sep="\t")
combined_ad_and_coverage=combined_ad_and_coverage(snps,3,snps_header)
snps_tmp_filter_table=adding_info_to_tables(combined_ad_and_coverage,0.2,0.8,3,3,8)
snps_final_filtered_table=ompg_filtering(snps_tmp_filter_table)
View(snps_final_filtered_table)







