source("ad_table_ompg_functions.R")
###########################################HB3######################################################
####################################################################################################
mh2d_snps<-read.table("r_ready_snps.vcf", sep="")
mh2d_snps_header<-read.table("snps_header.csv", sep="\t")
mh2d_snps_assembled_ad_table=ad_table(mh2d_snps,3,mh2d_snps_header)
mh2d_snps_assembled_coverage_table=coverage_table(mh2d_snps,3,mh2d_snps_header)
mh2d_snps_alt_ref<-read.table("pass_snps_alt_ref.csv", sep="")
colnames(mh2d_snps_alt_ref)<-c("CHROM","POS","REF","ALT")

mh2d_indels<-read.table("r_ready_indels.vcf", sep="")
mh2d_indels_header<-read.table("indels_header.csv", sep="\t")
mh2d_indels_assembled_ad_table=ad_table(mh2d_indels,3,mh2d_indels_header)
mh2d_indels_assembled_coverage_table=coverage_table(mh2d_indels,3,mh2d_indels_header)
mh2d_indels_alt_ref<-read.table("pass_indels_alt_ref.csv", sep="")

#####################################################################################################

###mh2d_p to mh2d_1a:1h##
mh2d_p_snps_assembled_ad_and_coverage_table=cbind(mh2d_snps_assembled_ad_table[,c(1:11)],mh2d_snps_assembled_coverage_table[,c(3:11)])
mh2d_p_snps_tmp_filter_table=adding_info_to_tables(mh2d_p_snps_assembled_ad_and_coverage_table,0.2,0.8,3,9,10)
mh2d_p_snps_tmp_filter_table=cbind(mh2d_p_snps_tmp_filter_table,mh2d_snps_alt_ref[3],mh2d_snps_alt_ref[4])
mh2d_p_snps_final_filtered_table=ompg_filtering(mh2d_p_snps_tmp_filter_table)

mh2d_p_indels_assembled_ad_and_coverage_table=cbind(mh2d_indels_assembled_ad_table[,c(1:11)],mh2d_indels_assembled_coverage_table[,c(3:11)])
mh2d_p_indels_tmp_filter_table=adding_info_to_tables(mh2d_p_indels_assembled_ad_and_coverage_table,0.2,0.6,3,9,10)
mh2d_p_indels_tmp_filter_table=cbind(mh2d_p_indels_tmp_filter_table,mh2d_indels_alt_ref[3],mh2d_indels_alt_ref[4])
mh2d_p_indels_final_filtered_table=ompg_filtering(mh2d_p_indels_tmp_filter_table)

write.csv(mh2d_p_indels_final_filtered_table,"mh2d_p_indels_final_filtered_table.csv")
write.csv(mh2d_p_indels_tmp_filter_table,"mh2d_p_indels_tmp_filter_table.csv")
write.csv(mh2d_p_snps_final_filtered_table,"mh2d_p_snps_final_filtered_table.csv")
write.csv(mh2d_p_snps_tmp_filter_table,"mh2d_p_snps_tmp_filter_table.csv")

###mh2d_1b to mh2d_2a;2n##
mh2d_1b_snps_assembled_ad_and_coverage_table=cbind(mh2d_snps_assembled_ad_table[,c(1:2,5,12:25)],mh2d_snps_assembled_coverage_table[,c(5,12:25)])
mh2d_1b_snps_tmp_filter_table=adding_info_to_tables(mh2d_1b_snps_assembled_ad_and_coverage_table,0.2,0.8,3,15,10)
mh2d_1b_snps_tmp_filter_table=cbind(mh2d_1b_snps_tmp_filter_table,mh2d_snps_alt_ref[3],mh2d_snps_alt_ref[4])
mh2d_1b_snps_final_filtered_table=ompg_filtering(mh2d_1b_snps_tmp_filter_table)

mh2d_1b_indels_assembled_ad_and_coverage_table=cbind(mh2d_indels_assembled_ad_table[,c(1:2,5,12:25)],mh2d_indels_assembled_coverage_table[,c(5,12:25)])
mh2d_1b_indels_tmp_filter_table=adding_info_to_tables(mh2d_1b_indels_assembled_ad_and_coverage_table,0.2,0.6,3,15,10)
mh2d_1b_indels_tmp_filter_table=cbind(mh2d_1b_indels_tmp_filter_table,mh2d_indels_alt_ref[3],mh2d_indels_alt_ref[4])
mh2d_1b_indels_final_filtered_table=ompg_filtering(mh2d_1b_indels_tmp_filter_table)

write.csv(mh2d_1b_indels_final_filtered_table,"mh2d_1b_indels_final_filtered_table.csv")
write.csv(mh2d_1b_indels_tmp_filter_table,"mh2d_1b_indels_tmp_filter_table.csv")
write.csv(mh2d_1b_snps_final_filtered_table,"mh2d_1b_snps_final_filtered_table.csv")
write.csv(mh2d_1b_snps_tmp_filter_table,"mh2d_1b_snps_tmp_filter_table.csv")

###mh2d_2a to mh2d_3a;3n##
mh2d_2a_snps_assembled_ad_and_coverage_table=cbind(mh2d_snps_assembled_ad_table[,c(1:2,12,26:39)],mh2d_snps_assembled_coverage_table[,c(12,26:39)])
mh2d_2a_snps_tmp_filter_table=adding_info_to_tables(mh2d_2a_snps_assembled_ad_and_coverage_table,0.2,0.8,3,15,10)
mh2d_2a_snps_tmp_filter_table=cbind(mh2d_2a_snps_tmp_filter_table,mh2d_snps_alt_ref[3],mh2d_snps_alt_ref[4])
mh2d_2a_snps_final_filtered_table=ompg_filtering(mh2d_2a_snps_tmp_filter_table)

mh2d_2a_indels_assembled_ad_and_coverage_table=cbind(mh2d_indels_assembled_ad_table[,c(1:2,12,26:39)],mh2d_indels_assembled_coverage_table[,c(12,26:39)])
mh2d_2a_indels_tmp_filter_table=adding_info_to_tables(mh2d_2a_indels_assembled_ad_and_coverage_table,0.2,0.6,3,15,10)
mh2d_2a_indels_tmp_filter_table=cbind(mh2d_2a_indels_tmp_filter_table,mh2d_indels_alt_ref[3],mh2d_indels_alt_ref[4])
mh2d_2a_indels_final_filtered_table=ompg_filtering(mh2d_2a_indels_tmp_filter_table)

write.csv(mh2d_2a_indels_final_filtered_table,"mh2d_2a_indels_final_filtered_table.csv")
write.csv(mh2d_2a_indels_tmp_filter_table,"mh2d_2a_indels_tmp_filter_table.csv")
write.csv(mh2d_2a_snps_final_filtered_table,"mh2d_2a_snps_final_filtered_table.csv")
write.csv(mh2d_2a_snps_tmp_filter_table,"mh2d_2a_snps_tmp_filter_table.csv")

###mh2d_3c to mh2d_4a;4j##
mh2d_3c_snps_assembled_ad_and_coverage_table=cbind(mh2d_snps_assembled_ad_table[,c(1:2,28,41:49)],mh2d_snps_assembled_coverage_table[,c(28,41:49)])
mh2d_3c_snps_tmp_filter_table=adding_info_to_tables(mh2d_3c_snps_assembled_ad_and_coverage_table,0.2,0.8,3,11,10)
mh2d_3c_snps_tmp_filter_table=cbind(mh2d_3c_snps_tmp_filter_table,mh2d_snps_alt_ref[3],mh2d_snps_alt_ref[4])
mh2d_3c_snps_final_filtered_table=ompg_filtering(mh2d_3c_snps_tmp_filter_table)

mh2d_3c_indels_assembled_ad_and_coverage_table=cbind(mh2d_indels_assembled_ad_table[,c(1:2,28,41:49)],mh2d_indels_assembled_coverage_table[,c(28,41:49)])
mh2d_3c_indels_tmp_filter_table=adding_info_to_tables(mh2d_3c_indels_assembled_ad_and_coverage_table,0.2,0.6,3,11,10)
mh2d_3c_indels_tmp_filter_table=cbind(mh2d_3c_indels_tmp_filter_table,mh2d_indels_alt_ref[3],mh2d_indels_alt_ref[4])
mh2d_3c_indels_final_filtered_table=ompg_filtering(mh2d_3c_indels_tmp_filter_table)

write.csv(mh2d_3c_indels_final_filtered_table,"mh2d_3c_indels_final_filtered_table.csv")
write.csv(mh2d_3c_indels_tmp_filter_table,"mh2d_3c_indels_tmp_filter_table.csv")
write.csv(mh2d_3c_snps_final_filtered_table,"mh2d_3c_snps_final_filtered_table.csv")
write.csv(mh2d_3c_snps_tmp_filter_table,"mh2d_3c_snps_tmp_filter_table.csv")

###mh2d_4a to Hb3_5a;5t##
mh2d_4a_snps_assembled_ad_and_coverage_table=cbind(mh2d_snps_assembled_ad_table[,c(1:2,40,50:69)],mh2d_snps_assembled_coverage_table[,c(40,50:69)])
mh2d_4a_snps_tmp_filter_table=adding_info_to_tables(mh2d_3c_snps_assembled_ad_and_coverage_table,0.2,0.8,3,21,10)
mh2d_4a_snps_tmp_filter_table=cbind(mh2d_4a_snps_tmp_filter_table,mh2d_snps_alt_ref[3],mh2d_snps_alt_ref[4])
mh2d_4a_snps_final_filtered_table=ompg_filtering(mh2d_3c_snps_tmp_filter_table)

mh2d_4a_indels_assembled_ad_and_coverage_table=cbind(mh2d_indels_assembled_ad_table[,c(1:2,40,50:69)],mh2d_indels_assembled_coverage_table[,c(40,50:69)])
mh2d_4a_indels_tmp_filter_table=adding_info_to_tables(mh2d_4a_indels_assembled_ad_and_coverage_table,0.2,0.6,3,21,10)
mh2d_4a_indels_tmp_filter_table=cbind(mh2d_4a_indels_tmp_filter_table,mh2d_indels_alt_ref[3],mh2d_indels_alt_ref[4])
mh2d_4a_indels_final_filtered_table=ompg_filtering(mh2d_4a_indels_tmp_filter_table)

write.csv(mh2d_4a_indels_final_filtered_table,"mh2d_4a_indels_final_filtered_table.csv")
write.csv(mh2d_4a_indels_tmp_filter_table,"mh2d_4a_indels_tmp_filter_table.csv")
write.csv(mh2d_4a_snps_final_filtered_table,"mh2d_4a_snps_final_filtered_table.csv")
write.csv(mh2d_4a_snps_tmp_filter_table,"mh2d_4a_snps_tmp_filter_table.csv")

###mh2d_5c to mh2d_6a;6w##
mh2d_5c_snps_assembled_ad_and_coverage_table=cbind(mh2d_snps_assembled_ad_table[,c(1:2,52,70:92)],mh2d_snps_assembled_coverage_table[,c(52,70:92)])
mh2d_5c_snps_tmp_filter_table=adding_info_to_tables(mh2d_5c_snps_assembled_ad_and_coverage_table,0.2,0.8,3,24,10)
mh2d_5c_snps_tmp_filter_table=cbind(mh2d_5c_snps_tmp_filter_table,mh2d_snps_alt_ref[3],mh2d_snps_alt_ref[4])
mh2d_5c_snps_final_filtered_table=ompg_filtering(mh2d_5c_snps_tmp_filter_table)

mh2d_5c_indels_assembled_ad_and_coverage_table=cbind(mh2d_indels_assembled_ad_table[,c(1:2,52,70:92)],mh2d_indels_assembled_coverage_table[,c(52,72:92)])
mh2d_5c_indels_tmp_filter_table=adding_info_to_tables(mh2d_5c_indels_assembled_ad_and_coverage_table,0.2,0.6,3,24,10)
mh2d_5c_indels_tmp_filter_table=cbind(mh2d_5c_indels_tmp_filter_table,mh2d_indels_alt_ref[3],mh2d_indels_alt_ref[4])
mh2d_5c_indels_final_filtered_table=ompg_filtering(mh2d_5c_indels_tmp_filter_table)

write.csv(mh2d_5c_indels_final_filtered_table,"mh2d_5c_indels_final_filtered_table.csv")
write.csv(mh2d_5c_indels_tmp_filter_table,"mh2d_5c_indels_tmp_filter_table.csv")
write.csv(mh2d_5c_snps_final_filtered_table,"mh2d_5c_snps_final_filtered_table.csv")
write.csv(mh2d_5c_snps_tmp_filter_table,"mh2d_5c_snps_tmp_filter_table.csv")














