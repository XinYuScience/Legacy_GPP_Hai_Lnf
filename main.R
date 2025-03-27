#Author: Xin Yu
#Email: xyu@bgc-jena.mpg.de

cluster<-T
if(cluster){
  rootpath<-'/Net/Groups/BGI/people/xyu/1st_MS_code'
}else{
  rootpath<-'X:/1st_MS_code'
}

#source functions code
source(paste0(rootpath,'/Functions_BG.R'))
#load packages
library(ggplot2)
library(dplyr)
library(lubridate)
library(stringi)
library(Kendall)
library(zoo)
library(RColorBrewer)


#path and folders
Datapath<-paste0(rootpath,'/data')
  
folders<-list.dirs(path = Datapath,full.names = TRUE, recursive = TRUE)

for (span in c(7)) {
  for (sd in c(1)) {
    #loop for different sites
    for (p in 2:3) {
      
      #read data  
      target_folder<-folders[p]
      
      DD_file<-list.files(full.names = FALSE, path = target_folder,pattern = 'FULLSET_DD') # daily data's filename
      
      site<-strsplit(DD_file,split = '_',fixed=TRUE)[[1]][2]
      
      df<-read.csv(glue::glue('{target_folder}/{DD_file}')) # read daily data
      
      #create DateTime, doy, month, year columns and order
      df<-create_doy(data=df)
      
      #select needed variables and rename
      df<-df %>%
        rename(
          GPP=GPP_NT_VUT_REF,
          SW_IN=SW_IN_F_MDS,
          TA=TA_F_MDS,
          VPD=VPD_F_MDS,
          P=P_F,
          LE=LE_CORR,# use LE with energy balance correction
          GPP_QC=NEE_VUT_REF_QC,
          SW_IN_QC=SW_IN_F_MDS_QC,
          TA_QC=TA_F_MDS_QC,
          VPD_QC=VPD_F_MDS_QC,
          LE_QC=LE_F_MDS_QC
        )
      
      # cumulative water deficit and water availability calculation
      df$LEperday<-df$LE*3600*24 # use LE with energy balance correction
      df$WdefCum<-computeWdefCum(df$LEperday,df$P) # cumulative water deficit
      awc<-abs(min(df$WdefCum,na.rm = T)) # the bucket size
      WAI<-simpleWAI_ET(precip = df$P, LE=df$LEperday,awc = awc) # water avaibility index
      df$WAI<-WAI$WAI
      

      #interested variables
      vars<-c('GPP','LE','SW_IN','TA','VPD')
      vars_rf_GPP<-c('GPP','SW_IN','TA','VPD')
      vars_qc<-c('GPP_QC','LE_QC','SW_IN_QC','TA_QC','VPD_QC')
      vars_time<-c('Date','year','month','week','doy','TIMESTAMP')
      
      #select variables in df
      df<-df %>% select(vars_time,vars,vars_qc,'SW_IN_POT','WAI')
      
      #quality control
      df_QC<-quality_control(data = df,site = site,vars = c(vars))
      
      if(site=='DE-Lnf'){
        df_QC<-df_QC %>% filter(!(year %in% c(2007,2008,2009)))
      }
      
      #growing season filter
      df_QC_GS<-gs_filter(data = df_QC,vars = c(vars,'WAI'))
      
      #calculate anomalies of intrested variables
      df_QC_GS<-anomaly_all_smooth(data = df_QC_GS,var=c(vars,'WAI')
                            ,site = site,figures_folder = '',span = span,time_scale='daily',sd=sd)
      
      write.csv(df_QC_GS,file = paste0(rootpath,'/results/',site,'_df_QC_GS_WAI_,',span,'_',sd,'.csv'))
           
      #define drought years  
      years<-define_drought(data=df_QC_GS,site=site)# drought year
      
      for (set in c('WAI')) {
        vars_set <- switch(  
          set,
          'WAI'='WAI'
        )  
        #calculate residuals
        diff_normal_GPP<-NA
        diff_legacy_GPP<-NA
        diff_drought_GPP<-NA
        var_imp_GPP<-NA
        var_ex_GPP<-rep(NA,length(years$normal_year))
        for (i in 1:length(c(years$normal_year,years$drought_year))) {
          random_normal_year<-c(years$normal_year,years$drought_year)[i]
          diff_GPP<-model_uncertainty(data=df_QC_GS,random_normal_year = random_normal_year,
                                               figures_folder = figures_folder,site=site,vars = c(vars_rf_GPP,vars_set),EVI_flag = 0)
          var_ex_GPP[i]<-diff_GPP$var_explained
          if (i==1){
            diff_normal_GPP<-diff_GPP$normal
            diff_legacy_GPP<-diff_GPP$legacy
            diff_drought_GPP<-diff_GPP$drought
            var_imp_GPP<-diff_GPP$var_importance
          }else{
            diff_normal_GPP<-rbind(diff_normal_GPP,diff_GPP$normal)
            diff_legacy_GPP<-rbind(diff_legacy_GPP,diff_GPP$legacy)
            diff_drought_GPP<-rbind(diff_drought_GPP,diff_GPP$drought)
            var_imp_GPP<-rbind(var_imp_GPP,diff_GPP$var_importance)
          }
        }
        
        var_ex<-data.frame(normal_year=c(years$normal_year,years$drought_year),
                           var_ex_GPP=var_ex_GPP)
        write.csv(var_ex,file = paste0(rootpath,'/results/',site,'_var_ex_abs_',set,'_',span,'_',sd,'.csv'))
        write.csv(var_imp_GPP,file = paste0(rootpath,'/results/',site,'_var_imp_GPP_abs_',set,'_',span,'_',sd,'.csv'))
        
        diff_normal_GPP<-create_doy(diff_normal_GPP)
        diff_legacy_GPP<-create_doy(diff_legacy_GPP)
        
        write.csv(diff_normal_GPP,file = paste0(rootpath,'/results/',site,'_diff_normal_abs_',set,'_',span,'_',sd,'.csv'))
        write.csv(diff_legacy_GPP,file = paste0(rootpath,'/results/',site,'_diff_legacy_abs_',set,'_',span,'_',sd,'.csv'))
      }
      
    }
  }
}