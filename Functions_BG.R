#Cumulative water deficit
computeWdefCum <- function(LE, precip) {
  
  n <- length(LE)
  ET <- LE / 2.45E+6    # J m-2 ==> kg m-2 (==mm)
  wdefCum <- rep(NA, n)
  wdefCum[1]  <- 0
  
  for (i in 2:n) {
    wdefCum[i] <- min(wdefCum[i-1]+precip[i]-ET[i],0)
    if (is.na(wdefCum[i])) {wdefCum[i] <- wdefCum[i-1] }
  } 
  
  wdefCum
  
}

#Water availability index
simpleWAI_ET  <- function(precip=NULL, LE=NULL, awc=awc) {
  
  n <- length(precip)
  ET <- LE / 2.45E+6    # J m-2 ==> kg m-2 (==mm)
  m<-365*5
  WAI_0 <- rep(NA, m)
  WAI_0[1]  <- awc*0.5
  P_0<-rep(precip[1:365],5)
  ET_0<-rep(ET[1:365],5)
  
  input_0  <- rep(NA,m)
  
  for (i in 2:m) {
    input_0[i] <- min(P_0[i], awc - WAI_0[i-1])
    WAI_0[i] <- max(WAI_0[i-1] + input_0[i] - ET_0[i],0)
    
    if (is.na(WAI_0[i])) {WAI_0[i] <- WAI_0[i-1] }
    
  } 
  
  WAI <- rep(NA, n)
  WAI[1]  <- WAI_0[length(WAI_0)]

  
  input  <- rep(NA, n)
  
  for (i in 2:n) {
    input[i] <- min(precip[i], awc - WAI[i-1])
    WAI[i] <- max(WAI[i-1] + input[i] - ET[i],0)
    if (is.na(WAI[i])) {WAI[i] <- WAI[i-1] }
  } 
  
  data.frame(WAI=WAI, WAIinput=input)
}

#quantify legacy effects with uncertainty
model_uncertainty<-function(data=data,random_normal_year=normal_year,
                                 figures_folder=figures_folder,site=site,
                                 vars=vars,EVI_flag=EVI_flag){
  
  #select training data
  training_data<-data %>% filter(!(year %in% c(years$legacy_year,random_normal_year)))
  legacy_data<-data %>% filter(year %in% years$legacy_year)
  random_normal_data<-data %>% filter(year %in% random_normal_year)
  
  if(EVI_flag==1){
    vars<-c(vars,'EVI')
  }
  
  #Random Forest
  library(randomForest)
  # for reproduciblity
  set.seed(123)
  vars<-paste(vars,'_Anom_smooth',sep='')
  train<-training_data %>% select(vars,'SW_IN_POT')
  rf<-randomForest(formula = GPP_Anom_smooth ~ .,data = train,ntree=400,na.action = na.exclude,
                   importance=T,mtry=5,nodesize=5)
  
  var_importance<-as.data.frame(importance(rf))
  var_importance$normal_year<-random_normal_year
  var_explained<-rf$rsq[length(rf$rsq)]
  
  #normalization and residuals calculation
  random_normal_data$GPP_Anom_smooth_rf<-predict(rf,newdata = random_normal_data)
  random_normal_data$GPP_Anom_smooth_rf_diff<-random_normal_data$GPP_Anom_smooth-random_normal_data$GPP_Anom_smooth_rf
  
  legacy_data$GPP_Anom_smooth_rf<-predict(rf,newdata = legacy_data)
  legacy_data$GPP_Anom_smooth_rf_diff<-legacy_data$GPP_Anom_smooth-legacy_data$GPP_Anom_smooth_rf
  legacy_data$normal_year<-rep(random_normal_year,length(legacy_data$Date))
  
 
  #return diff
  normal<-random_normal_data %>% select(Date,year,doy,GPP_Anom_smooth,GPP_Anom_smooth_rf,GPP_Anom_smooth_rf_diff)
  legacy<-legacy_data %>% select(Date,year,doy,GPP_Anom_smooth,GPP_Anom_smooth_rf,GPP_Anom_smooth_rf_diff,normal_year)
  
  
  diff<-list(normal=normal,legacy=legacy,var_importance=var_importance,var_explained=var_explained)
  return(diff)
}


#growing season filter
gs_filter<-function(data=data,vars=vars,quan_GPP=0.25){
  #define growing season
  mean_GPP<-data %>% 
    select(doy,GPP) %>%
    dplyr::group_by(doy) %>%
    summarize(seasonal = mean(GPP,na.rm = T))
  #calculate seasonal cycle
  
  span<-7
  fit <- with(mean_GPP,
              ksmooth(doy,seasonal, kernel = "normal", bandwidth = span))
  mean_GPP$seasonal_smooth<-fit$y
  mean_GPP_GS<-mean_GPP %>% filter(mean_GPP$seasonal_smooth>=quan_GPP*max(mean_GPP$seasonal_smooth,na.rm = T))
  doy_GS<-mean_GPP_GS$doy # doy in growing season
  
  
  # data filter of growing season
  for (i in 1:length(vars)){
    if (vars[i] %in% colnames(data)){
      data[[vars[i]]][data$doy<doy_GS[1] | data$doy>doy_GS[length(doy_GS)]]<-NA
    }
  }
  return(data)
}

#create DateTime, doy, month, year columns and order
create_doy<-function(data=data){
  if('TIMESTAMP' %in% colnames(data)){
    data$Date<- ymd(data$TIMESTAMP)
  }
  data$Date<-as.Date(data$Date)
  #browser()
  data<- mutate(data,
    doy=yday(Date),
    week=week(Date),
    month=month(Date),
    year=year(Date)
  ) 
  
  data<-select(data,Date,year,month,week,doy,everything())
  return(data)
}

#calculate anomalies of all variables
anomaly_all_smooth<-function(data=data,vars=c('SWC_1','SWC_2','TA','SW_IN','TS_1','TS_2','VPD','GPP_NT')
                      ,site=site,figures_folder=figures_folder,span=span,time_scale=time_scale,sd=sd){
  for (i in 1:length(vars)){
    if (vars[i] %in% colnames(data)){
      if (time_scale=='daily'){
        data<-anomaly_smooth(data=data,var = vars[i],vars=vars,site = site,figures_folder = figures_folder,span=span,sd=sd)
      }else if(time_scale=='weekly'){
        data<-anomaly_weekly(data=data,var = vars[i],vars=vars,site = site,figures_folder = figures_folder,sd=sd)
      }
    }
  }
  return(data)
}

#calculate anomaly
anomaly_smooth<-function(data = df,var='SWC_1',vars=vars,site = site,figures_folder=figures_folder,span=span,sd=sd){
 
  data_QC<-data
  
  #remove long-term trend
  #mk test
  mk<-MannKendall(data_QC[[var]])
  if (mk$sl[1]<=0.05){
    linearMod <- lm(get(var) ~ as.Date(Date), data=data_QC,na.action=na.exclude)
    data_QC[[paste(var,'_detrend', sep="")]]<-residuals(linearMod)
  }else{
    data_QC[[paste(var,'_detrend', sep="")]]<-data_QC[[var]]
  }
  
  #calculate seasonal cycle
  #mean
  mean_data<-data_QC %>% group_by(doy) %>%
    summarize(seasonal = mean(get(paste(var,'_detrend', sep="")),na.rm = T))
  
  
  #add mean seasonal cycle to existing data using matching value 'doy'
  data_QC[[paste(var,'_detrend_seasonal', sep="")]]<- mean_data$seasonal[match(data_QC$doy,mean_data$doy)]
  
  #calculate anomaly
  data_QC[[paste(var,'_Anom', sep="")]]<-data_QC[[paste(var,'_detrend', sep="")]]-data_QC[[paste(var,'_detrend_seasonal', sep="")]]
  
  data_QC[[paste(var,'_Anom_smooth', sep="")]]<-rollapply(data_QC[[paste(var,'_Anom', sep="")]]
                                                          , width = span, by = 1, FUN = mean, na.rm = T,fill=NA)
  
  return(data_QC)
}

#quality control
quality_control<-function(data=data,vars=vars,site=site,flag=0.7){
  #browser()
  for (i in 1:length(vars)){
    if (vars[i] %in% colnames(data)){
        qc_var<-paste(vars[i],'_QC',sep = '')
        data[[qc_var]][data[[qc_var]]==-9999]<-1
        data[[vars[i]]][data[[qc_var]]<flag]<-NA
        data[[vars[i]]][data[[vars[i]]]==-9999]<-NA
        if(vars[i] == 'GPP'){
          data[[vars[i]]][data[[vars[i]]]<0]<-NA
        }
    }
  }
  return(data)
}

#define drought years
define_drought<-function(data=data,legacy_length=2,site=site){
  if(site=='DE-Hai'){
    drought_year<-c(2003,2018)
  }else{
	  drought_year<-c(2003)
  }
  if(legacy_length==2){
    drought_legacy_year<-sort(unique(c(drought_year,drought_year+1,drought_year+2)))# drought and legacy year
  }else if(legacy_length==3){
    drought_legacy_year<-sort(unique(c(drought_year,drought_year+1,drought_year+2,drought_year+3)))# drought and legacy year
  }
  drought_legacy_year<-drought_legacy_year[drought_legacy_year<=max(data$year)]# remove years beyond the period
  legacy_year<-sort(drought_legacy_year[is.na(pmatch(drought_legacy_year,drought_year))])# select pure legacy year
  normal_year<-unique(data$year)[is.na(pmatch(unique(data$year),drought_legacy_year))]
  normal_drought_year<-c(normal_year,drought_year)
  
  years<-list(drought_year=drought_year,legacy_year=legacy_year,drought_legacy_year=drought_legacy_year
              ,normal_year=normal_year, normal_drought_year= normal_drought_year)
  return(years)
}
