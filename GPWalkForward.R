## ---------------------------
##
## Script name: Walk Forward
##
## Purpose of script: Evaluate Time Series Models with Walkforward validation
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 06-04-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

################################################################################

### Load Libs

################################################################################

library(readxl)
library(ggplot2)
library(roll)
library(gridExtra)
library(readr)
#library(tidyverse)
library(lattice)
library(readr)

################################################################################

###
###                               Functions
###

################################################################################

# ------------------------------------------------------------------------------

## Source Main Supporting Funcs

#source('./support/functions.R', echo=F)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

# Evaludate Regression Model

eval.model = function(f = 3,col,col_g){
  
  df = data.frame(matrix(ncol = 0, nrow = length(col)))
  res.df = data.frame(matrix(ncol = 0, nrow = f))
  
  for (i in 1:(nrow(df))) {
    #row growth rate
    g = col_g[i]
    #row current cases
    cases = col[i]
    #next pred values
    forecasts = cases * (g ^ (1:f))
    #actual vals
    vals = col[(i+1):(i+f)]
    #residuals
    res = forecasts - vals
    #Save forecasts
    df[(i+1):(i+f),(ncol(df)+1)] = forecasts
    #save residuals
    res.df[,i] = res
  }
  
  #Remove extra Rows
  df = df[1:nrow(boot.summary),]
  # Remnove Columns with NA predictions
  df = Filter(function(x)!all(is.na(x)), df)
  
  res.df = Filter(function(x)!all(is.na(x)), res.df)
  
  # Compute Mean squaured error (we already have the error)
  res.vec = as.vector(as.matrix(res.df))
  res.vec = res.vec[!is.na(res.vec)]
  res.vec = res.vec ^ 2
  MSE = sum(res.vec) / length(res.vec)
  
  return(list("Forecasts" = df,
              "Error" = MSE))
  
}


# ------------------------------------------------------------------------------

# Function to Compute growth and MSE over param ranges

compute_error = function(x,w,f,cases,lower = 7,FUN = mean){
  
  #check if observations to few
  if(nrow(x)<=lower | nrow(x) == nrow(cases)){
    return(NULL)
  }
  
  future_case = cases$cases[-(1:nrow(x))]
  res <- expand.grid(window=w, forecast=f)
  res$growth = rep(0,nrow(res))
  res$mse = rep(0,nrow(res))
  
  for (i in 1:nrow(res)) {
    # Set Loop params for w & f to iterate with
    w_loop = res$window[i]
    f_loop = res$forecast[i]
    
    #compute growth rate for as the max over the past w_loop days
    #Take max of the first obersvation to keep indexs within positive range
    res$growth[i] =  FUN(exp(x$gr[max(1,(length(x$gr)+1-w_loop)):length(x$gr)]))
    
    print(res$growth[i])
    
    # For the last day in GP dataset forecast values for the next f
    l_est = exp(x$`log(OD)`)[nrow(x)]
    forecasts = l_est * (res$growth[i] ^ (1:f_loop))
    
    #Compute Residuals Between Forecasts and actual
    # Take the mimimum length of both vectors
    comp.index = mean(length(forecasts),length(future_case))
    res.vec = forecasts[1:comp.index] - future_case[1:comp.index]
    
    #Compute MSE and MSE with extra bias against neg res
    res$mse[i] = mean(res.vec * res.vec)
    res.vec[which(res.vec < 0)] = res.vec[which(res.vec < 0)] ^ 2
    res$mse2[i] = mean(res.vec * res.vec)
  }
  return(res)
}


# ------------------------------------------------------------------------------

# Reshape Residual df for plotting

exp.res.df = function(win_r,f_range,gp.df){
  res_b <- expand.grid(window=win_r, forecast=f_range)
  res_b$mse = rep(0,nrow(res_b))
  res_b$mse2 = rep(0,nrow(res_b))
  
  for (i in 1:length(gp.df)) {
    res_b$mse = res_b$mse + gp.df[[1]]$mse
    res_b$mse2 = res_b$mse2 + gp.df[[1]]$mse2
  }
  #take the mean after summing
  res_b$mse = res_b$mse/length(gp.df)
  res_b$mse2 = res_b$mse2/length(gp.df)
  
  res_b$logmse = log(res_b$mse)
  res_b$logmse2 = log(res_b$mse2)
  
  return(res_b)
}

# ------------------------------------------------------------------------------

# Function to Compute growth and MSE over param ranges

compute_pred = function(x,w,f,cases,lower = 7,FUN = mean){
  print(deparse(substitute(x)))
  #check if observations to few
  if(nrow(x)<=lower | nrow(x) == nrow(cases)){
    return(NULL)
  }
  
  #compute growth rate for as the max over the past w_loop days
  #Take max of the first obersvation to keep indexs within positive range
  g =  FUN(exp(x$gr[max(1,(length(x$gr)+1-w)):length(x$gr)]))
  
  
  # For the last day in GP dataset forecast values for the next f
  l_est = exp(x$`log(OD)`)[nrow(x)]
  forecasts = l_est * (g ^ (1:f))
  
  forecast.index = (nrow(x) + 1):(nrow(x) + f)
  
  ret.df = data.frame("Forecasts" = forecasts,
                      "Index" = forecast.index)
  
  return(ret.df)
}

# ------------------------------------------------------------------------------

# Get Dataframe for plotting projections 

visualise.forecast = function(gpList,win,f,FUN = mean){
  names(gpList) = parse_number(names(gpList))
  #compute mean growth rate and MSE for each GP fit period
  gp.forecast = lapply(gpList, compute_pred,win,f,raw_cases_input,7,max)
  gp.forecast[sapply(gp.forecast, is.null)] <- NULL
  
  #res.long.mean = exp.res.df(win_r,f_range,gp.perf.mean)
  
  
  forecast.df = data.frame("Forecasts" =numeric(), "Index" =numeric(),"Group" = factor())
  
  for (i in 1:length(gp.forecast)) {
    loop.df = gp.forecast[[i]]
    loop.df$Group = rep(names(gp.forecast)[i],nrow(loop.df))
    forecast.df = rbind(forecast.df,loop.df)
  }
  
  forecast.df$Data = rep("Prediction",nrow(forecast.df))
  
  true.dat = data.frame("Forecasts" = raw_cases_input$cases, "Index" = raw_cases_input$date,
                        "Group" = rep("True",nrow(raw_cases_input)), "Data" = rep("True",nrow(raw_cases_input)))
  
  comd.df = rbind(forecast.df,true.dat)
  
  return(comd.df)
  
}


################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################


##
##  1. Get Data
##

raw_cases_input <- read_csv("./Python GP Code/Data/y_true_20_04_2020.csv")[,c(1,3)]
colnames(raw_cases_input) = c("date","cases")



# ------------------------------------------------------------------------------

##
##  2. Load in GPs

temp = list.files(pattern="*.csv")
temp = temp[seq(1,length(temp),2)] # only read in main data
gpList = lapply(temp, read_csv)
names(gpList ) = temp


# ------------------------------------------------------------------------------

##
##  3. Grid Search over params and MSE and MSE2
##

win_r <- 2:16 # Window Range for estimating growth rate
f_range <- 4:14 # forecast days


#compute mean growth rate and MSE for each GP fit period
gp.perf.mean = lapply(gpList, compute_error,win_r,f_range,raw_cases_input,7,mean)
gp.perf.mean[sapply(gp.perf.mean, is.null)] <- NULL
res.long.mean = exp.res.df(win_r,f_range,gp.perf.mean)

grid.mean1 = levelplot(logmse ~ window*forecast, 
                       data=res.long.mean,
                       xlab="Window Length",
                       ylab="Forecast Period",
                       main="Natural Log of MSE for Param Grid Mean")


grid.mean2 = levelplot(logmse2 ~ window*forecast, 
                       data=res.long.mean,
                       xlab="Window Length",
                       ylab="Forecast Period",
                       main="Natural Log of MSE2 for Param Grid Mean")


#compute max growth rate and MSE for each GP fit period
gp.perf.max = lapply(gpList, compute_error,win_r,f_range,raw_cases_input,7,max)
gp.perf.max[sapply(gp.perf.max, is.null)] <- NULL
res.long.max = exp.res.df(win_r,f_range,gp.perf.max)

## Try lattice plot
grid.max1 = levelplot(logmse ~ window*forecast, 
                      data=res.long.max,
                      xlab="Window Length",
                      ylab="Forecast Period",
                      main="Natural Log of MSE for Param Grid Max")


grid.max2 = levelplot(logmse2 ~ window*forecast, 
                      data=res.long.max,
                      xlab="Window Length",
                      ylab="Forecast Period",
                      main="Natural Log of MSE2 for Param Grid Max")



#compute min growth rate and MSE for each GP fit period
gp.perf.min = lapply(gpList, compute_error,win_r,f_range,raw_cases_input,7,min)
gp.perf.min[sapply(gp.perf.min, is.null)] <- NULL
res.long.min = exp.res.df(win_r,f_range,gp.perf.min)

## Try lattice plot
grid.min1 = levelplot(logmse ~ window*forecast, 
                      data=res.long.min,
                      xlab="Window Length",
                      ylab="Forecast Period",
                      main="Natural Log of MSE for Param Grid Min")


grid.min2 = levelplot(logmse2 ~ window*forecast, 
                      data=res.long.min,
                      xlab="Window Length",
                      ylab="Forecast Period",
                      main="Natural Log of MSE2 for Param Grid Min")


library(gridExtra)
grid.arrange(grid.mean1, grid.mean2,grid.max1,grid.max2,grid.min1,grid.min2,  nrow = 3)



# ------------------------------------------------------------------------------

##
##  3. Visualise Projection
##

# Change point shapes and colors
forecast.plot.mean = ggplot(visualise.forecast(gpList,win = 7, f = 7, FUN = mean), aes(x=Index, y=Forecasts, group=Group, colour = Data)) +
  geom_line(size = 1.2)+ 
  #scale_y_log10() +
  theme_bw() +
  labs(title = "Daily New Inpatients GP Mean",
       x = "t", y = "Case Counts")



forecast.plot.max = ggplot(visualise.forecast(gpList,win = 7, f = 7, FUN = max), aes(x=Index, y=Forecasts, group=Group, colour = Data)) +
  geom_line(size = 1.2)+ 
  #scale_y_log10() +
  theme_bw() +
  labs(title = "Daily New Inpatients GP Max",
       x = "t", y = "Case Counts")



forecast.plot.min = ggplot(visualise.forecast(gpList,win = 7, f = 7, FUN = min), aes(x=Index, y=Forecasts, group=Group, colour = Data)) +
  geom_line(size = 1.2)+ 
  #scale_y_log10() +
  theme_bw() +
  labs(title = "Daily New Inpatients GP Min",
       x = "t", y = "Case Counts")


grid.arrange(forecast.plot.mean, forecast.plot.max,forecast.plot.min,  nrow = 3)




# ------------------------------------------------------------------------------

##
##  4. 




# Load the lattice package


# Dummy data



# ------------------------------------------------------------------------------

##
##  5. 

stop("End of Script")






