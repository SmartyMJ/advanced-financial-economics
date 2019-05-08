clear
clc

% An application of forecast combination

% Quarterly Data from St. Louis FRED website w/ Seas Adj, Real GDP
% Real GDP,  Log - Quarterly First Difference                       
% Monthly Data from FRED, Total Employment Non-Farm Payrolls        
% Total Non-Farm Payrolls, Log - Monthly First Difference           
% This code runs various examples from                              
% Armesto, M.T., K.M. Engemann, and M.T. Owyang, 2010, Forecasting  
% with mixed frequencies, Federal Reserve Bank of St. Louis Review  
% 92, 521?36.

% Specify lag structure and sample size 
Xlag = 12;
Ylag = 1;
Horizon = 3;
EstStart = '1975-07-01';
EstEnd = '2009-07-01';
Method = 'recursive';

% Load first dataset
[DataY,DataYdate] = xlsread('mydata.xlsx','sheet1');
DataYdate = DataYdate(2:end,1);
[DataX,DataXdate] = xlsread('mydata.xlsx','sheet2');
DataXdate = DataXdate(2:end,1);

DataXgrowth = log(DataX(2:end)./DataX(1:end-1))*100;
DataYgrowth = log(DataY(2:end)./DataY(1:end-1))*100;
DataX=DataXgrowth;
DataY=DataYgrowth;
DataYdate = DataYdate(2:end);
DataXdate = DataXdate(2:end);

% Forecast with Beta polynomials
OutputForecast1 = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);

% Load second dataset
[DataY,DataYdate] = xlsread('mydata.xlsx','sheet1');
DataYdate = DataYdate(2:end,1);
[DataX,DataXdate] = xlsread('mydata.xlsx','sheet3');
DataXdate = DataXdate(2:end,1);

DataXgrowth = log(DataX(2:end)./DataX(1:end-1))*100;
DataYgrowth = log(DataY(2:end)./DataY(1:end-1))*100;
DataX=DataXgrowth;
DataY=DataYgrowth;
DataYdate = DataYdate(2:end);
DataXdate = DataXdate(2:end);

% Forecast with Beta polynomials
OutputForecast2 = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);

YfMSFE = ForecastCombine(OutputForecast1,OutputForecast2);
YfDMSFE = ForecastCombine(OutputForecast1,OutputForecast2,'DMSFE');
YfAIC = ForecastCombine({OutputForecast1,OutputForecast2},'aic');
YfBIC = ForecastCombine({OutputForecast1,OutputForecast2},'bic');
YfFlat = ForecastCombine(OutputForecast1,OutputForecast2,'flat');

disp('Forecast by Model 1')
disp(OutputForecast1.Yf')
disp('Forecast by Model 2')
disp(OutputForecast2.Yf')
disp('Combined forecast by MSFE')
disp(YfMSFE')
disp('Combined forecast by DMSFE')
disp(YfDMSFE')
disp('Combined forecast by AIC')
disp(YfAIC')
disp('Combined forecast by BIC')
disp(YfBIC')
disp('Combined forecast by equal weight')
disp(YfFlat')