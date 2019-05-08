clear
clc

% This is a demo of MIDAS quantile regression

% load the price data and form returns
P = xlsread('DataQuantile.xlsx','Sheet1','B2:B4437');
[~,dates] = xlsread('DataQuantile.xlsx','Sheet1','A2:A4437');
dates = datenum(dates);
y = log(P(2:end)./P(1:end-1));
dates = dates(2:end);

% Call the function with default settings
[estParams,condQuantile1,yLowFreq,xHighFreq,yDates] = MidasQuantile(y,'Dates',dates); %#ok<*ASGLU>
subplot(2,2,1);plot(yDates,condQuantile1);dateaxis;

% Estimate quantiles for quarterly returns instead of monthly returns
[estParams,condQuantile2,yLowFreq,xHighFreq,yDates] = MidasQuantile(y,'Dates',dates,'Period',66);
subplot(2,2,2);plot(yDates,condQuantile2);dateaxis;

% Fit 25% conditional quantile
[estParams,condQuantile3,yLowFreq,xHighFreq,yDates] = MidasQuantile(y,'Dates',dates,'Period',66,'Quantile',0.25); 
subplot(2,2,3);plot(yDates,condQuantile3);dateaxis;

% Estimation and forecast
yEst = y(1:4200);
estParams = MidasQuantile(yEst);
[~,condQuantile4,yLowFreq,xHighFreq,yDates] = MidasQuantile(y,'Dates',dates,'Params',estParams);
subplot(2,2,4);plot(yDates,condQuantile4);dateaxis;
