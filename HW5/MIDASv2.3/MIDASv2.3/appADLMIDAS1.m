clear
clc

% Using MIDAS to forecast GDP growth with monthly employment growth.
% This application is comparable to the "Example_fixed_forecast"
% in the original MIDAS package.

% Quarterly Data from St. Louis FRED website w/ Seas Adj, Real GDP
% Real GDP,  Log - Quarterly First Difference                       
% Monthly Data from FRED, Total Employment Non-Farm Payrolls        
% Total Non-Farm Payrolls, Log - Monthly First Difference           
% This code runs various examples from                              
% Armesto, M.T., K.M. Engemann, and M.T. Owyang, 2010, Forecasting  
% with mixed frequencies, Federal Reserve Bank of St. Louis Review  
% 92, 521-36.                                                      

% Load data
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

% Specify lag structure and sample size 
Xlag = 3;
Ylag = 2;
% Horizon = 1;
Horizon = 2;
EstStart = '1985-01-01';
EstEnd = '2009-01-01';
Method = 'fixedWindow';

% Forecast with various weight polynomials
[OutputForecast1,OutputEstimate1,MixedFreqData,ExtendedForecast] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
%% Break

[OutputForecast2,OutputEstimate2] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','betaNN','Method',Method,'Display','estimate');
[OutputForecast3,OutputEstimate3] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','expAlmon','Method',Method,'Display','estimate');
[OutputForecast4,OutputEstimate4] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','umidas','Method',Method,'Display','estimate');
[OutputForecast5,OutputEstimate5] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','step','Method',Method,'Display','estimate');
[OutputForecast6,OutputEstimate6] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','Almon','Method',Method,'Display','estimate');

fprintf('RMSE Beta:          %5.4f\n',OutputForecast1.RMSE);
fprintf('RMSE Beta Non-Zero: %5.4f\n',OutputForecast2.RMSE);
fprintf('RMSE Exp Almon:     %5.4f\n',OutputForecast3.RMSE);
fprintf('RMSE U-MIDAS:       %5.4f\n',OutputForecast4.RMSE);
fprintf('RMSE Stepfun:       %5.4f\n',OutputForecast5.RMSE);
fprintf('RMSE Almon:         %5.4f\n',OutputForecast6.RMSE);
disp(' ')

% Generate a plot of weights
Xlag = MixedFreqData.Xlag;
nModel = 6;
nrow = ceil(sqrt(nModel));
ncol = ceil(nModel./nrow);
figure(1)
for m = 1:nModel
    weights = eval(['OutputEstimate',num2str(m),'.estWeights']);
    subplot(nrow,ncol,m);plot(1:Xlag,weights);
    title(eval(['OutputEstimate',num2str(m),'.model(7:end)']));
end

% Adjust the leads and lags of the high frequency variable

% Reg Y(01/01/85) on Y(10/01/84),X(10/01/84),X(09/01/84),...,X(02/01/84)
MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,'EstStart',EstStart,'EstEnd',EstEnd,'Horizon',3);

% Reg Y(01/01/85) on Y(10/01/84),X(11/01/84),X(10/01/84),...,X(03/01/84)
MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,'EstStart',EstStart,'EstEnd',EstEnd,'Horizon',2);

% Reg Y(01/01/85) on Y(10/01/84),X(12/01/84),X(11/01/84),...,X(04/01/84)
MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,'EstStart',EstStart,'EstEnd',EstEnd,'Horizon',1);
