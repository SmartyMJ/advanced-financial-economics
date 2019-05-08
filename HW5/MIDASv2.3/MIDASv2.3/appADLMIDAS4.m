clear
clc

% This demo illustrate how to include exogenous regressors in MIDAS
% The model specification is the same as that in appADLMIDAS1

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
Xlag = 9;
Horizon = 3;
EstStart = '1985-01-01';
EstEnd = '2009-01-01';
Method = 'fixedWindow';

% We would like to include Y(t-1) and Y(t-3) in the regression
Ylag = {1,3};
MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);

% Alternatively, we manually put Y(t-1) and Y(t-3) as exogenous regressors
DataYcut = DataY(4:end);
DataYdateCut = DataYdate(4:end);
ExoReg = [DataY(3:end-1),DataY(1:end-3)];
ExoRegDate = DataYdate(4:end);
MIDAS_ADL(DataYcut,DataYdateCut,DataX,DataXdate,'Xlag',Xlag,'Ylag',0,'ExoReg',ExoReg,'ExoRegDate',ExoRegDate,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);

