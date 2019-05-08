clear
clc

% Purpose:
%
% This demo replicates the empirical results of Colacito et al. (2011), 
% and illustrates a tri-variate DCC-MIDAS model with stock returns, 
% exchange rates and bond yields.
%
% Reference:
%
% Colacito,R., Engle,R.F. and Ghysels, E. (2011), A Component Model for
% Dynamic Correlations, Journal of Econometrics, 164, 45-59.
%
% Written by Hang Qian
% Contact: matlabist@gmail.com

% In order to reproduce the empirical results of Colacito et al. (2011)
% Users must obtain their original data. The following source may help:
% At the website of Kenneth French, "10 Industry Portfolios", Columns Enrgy, HiTec 
reproduceFlag = exist('energy-hitech-10yrbond_2000toend.txt','file') == 2;

if reproduceFlag
    % Load sample data: Energy, Hi-Tech and Bond portfolio, 1971-2006
    Data = load('energy-hitech-10yrbond_2000toend.txt');    
else
    % Load alternative data
    % NASDAQ Composite Index, percentage change, 1971-2006
    % Japan / U.S. Foreign Exchange Rate, percentage change, 1971-2006
    % 10-Year Treasury Constant Maturity Rate, percentage change, 1971-2006
    % Data Source: FRED Economic Data
    % https://research.stlouisfed.org/fred2/series/NASDAQCOM (and DEXJPUS, DGS10)
    y1 = xlsread('NASDAQCOM.xls','B22:B9257');
    y2 = xlsread('DEXJPUS.xls','B14:B9249');
    y3 = xlsread('DGS10.xls','B16:B9251');    
    Data = [y1,y2,y3];
end
nobs = size(Data,1);

% Estimate the DCC-MIDAS model
options = optimoptions('fmincon','Algorithm','active-set');
[estParamsStep1,~,estParamsStep2,~,Variance,LongRunVar,CorrMatrix,LongRunCorrMatrix]...
    = DccMidas(Data,'Period',20,'NumLagsVar',36,'NumLagsCorr',144,'options',options,'ZeroLogL',1:3600,'mu0',0.001);
CorrMatrix = reshape(CorrMatrix,9,nobs)';
LongRunCorrMatrix = reshape(LongRunCorrMatrix,9,nobs)';

% Reproduce Table 1 of Colacito et al. (2011)
if reproduceFlag
    RowNames = {'Energy';'Hi-Tech';'Bond'};
else
    RowNames = {'NASDAQ';'ExchangeRate';'Bond'};
end
VariableNames = {'mu','alpha','beta','theta','w','m'};
Table1a = array2table(estParamsStep1','RowNames',RowNames,'VariableNames',VariableNames);
RowNames = {'DCC-MIDAS'};
VariableNames = {'a','b','w'};
Table1b = array2table(estParamsStep2','RowNames',RowNames,'VariableNames',VariableNames);
disp(Table1a)
disp(Table1b)

% Reproduce Figure 1 of Colacito et al. (2011)
seq = (3601:nobs)';
dates = linspace(1970.5,2006.5,nobs);
figure(1)
subplot(3,3,1); 
plot(dates(seq),Variance(seq,1),'g--','LineWidth',1);
xlim([1985,2006.5])
if reproduceFlag; ylim([0,60]); else ylim([0,40]); end
hold on
plot(dates(seq),LongRunVar(seq,1),'b-','LineWidth',2);
hold off

subplot(3,3,5); 
plot(dates(seq),Variance(seq,2),'g--','LineWidth',1);
xlim([1985,2006.5])
if reproduceFlag; ylim([0,60]); else ylim([0,5]); end
hold on
plot(dates(seq),LongRunVar(seq,2),'b-','LineWidth',2);
hold off

subplot(3,3,9); 
plot(dates(seq),Variance(seq,3),'g--','LineWidth',1);
xlim([1985,2006.5])
if reproduceFlag; ylim([0,4]); end
hold on; plot(dates(seq),LongRunVar(seq,3),'b-','LineWidth',2); hold off

subplot(3,3,2); 
plot(dates(seq),CorrMatrix(seq,2),'g--','LineWidth',1);
xlim([1985,2006.5])
if reproduceFlag; ylim([-0.2,0.8]); end
hold on; plot(dates(seq),LongRunCorrMatrix(seq,2),'b-','LineWidth',2); hold off

subplot(3,3,3); 
plot(dates(seq),CorrMatrix(seq,3),'g--','LineWidth',1);
xlim([1985,2006.5])
if reproduceFlag; ylim([-0.4,0.4]); end
hold on; plot(dates(seq),LongRunCorrMatrix(seq,3),'b-','LineWidth',2); hold off

subplot(3,3,6); 
plot(dates(seq),CorrMatrix(seq,6),'g--','LineWidth',1);
xlim([1985,2006.5])
if reproduceFlag; ylim([-0.5,0.5]); end
hold on; plot(dates(seq),LongRunCorrMatrix(seq,6),'b-','LineWidth',2); hold off


