%Purpose:
%    Econ 525-Spring2019
%Authors:
%    Malik Jabati, Allison Tormey, Malik Jabati — 7May2019
%    UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Note:
%   This .m file is dependent upon the 'ESGdata440.csv', 'Fama French 5
%   factors.csv','440MonPrice.csv', '440MonShare.csv', 'DGS10.csv',
%   and '440Ticker.mat' files.
%Assumptions:
%    Databases = RepRisk database, CRSP returns, Fama-French five factors
%    Assets = Balanced panel of 440 firms from S&P 500 (June 2007 to June 2018)
%    Frequency = Monthly


%% Load in, clean, and prepare data
clear all; close all; clc

load('440Ticker.mat');

ESGdata = readtable('ESGdata440.csv');

ESGdata = removevars(ESGdata,{'RepRisk_ID','ISIN','country_sector_average'});

ESGdata_unstacked = unstack(ESGdata,'current_RRI','Tickers');

Date = cellstr(datestr(ESGdata_unstacked.date));

RRI = table2array(ESGdata_unstacked(:,2:end));

Fivefactors_raw = readtable('Fama French 5 factors.csv');

Mkt_RF = Fivefactors_raw.Mkt_RF;
SMB = Fivefactors_raw.SMB;
HML = Fivefactors_raw.HML;
RMW = Fivefactors_raw.RMW;
CMA = Fivefactors_raw.CMA;

Rf_raw = readtable('DGS10.csv');

%Set risk-free rate
Rf = Rf_raw.DGS10 ./ 100;

rawpricedata = readtable('440MonPrice.csv', 'TreatAsEmpty', '#N/A');

prices = table2array(rawpricedata(:,2:end));


returns = tick2ret(prices);

rawsharedata = readtable('440MonShare.csv');

prices140 = prices(2:end,:);

shares = table2array(rawsharedata(2:end,2:end));

mktcap = prices140 .* shares;

sizetable = array2table(mktcap,'VariableNames',Ticker_unq,'RowNames',Date);

excessreturns = returns - Rf;

%Construct statistics table of stock returns
STOCKmean = nanmean(returns, 'all');
STOCKstd = nanstd(returns, 1, 'all');
STOCKskw = skewness(returns, 1, 'all');
STOCKkrt = kurtosis(returns, 1, 'all');
STOCKmin = min(returns, [], 'all');
STOCKq25 = prctile(returns, 25, 'all');
STOCKq50 = prctile(returns, 50, 'all');
STOCKq75 = prctile(returns, 75, 'all');
STOCKmax = max(returns,[], 'all');

STOCK_stats_table = array2table([STOCKmean; STOCKstd; STOCKskw; STOCKkrt;...
    STOCKmin; STOCKq25; STOCKq50; STOCKq75; STOCKmax],'VariableNames',{'Stocks'});
STOCK_stats_table.Properties.RowNames = {'Mean','Std','Skew','Kurt','Min','Q25','Q50','Q75','Max'};
disp(STOCK_stats_table);

%Construct statistics table of RRI scores
RRImean = nanmean(RRI, 'all');
RRIstd = nanstd(RRI, 1, 'all');
RRIskw = skewness(RRI, 1, 'all');
RRIkrt = kurtosis(RRI, 1, 'all');
RRImin = min(RRI, [], 'all');
RRIq25 = prctile(RRI, 25, 'all');
RRIq50 = prctile(RRI, 50, 'all');
RRIq75 = prctile(RRI, 75, 'all');
RRImax = max(RRI,[], 'all');

RRI_stats_table = array2table([RRImean; RRIstd; RRIskw; RRIkrt;...
    RRImin; RRIq25; RRIq50; RRIq75; RRImax],'VariableNames',{'RRI'});
RRI_stats_table.Properties.RowNames = {'Mean','Std','Skew','Kurt','Min','Q25','Q50','Q75','Max'};
disp(RRI_stats_table);


%% Construct ESG factor

%idxRRI = zeros(size(currentRRI,1),1);
%for i = 1:size(currentRRI,1);
%    if currentRRI(i,:) <= csaRRI(i,:);
%        idxRRI(i,1) = 1;
%    elseif currentRRI(i,:) > csaRRI(i,:);
%        idxRRI(i,1) = 0;
%    end
%end
%
%idxESG = reshape(idxRRI,[140,440]);
%ESG = zeros(size(idxESG,1),1);
%
%for i = 1:size(idxESG,1)
%    ESGgoodret = returns(i,logical(idxESG(i,:)));
%    ESGgoodret((isnan(ESGgoodret))) = [];
%    ESGbadret = returns(i,logical(~idxESG(i,:)));
%    ESGbadret(isnan(ESGbadret)) = [];
%    ESGwgt = [repelem(1/size(ESGgoodret,2),size(ESGgoodret,2)),repelem(-1/size(ESGbadret,2),size(ESGbadret,2))];
%    ESGassret = [ESGgoodret,ESGbadret];
%    ESG(i,1) = sum(([ESGgoodret,ESGbadret] .* ESGwgt),2);
%end

% Build ESG factor

%RRIhigh = ESG.current_RRI<ESG.country_sector_average;
%RRIhigh = double(RRIhigh);

%RRIhigh = reshape(RRIhigh,[140,440]);
%ESG = zeros(size(RRIhigh,1),1);
%for i = 1:size(RRIhigh,1)
%    ESGgoodret = ret(i,logical(RRIhigh(i,:)));
%    ESGgoodret((isnan(ESGgoodret))) = [];
%    ESGbadret = ret(i,logical(~RRIhigh(i,:)));
%    ESGbadret(isnan(ESGbadret)) = [];
%    ESGwgt = [repelem(1/size(ESGgoodret,2),size(ESGgoodret,2)),repelem(-1/size(ESGbadret,2),size(ESGbadret,2))];
%    ESGassret = [ESGgoodret,ESGbadret];
%    ESG(i,1) = sum(([ESGgoodret,ESGbadret] .* ESGwgt),2);
%end


%% 1. Is ESG related to the covariance of returns?
%cov = (dm_ret*dm_ret')./381;

%calculate eigenvalues/vectors
%[V,D] = eig(cov);

%perform PCA
%[coeff,latent,explained] = pcacov(cov);
%sum90 = sum(explained(1:61,1));
%L = 61;
%VL = V(:,1:61);
%Factors = [Mkt_RF,HML,SMB,RMW,CMA,ESG];

%Calculate Canonical Correlation
%[A,B,r,U,V,stats] = canoncorr(VL,Factors);


%% 2. Is the ESG Factor priced?

% First pass over columns of returns
betas=zeros(6,size(returns,2));
for ii = 1:440
    [b,bint,r,rint,stats] = regress(excessreturns(:,ii),[ones(size(returns,1),1),Mkt_RF,SMB,HML,RMW,CMA,RRI(:,ii)]); %conduct the multivariate regression
    betas(:,ii) = b(2:7,1); %store the betas
end

%warning: X is rank deficient to within machine precision.
betatable = array2table(betas,'VariableNames',Ticker_unq);
betatable.Properties.RowNames = {'Mkt_RF','SMB','HML','RMW','CMA','RRI'};

ESGbetatable = array2table([Ticker_unq,num2cell(transpose(betas(end,:)))],'VariableNames',{'Ticker','beta'});

ESGbeta = cell2mat(ESGbetatable.beta);

ESGbeta2pass = zeros(size(ESGbeta,1),12);

betaprctile10 = prctile(cell2mat(ESGbetatable.beta),10);
betaprctile20 = prctile(cell2mat(ESGbetatable.beta),20);
betaprctile30 = prctile(cell2mat(ESGbetatable.beta),30);
betaprctile40 = prctile(cell2mat(ESGbetatable.beta),40);
betaprctile50 = prctile(cell2mat(ESGbetatable.beta),50);
betaprctile60 = prctile(cell2mat(ESGbetatable.beta),60);
betaprctile70 = prctile(cell2mat(ESGbetatable.beta),70);
betaprctile80 = prctile(cell2mat(ESGbetatable.beta),80);
betaprctile90 = prctile(cell2mat(ESGbetatable.beta),90);

%Second Pass
%EIV technique: double sort size and the estimated beta

%Preset the portfolio cell arrays
Port1 = (zeros(size(RRI,2),1));
Port2 = (zeros(size(RRI,2),1));
Port3 = (zeros(size(RRI,2),1));
Port4 = (zeros(size(RRI,2),1));
Port5 = (zeros(size(RRI,2),1));
Port6 = (zeros(size(RRI,2),1));
Port7 = (zeros(size(RRI,2),1));
Port8 = (zeros(size(RRI,2),1));
Port9 = (zeros(size(RRI,2),1));
Port10 = (zeros(size(RRI,2),1));
Port11 = (zeros(size(RRI,2),1));
Port12 = (zeros(size(RRI,2),1));
Port13 = (zeros(size(RRI,2),1));
Port14 = (zeros(size(RRI,2),1));
Port15 = (zeros(size(RRI,2),1));
Port16 = (zeros(size(RRI,2),1));
Port17 = (zeros(size(RRI,2),1));
Port18 = (zeros(size(RRI,2),1));
Port19 = (zeros(size(RRI,2),1));
Port20 = (zeros(size(RRI,2),1));
Port21 = (zeros(size(RRI,2),1));
Port22 = (zeros(size(RRI,2),1));
Port23 = (zeros(size(RRI,2),1));
Port24 = (zeros(size(RRI,2),1));
Port25 = (zeros(size(RRI,2),1));
Port26 = (zeros(size(RRI,2),1));
Port27 = (zeros(size(RRI,2),1));
Port28 = (zeros(size(RRI,2),1));
Port29 = (zeros(size(RRI,2),1));
Port30 = (zeros(size(RRI,2),1));
Port31 = (zeros(size(RRI,2),1));
Port32 = (zeros(size(RRI,2),1));
Port33 = (zeros(size(RRI,2),1));
Port34 = (zeros(size(RRI,2),1));
Port35 = (zeros(size(RRI,2),1));
Port36 = (zeros(size(RRI,2),1));
Port37 = (zeros(size(RRI,2),1));
Port38 = (zeros(size(RRI,2),1));
Port39 = (zeros(size(RRI,2),1));
Port40 = (zeros(size(RRI,2),1));
Port41 = (zeros(size(RRI,2),1));
Port42 = (zeros(size(RRI,2),1));
Port43 = (zeros(size(RRI,2),1));
Port44 = (zeros(size(RRI,2),1));
Port45 = (zeros(size(RRI,2),1));
Port46 = (zeros(size(RRI,2),1));
Port47 = (zeros(size(RRI,2),1));
Port48 = (zeros(size(RRI,2),1));
Port49 = (zeros(size(RRI,2),1));
Port50 = (zeros(size(RRI,2),1));
Port51 = (zeros(size(RRI,2),1));
Port52 = (zeros(size(RRI,2),1));
Port53 = (zeros(size(RRI,2),1));
Port54 = (zeros(size(RRI,2),1));
Port55 = (zeros(size(RRI,2),1));
Port56 = (zeros(size(RRI,2),1));
Port57 = (zeros(size(RRI,2),1));
Port58 = (zeros(size(RRI,2),1));
Port59 = (zeros(size(RRI,2),1));
Port60 = (zeros(size(RRI,2),1));
Port61 = (zeros(size(RRI,2),1));
Port62 = (zeros(size(RRI,2),1));
Port63 = (zeros(size(RRI,2),1));
Port64 = (zeros(size(RRI,2),1));
Port65 = (zeros(size(RRI,2),1));
Port66 = (zeros(size(RRI,2),1));
Port67 = (zeros(size(RRI,2),1));
Port68 = (zeros(size(RRI,2),1));
Port69 = (zeros(size(RRI,2),1));
Port70 = (zeros(size(RRI,2),1));
Port71 = (zeros(size(RRI,2),1));
Port72 = (zeros(size(RRI,2),1));
Port73 = (zeros(size(RRI,2),1));
Port74 = (zeros(size(RRI,2),1));
Port75 = (zeros(size(RRI,2),1));
Port76 = (zeros(size(RRI,2),1));
Port77 = (zeros(size(RRI,2),1));
Port78 = (zeros(size(RRI,2),1));
Port79 = (zeros(size(RRI,2),1));
Port80 = (zeros(size(RRI,2),1));
Port81 = (zeros(size(RRI,2),1));
Port82 = (zeros(size(RRI,2),1));
Port83 = (zeros(size(RRI,2),1));
Port84 = (zeros(size(RRI,2),1));
Port85 = (zeros(size(RRI,2),1));
Port86 = (zeros(size(RRI,2),1));
Port87 = (zeros(size(RRI,2),1));
Port88 = (zeros(size(RRI,2),1));
Port89 = (zeros(size(RRI,2),1));
Port90 = (zeros(size(RRI,2),1));
Port91 = (zeros(size(RRI,2),1));
Port92 = (zeros(size(RRI,2),1));
Port93 = (zeros(size(RRI,2),1));
Port94 = (zeros(size(RRI,2),1));
Port95 = (zeros(size(RRI,2),1));
Port96 = (zeros(size(RRI,2),1));
Port97 = (zeros(size(RRI,2),1));
Port98 = (zeros(size(RRI,2),1));
Port99 = (zeros(size(RRI,2),1));
Port100 = (zeros(size(RRI,2),1));


% Second pass over columns of betas (10x10 portfolios)
size_June = mktcap(6:12:138,:);
for j = 1:size(size_June,1)
    size_June_table = array2table([Ticker_unq,num2cell(transpose(size_June(j,:)))],'VariableNames',{'Ticker','size'});
    size_temp = cell2mat(size_June_table.size);
    
    sizeprctile10 = prctile(cell2mat(size_June_table.size),10);
    sizeprctile20 = prctile(cell2mat(size_June_table.size),20);
    sizeprctile30 = prctile(cell2mat(size_June_table.size),30);
    sizeprctile40 = prctile(cell2mat(size_June_table.size),40);
    sizeprctile50 = prctile(cell2mat(size_June_table.size),50);
    sizeprctile60 = prctile(cell2mat(size_June_table.size),60);
    sizeprctile70 = prctile(cell2mat(size_June_table.size),70);
    sizeprctile80 = prctile(cell2mat(size_June_table.size),80);
    sizeprctile90 = prctile(cell2mat(size_June_table.size),90);
    
    for k = 1:size(mktcap,2)
        if (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile90)
            Port1(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile90)
            Port2(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile90)
            Port3(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile90)
            Port4(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile90)
            Port5(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile90)
            Port6(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile90)
            Port7(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile90)
            Port8(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile90)
            Port9(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile90)
            Port10(k,1) = 1;
            
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port11(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port12(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port13(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port14(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port15(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port16(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port17(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port18(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port19(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile80) && (ESGbeta(k,1) < betaprctile90)
            Port20(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port21(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port22(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port23(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port24(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port25(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port26(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port27(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port28(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port29(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile70) && (ESGbeta(k,1) < betaprctile80)
            Port30(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port31(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port32(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port33(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port34(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port35(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port36(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port37(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port38(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port39(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile60) && (ESGbeta(k,1) < betaprctile70)
            Port40(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port41(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port42(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port43(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port44(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port45(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port46(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port47(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port48(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port49(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile50) && (ESGbeta(k,1) < betaprctile60)
            Port50(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port51(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port52(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port53(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port54(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port55(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port56(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port57(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port58(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port59(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile40) && (ESGbeta(k,1) < betaprctile50)
            Port60(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port61(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port62(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port63(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port64(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port65(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port66(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port67(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port68(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port69(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile30) && (ESGbeta(k,1) < betaprctile40)
            Port70(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port71(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port72(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port73(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port74(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port75(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port76(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port77(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port78(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port79(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile20) && (ESGbeta(k,1) < betaprctile30)
            Port80(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port81(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port82(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port83(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port84(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port85(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port86(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port87(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port88(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port89(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) >= betaprctile10) && (ESGbeta(k,1) < betaprctile20)
            Port90(k,1) = 1;
            
        elseif (size_temp(k,1) >= sizeprctile90) &&  (ESGbeta(k,1) < betaprctile10)
            Port91(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile80) && (size_temp(k) < sizeprctile90) && (ESGbeta(k,1) < betaprctile10)
            Port92(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile70) && (size_temp(k) < sizeprctile80) && (ESGbeta(k,1) < betaprctile10)
            Port93(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile60) && (size_temp(k) < sizeprctile70) && (ESGbeta(k,1) < betaprctile10)
            Port94(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile50) && (size_temp(k) < sizeprctile60) && (ESGbeta(k,1) < betaprctile10)
            Port95(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile40) && (size_temp(k) < sizeprctile50) && (ESGbeta(k,1) < betaprctile10)
            Port96(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile30) && (size_temp(k) < sizeprctile40) && (ESGbeta(k,1) < betaprctile10)
            Port97(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile20) && (size_temp(k) < sizeprctile30) && (ESGbeta(k,1) < betaprctile10)
            Port98(k,1) = 1;
        elseif (size_temp(k) >= sizeprctile10) && (size_temp(k) < sizeprctile20) && (ESGbeta(k,1) < betaprctile10)
            Port99(k,1) = 1;
        elseif (size_temp(k) < sizeprctile10) && (ESGbeta(k,1) < betaprctile10)
            Port100(k,1) = 1;
            
            
        end
    end
    
    %Compute statistics
    Port1mean = mean(ESGbeta(logical(Port1),1));
    ESGbeta2pass(logical(Port1),j) = Port1mean;
    
    Port2mean = mean(ESGbeta(logical(Port2),1));
    ESGbeta2pass(logical(Port2),j) = Port2mean;
    
    Port3mean = mean(ESGbeta(logical(Port3),1));
    ESGbeta2pass(logical(Port3),j) = Port3mean;
    
    Port4mean = mean(ESGbeta(logical(Port4),1));
    ESGbeta2pass(logical(Port4),j) = Port4mean;
    
    Port5mean = mean(ESGbeta(logical(Port5),1));
    ESGbeta2pass(logical(Port5),j) = Port5mean;
    
    Port6mean = mean(ESGbeta(logical(Port6),1));
    ESGbeta2pass(logical(Port6),j) = Port6mean;
    
    Port7mean = mean(ESGbeta(logical(Port7),1));
    ESGbeta2pass(logical(Port7),j) = Port7mean;
    
    Port8mean = mean(ESGbeta(logical(Port8),1));
    ESGbeta2pass(logical(Port8),j) = Port8mean;
    
    Port9mean = mean(ESGbeta(logical(Port9),1));
    ESGbeta2pass(logical(Port9),j) = Port9mean;
    
    Port10mean = mean(ESGbeta(logical(Port10),1));
    ESGbeta2pass(logical(Port10),j) = Port10mean;
    
    
    Port11mean = mean(ESGbeta(logical(Port11),1));
    ESGbeta2pass(logical(Port11),j) = Port11mean;
    
    Port12mean = mean(ESGbeta(logical(Port12),1));
    ESGbeta2pass(logical(Port12),j) = Port12mean;
    
    Port13mean = mean(ESGbeta(logical(Port13),1));
    ESGbeta2pass(logical(Port13),j) = Port13mean;
    
    Port14mean = mean(ESGbeta(logical(Port14),1));
    ESGbeta2pass(logical(Port14),j) = Port14mean;
    
    Port15mean = mean(ESGbeta(logical(Port15),1));
    ESGbeta2pass(logical(Port15),j) = Port15mean;
    
    Port16mean = mean(ESGbeta(logical(Port16),1));
    ESGbeta2pass(logical(Port16),j) = Port16mean;
    
    Port17mean = mean(ESGbeta(logical(Port17),1));
    ESGbeta2pass(logical(Port17),j) = Port17mean;
    
    Port18mean = mean(ESGbeta(logical(Port18),1));
    ESGbeta2pass(logical(Port18),j) = Port18mean;
    
    Port19mean = mean(ESGbeta(logical(Port19),1));
    ESGbeta2pass(logical(Port19),j) = Port19mean;
    
    Port20mean = mean(ESGbeta(logical(Port20),1));
    ESGbeta2pass(logical(Port20),j) = Port20mean;
    
    
    Port21mean = mean(ESGbeta(logical(Port21),1));
    ESGbeta2pass(logical(Port21),j) = Port21mean;
    
    Port22mean = mean(ESGbeta(logical(Port22),1));
    ESGbeta2pass(logical(Port22),j) = Port22mean;
    
    Port23mean = mean(ESGbeta(logical(Port23),1));
    ESGbeta2pass(logical(Port23),j) = Port23mean;
    
    Port24mean = mean(ESGbeta(logical(Port24),1));
    ESGbeta2pass(logical(Port24),j) = Port24mean;
    
    Port25mean = mean(ESGbeta(logical(Port25),1));
    ESGbeta2pass(logical(Port25),j) = Port25mean;
    
    Port26mean = mean(ESGbeta(logical(Port26),1));
    ESGbeta2pass(logical(Port26),j) = Port26mean;
    
    Port27mean = mean(ESGbeta(logical(Port27),1));
    ESGbeta2pass(logical(Port27),j) = Port27mean;
    
    Port28mean = mean(ESGbeta(logical(Port28),1));
    ESGbeta2pass(logical(Port28),j) = Port28mean;
    
    Port29mean = mean(ESGbeta(logical(Port29),1));
    ESGbeta2pass(logical(Port29),j) = Port29mean;
    
    Port30mean = mean(ESGbeta(logical(Port30),1));
    ESGbeta2pass(logical(Port30),j) = Port30mean;
    
    
    Port31mean = mean(ESGbeta(logical(Port31),1));
    ESGbeta2pass(logical(Port31),j) = Port31mean;
    
    Port32mean = mean(ESGbeta(logical(Port32),1));
    ESGbeta2pass(logical(Port32),j) = Port32mean;
    
    Port33mean = mean(ESGbeta(logical(Port33),1));
    ESGbeta2pass(logical(Port33),j) = Port33mean;
    
    Port34mean = mean(ESGbeta(logical(Port34),1));
    ESGbeta2pass(logical(Port34),j) = Port34mean;
    
    Port35mean = mean(ESGbeta(logical(Port35),1));
    ESGbeta2pass(logical(Port35),j) = Port35mean;
    
    Port36mean = mean(ESGbeta(logical(Port36),1));
    ESGbeta2pass(logical(Port36),j) = Port36mean;
    
    Port37mean = mean(ESGbeta(logical(Port37),1));
    ESGbeta2pass(logical(Port37),j) = Port37mean;
    
    Port38mean = mean(ESGbeta(logical(Port38),1));
    ESGbeta2pass(logical(Port38),j) = Port38mean;
    
    Port39mean = mean(ESGbeta(logical(Port39),1));
    ESGbeta2pass(logical(Port39),j) = Port39mean;
    
    Port40mean = mean(ESGbeta(logical(Port40),1));
    ESGbeta2pass(logical(Port40),j) = Port40mean;
    
    
    Port41mean = mean(ESGbeta(logical(Port41),1));
    ESGbeta2pass(logical(Port41),j) = Port41mean;
    
    Port42mean = mean(ESGbeta(logical(Port42),1));
    ESGbeta2pass(logical(Port42),j) = Port42mean;
    
    Port43mean = mean(ESGbeta(logical(Port43),1));
    ESGbeta2pass(logical(Port43),j) = Port43mean;
    
    Port44mean = mean(ESGbeta(logical(Port44),1));
    ESGbeta2pass(logical(Port44),j) = Port44mean;
    
    Port45mean = mean(ESGbeta(logical(Port45),1));
    ESGbeta2pass(logical(Port45),j) = Port45mean;
    
    Port46mean = mean(ESGbeta(logical(Port46),1));
    ESGbeta2pass(logical(Port46),j) = Port46mean;
    
    Port47mean = mean(ESGbeta(logical(Port47),1));
    ESGbeta2pass(logical(Port47),j) = Port47mean;
    
    Port48mean = mean(ESGbeta(logical(Port48),1));
    ESGbeta2pass(logical(Port48),j) = Port48mean;
    
    Port49mean = mean(ESGbeta(logical(Port49),1));
    ESGbeta2pass(logical(Port49),j) = Port49mean;
    
    Port50mean = mean(ESGbeta(logical(Port50),1));
    ESGbeta2pass(logical(Port50),j) = Port50mean;
    
    
    Port51mean = mean(ESGbeta(logical(Port51),1));
    ESGbeta2pass(logical(Port51),j) = Port51mean;
    
    Port52mean = mean(ESGbeta(logical(Port52),1));
    ESGbeta2pass(logical(Port52),j) = Port52mean;
    
    Port53mean = mean(ESGbeta(logical(Port53),1));
    ESGbeta2pass(logical(Port53),j) = Port53mean;
    
    Port54mean = mean(ESGbeta(logical(Port54),1));
    ESGbeta2pass(logical(Port54),j) = Port54mean;
    
    Port55mean = mean(ESGbeta(logical(Port55),1));
    ESGbeta2pass(logical(Port55),j) = Port55mean;
    
    Port56mean = mean(ESGbeta(logical(Port56),1));
    ESGbeta2pass(logical(Port56),j) = Port56mean;
    
    Port57mean = mean(ESGbeta(logical(Port57),1));
    ESGbeta2pass(logical(Port57),j) = Port57mean;
    
    Port58mean = mean(ESGbeta(logical(Port58),1));
    ESGbeta2pass(logical(Port58),j) = Port58mean;
    
    Port59mean = mean(ESGbeta(logical(Port59),1));
    ESGbeta2pass(logical(Port59),j) = Port59mean;
    
    Port60mean = mean(ESGbeta(logical(Port60),1));
    ESGbeta2pass(logical(Port60),j) = Port60mean;
    
    
    Port61mean = mean(ESGbeta(logical(Port61),1));
    ESGbeta2pass(logical(Port61),j) = Port61mean;
    
    Port62mean = mean(ESGbeta(logical(Port62),1));
    ESGbeta2pass(logical(Port62),j) = Port62mean;
    
    Port63mean = mean(ESGbeta(logical(Port63),1));
    ESGbeta2pass(logical(Port63),j) = Port63mean;
    
    Port64mean = mean(ESGbeta(logical(Port64),1));
    ESGbeta2pass(logical(Port64),j) = Port64mean;
    
    Port65mean = mean(ESGbeta(logical(Port65),1));
    ESGbeta2pass(logical(Port65),j) = Port65mean;
    
    Port66mean = mean(ESGbeta(logical(Port66),1));
    ESGbeta2pass(logical(Port66),j) = Port66mean;
    
    Port67mean = mean(ESGbeta(logical(Port67),1));
    ESGbeta2pass(logical(Port67),j) = Port67mean;
    
    Port68mean = mean(ESGbeta(logical(Port68),1));
    ESGbeta2pass(logical(Port68),j) = Port68mean;
    
    Port69mean = mean(ESGbeta(logical(Port69),1));
    ESGbeta2pass(logical(Port69),j) = Port69mean;
    
    Port70mean = mean(ESGbeta(logical(Port70),1));
    ESGbeta2pass(logical(Port70),j) = Port70mean;
    
    
    Port71mean = mean(ESGbeta(logical(Port71),1));
    ESGbeta2pass(logical(Port71),j) = Port71mean;
    
    Port72mean = mean(ESGbeta(logical(Port72),1));
    ESGbeta2pass(logical(Port72),j) = Port72mean;
    
    Port73mean = mean(ESGbeta(logical(Port73),1));
    ESGbeta2pass(logical(Port73),j) = Port73mean;
    
    Port74mean = mean(ESGbeta(logical(Port74),1));
    ESGbeta2pass(logical(Port74),j) = Port74mean;
    
    Port75mean = mean(ESGbeta(logical(Port75),1));
    ESGbeta2pass(logical(Port75),j) = Port75mean;
    
    Port76mean = mean(ESGbeta(logical(Port76),1));
    ESGbeta2pass(logical(Port76),j) = Port76mean;
    
    Port77mean = mean(ESGbeta(logical(Port77),1));
    ESGbeta2pass(logical(Port77),j) = Port77mean;
    
    Port78mean = mean(ESGbeta(logical(Port78),1));
    ESGbeta2pass(logical(Port78),j) = Port78mean;
    
    Port79mean = mean(ESGbeta(logical(Port79),1));
    ESGbeta2pass(logical(Port79),j) = Port79mean;
    
    Port80mean = mean(ESGbeta(logical(Port80),1));
    ESGbeta2pass(logical(Port80),j) = Port80mean;
    
    
    Port81mean = mean(ESGbeta(logical(Port81),1));
    ESGbeta2pass(logical(Port81),j) = Port81mean;
    
    Port82mean = mean(ESGbeta(logical(Port82),1));
    ESGbeta2pass(logical(Port82),j) = Port82mean;
    
    Port83mean = mean(ESGbeta(logical(Port83),1));
    ESGbeta2pass(logical(Port83),j) = Port83mean;
    
    Port84mean = mean(ESGbeta(logical(Port84),1));
    ESGbeta2pass(logical(Port84),j) = Port84mean;
    
    Port85mean = mean(ESGbeta(logical(Port85),1));
    ESGbeta2pass(logical(Port85),j) = Port85mean;
    
    Port86mean = mean(ESGbeta(logical(Port86),1));
    ESGbeta2pass(logical(Port86),j) = Port86mean;
    
    Port87mean = mean(ESGbeta(logical(Port87),1));
    ESGbeta2pass(logical(Port87),j) = Port87mean;
    
    Port88mean = mean(ESGbeta(logical(Port88),1));
    ESGbeta2pass(logical(Port88),j) = Port88mean;
    
    Port89mean = mean(ESGbeta(logical(Port89),1));
    ESGbeta2pass(logical(Port89),j) = Port89mean;
    
    Port90mean = mean(ESGbeta(logical(Port90),1));
    ESGbeta2pass(logical(Port90),j) = Port90mean;
    
    
    Port91mean = mean(ESGbeta(logical(Port91),1));
    ESGbeta2pass(logical(Port91),j) = Port91mean;
    
    Port92mean = mean(ESGbeta(logical(Port92),1));
    ESGbeta2pass(logical(Port92),j) = Port92mean;
    
    Port93mean = mean(ESGbeta(logical(Port93),1));
    ESGbeta2pass(logical(Port93),j) = Port93mean;
    
    Port94mean = mean(ESGbeta(logical(Port94),1));
    ESGbeta2pass(logical(Port94),j) = Port94mean;
    
    Port95mean = mean(ESGbeta(logical(Port95),1));
    ESGbeta2pass(logical(Port95),j) = Port95mean;
    
    Port96mean = mean(ESGbeta(logical(Port96),1));
    ESGbeta2pass(logical(Port96),j) = Port96mean;
    
    Port97mean = mean(ESGbeta(logical(Port97),1));
    ESGbeta2pass(logical(Port97),j) = Port97mean;
    
    Port98mean = mean(ESGbeta(logical(Port98),1));
    ESGbeta2pass(logical(Port98),j) = Port98mean;
    
    Port99mean = mean(ESGbeta(logical(Port99),1));
    ESGbeta2pass(logical(Port99),j) = Port99mean;
    
    Port100mean = mean(ESGbeta(logical(Port100),1));
    ESGbeta2pass(logical(Port100),j) = Port100mean;
    
    
end

ESGbeta2pass = transpose(ESGbeta2pass);

beta2pass = zeros(size(RRI,1)+4,size(RRI,2));
for l = 1:size(RRI,2)
    beta2pass(:,l) = repelem(ESGbeta2pass(:,l),12);
end

beta2pass = beta2pass(1:end-4,:);
beta2pass = transpose(beta2pass);

lambdas = zeros(6,size(returns,2));
stats2pass = zeros(6,size(returns,2));
Mkt_Rfbetas = transpose(betas(1,:));
SMBbetas = transpose(betas(2,:));
HMLbetas = transpose(betas(3,:));
RMWbetas = transpose(betas(4,:));
CMAbetas = transpose(betas(5,:));

avgret = transpose(mean(excessreturns));

for ii = 1:140
    [b,bint,r,rint,stats] = regress(avgret,[ones(size(avgret,1),1),Mkt_Rfbetas,SMBbetas,HMLbetas,RMWbetas,CMAbetas,beta2pass(:,ii)]); %conduct the multivariate regression
    lambdas(:,ii) = b(2:end,1); %store the betas
end
%% Reasonable Reward-to-Risk
    
    % Compute Fama-French-style portfolio returns for each factor
    % Compute zero-investment long-only portfolio returns combined with the market for each factor
    
    % Compute Sharpe Ratios for both types of portfolios

%% Construct the descriptive statistics table for the variables
Mkt_RFmean = mean(Mkt_RF);
SMBmean = mean(SMB);
HMLmean = mean(HML);
RMWmean = mean(RMW);
CMAmean = mean(CMA);

Mkt_RFstd = std(Mkt_RF);
SMBstd = std(SMB);
HMLstd = std(HML);
RMWstd = std(RMW);
CMAstd = std(CMA);

Mkt_RFskw = skewness(Mkt_RF);
SMBskw = skewness(SMB);
HMLskw = skewness(HML);
RMWskw = skewness(RMW);
CMAskw = skewness(CMA);

Mkt_RFkrt = kurtosis(Mkt_RF);
SMBkrt = kurtosis(SMB);
HMLkrt = kurtosis(HML);
RMWkrt = kurtosis(RMW);
CMAkrt = kurtosis(CMA);

Mkt_RFmin= min(Mkt_RF);
SMBmin = min(SMB);
HMLmin = min(HML);
RMWmin = min(RMW);
CMAmin = min(CMA);

Mkt_RFQ25= prctile(Mkt_RF,25);
SMBQ25 = prctile(SMB,25);
HMLQ25 = prctile(HML,25);
RMWQ25 = prctile(RMW,25);
CMAQ25 = prctile(CMA,25);

Mkt_RFQ50= prctile(Mkt_RF,50);
SMBQ50 = prctile(SMB,50);
HMLQ50 = prctile(HML,50);
RMWQ50 = prctile(RMW,50);
CMAQ50 = prctile(CMA,50);

Mkt_RFQ75= prctile(Mkt_RF,75);
SMBQ75 = prctile(SMB,75);
HMLQ75 = prctile(HML,75);
RMWQ75 = prctile(RMW,75);
CMAQ75 = prctile(CMA,75);

Mkt_RFmax= max(Mkt_RF);
SMBmax = max(SMB);
HMLmax = max(HML);
RMWmax = max(RMW);
CMAmax = max(CMA);

DescriptiveStats4FiveFactors = array2table([Mkt_RFmean,SMBmean,HMLmean,RMWmean,CMAmean;Mkt_RFstd,SMBstd,HMLstd,RMWstd,CMAstd;Mkt_RFskw,SMBskw,HMLskw,RMWskw,CMAskw;Mkt_RFkrt,SMBkrt,HMLkrt,RMWkrt,CMAkrt;Mkt_RFmin,SMBmin,HMLmin,RMWmin,CMAmin;Mkt_RFQ25,SMBQ25,HMLQ25,RMWQ25,CMAQ25;Mkt_RFQ50,SMBQ50,HMLQ50,RMWQ50,CMAQ50;Mkt_RFQ75,SMBQ75,HMLQ75,RMWQ75,CMAQ75;Mkt_RFmax,SMBmax,HMLmax,RMWmax,CMAmax]);
DescriptiveStats4FiveFactors.Properties.VariableNames = {'Mkt_Rf','SMB','HML','RMW','CMA'};
DescriptiveStats4FiveFactors.Properties.RowNames = {'Mean','Std','Skew','Kurt','Min','Q25','Q50','Q75','Max'};
disp(DescriptiveStats4FiveFactors);

lambdas07 = mean(lambdas(:,1:12),2)';
lambdas08 = mean(lambdas(:,13:24),2)';
lambdas09 = mean(lambdas(:,15:36),2)';
lambdas10 = mean(lambdas(:,37:48),2)';
lambdas11 = mean(lambdas(:,49:60),2)';
lambdas12 = mean(lambdas(:,61:72),2)';
lambdas13 = mean(lambdas(:,73:84),2)';
lambdas14 = mean(lambdas(:,85:96),2)';
lambdas15 = mean(lambdas(:,97:108),2)';
lambdas16 = mean(lambdas(:,109:120),2)';
lambdas17 = mean(lambdas(:,121:132),2)';
lambdas18 = mean(lambdas(:,133:140),2)';

tstat07 = mean(stats2pass(:,1:12),2)';
tstat08 = mean(stats2pass(:,13:24),2)';
tstat09 = mean(stats2pass(:,15:36),2)';
tstat10 = mean(stats2pass(:,37:48),2)';
tstat11 = mean(stats2pass(:,49:60),2)';
tstat12 = mean(stats2pass(:,61:72),2)';
tstat13 = mean(stats2pass(:,73:84),2)';
tstat14 = mean(stats2pass(:,85:96),2)';
tstat15 = mean(stats2pass(:,97:108),2)';
tstat16 = mean(stats2pass(:,109:120),2)';
tstat17 = mean(stats2pass(:,121:132),2)';
tstat18 = mean(stats2pass(:,133:140),2)';

FamaMacBethResults = array2table([lambdas07;tstat07;lambdas08;tstat08;lambdas09;tstat09;lambdas10;tstat10;lambdas11;tstat11;lambdas12;tstat12;lambdas13;tstat13;lambdas14;tstat14;lambdas15;tstat15;lambdas16;tstat16;lambdas17;tstat17;lambdas18;tstat18]);
FamaMacBethResults.Properties.VariableNames = {'Mkt_Rf','SMB','HML','RMW','CMA','ESG'};
FamaMacBethResults.Properties.RowNames = {'gamma1_07','tstat_07','gamma1_08','tstat_08','gamma1_09','tstat_09','gamma1_10','tstat_10','gamma1_11','tstat_11','gamma1_12','tstat_12','gamma1_13','tstat_13','gamma1_14','tstat_14','gamma1_15','tstat_15','gamma1_16','tstat_16','gamma1_17','tstat_17','gamma1_18','tstat_18'};
%disp(FamaMacBethResults);

Mkt_Rfhighb = zeros(size(Mkt_Rfbetas,1),1);
Mkt_Rflowb = zeros(size(Mkt_Rfbetas,1),1);
for i = 1:size(Mkt_Rfbetas,1)
    if Mkt_Rfbetas(i,1) >= prctile(Mkt_Rfbetas,70)
        Mkt_Rfhighb(i,1) = 1;
    elseif Mkt_Rfbetas(i,1) <= prctile(Mkt_Rfbetas,30)
        Mkt_Rflowb(i,1) = 1;
    end
    
end

SMBhighb = zeros(size(SMBbetas,1),1);
SMBlowb = zeros(size(SMBbetas,1),1);
for i = 1:size(SMBbetas,1)
    if SMBbetas(i,1) >= prctile(SMBbetas,70)
        SMBhighb(i,1) = 1;
    elseif SMBbetas(i,1) <= prctile(SMBbetas,30)
        SMBlowb(i,1) = 1;
    end
    
end

HMLhighb = zeros(size(HMLbetas,1),1);
HMLlowb = zeros(size(HMLbetas,1),1);
for i = 1:size(HMLbetas,1)
    if HMLbetas(i,1) >= prctile(HMLbetas,70)
        HMLhighb(i,1) = 1;
    elseif HMLbetas(i,1) <= prctile(HMLbetas,30)
        HMLlowb(i,1) = 1;
    end
    
end

RMWhighb = zeros(size(RMWbetas,1),1);
RMWlowb = zeros(size(RMWbetas,1),1);
for i = 1:size(RMWbetas,1)
    if RMWbetas(i,1) >= prctile(RMWbetas,70)
        RMWhighb(i,1) = 1;
    elseif SMBbetas(i,1) <= prctile(RMWbetas,30)
        RMWlowb(i,1) = 1;
    end
    
end

CMAhighb = zeros(size(CMAbetas,1),1);
CMAlowb = zeros(size(CMAbetas,1),1);
for i = 1:size(CMAbetas,1)
    if CMAbetas(i,1) >= prctile(CMAbetas,70)
        CMAhighb(i,1) = 1;
    elseif CMAbetas(i,1) <= prctile(CMAbetas,30)
        CMAlowb(i,1) = 1;
    end
    
end

ESGhighb = zeros(size(ESGbeta,1),1);
ESGlowb = zeros(size(ESGbeta,1),1);
for i = 1:size(ESGbeta,1)
    if ESGbeta(i,1) >= prctile(ESGbeta,70)
        ESGhighb(i,1) = 1;
    elseif ESGbeta(i,1) <= prctile(ESGbeta,30)
        ESGlowb(i,1) = 1;
    end
    
end

Mkt_Rflowbret = returns(:,logical(Mkt_Rflowb));
Mkt_Rflowbret(:,any(isnan(Mkt_Rflowbret))) = [];

Mkt_Rfhighbret = returns(:,logical(Mkt_Rfhighb));
Mkt_Rfhighbret(:,any(isnan(Mkt_Rfhighbret))) = [];

pm_ret = [Mkt_Rflowbret,Mkt_Rfhighbret];
pm_wgt = [repelem(-1/size(Mkt_Rflowbret,2),size(Mkt_Rflowbret,2)),repelem(1/size(Mkt_Rfhighbret,2),size(Mkt_Rfhighbret,2))];

pm_ret = sum(pm_ret .* pm_wgt,2);
pm_std = std(pm_ret);






