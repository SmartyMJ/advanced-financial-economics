%%COMMENTS:
%good work

%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon HW_7_Tickers.xlsx and TickerData.xlsx.
%Author:
    %Malik Jabati — 29Mar2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Assumptions:
    %Use the constituents of the DJIA-30 listed in the Excel file provided.
    %Use monthly return data from 01/01/2017 through 02/01/2019.
    %Set your expected returns, variances, and covariances to their historical averages.
    %Portfolio 1: long only, no leverage.
    %Portfolio 2: long only, no leverage, 6 ≤ Np ≤ 15, where Np is the number of assets in the portfolio.
%% Load in data

%All data was cleaned in Excel before importing

%HouseKeeping
    clear all; close all; clc
    
%Import tickers
    [num, txt, raw] = xlsread('HW_7_Tickers.xlsx');
    
    tickers = txt(2:end,1);
    
%Import ticker prices
    [num_p, txt_p, raw_p] = xlsread('TickerData.xlsx','Prices');
    
%Import ticker weights
    [num_w, txt_w, raw_w] = xlsread('TickerData.xlsx','Weights');

%Set dates    
    Dates = x2mdate(num_p(2:end,1));
    Dates = datetime(Dates,'ConvertFrom','datenum');
    
%Returns
    Returns = diff(log(num_p(:,2:end)));
    
%Benchmark weights
    BenchmarkWeights = num_w';

%Enter the API key
%    apikey = 'HxTGtomxL79TZzQg_Ey4'; 

%Establish database connection
%    c = quandl(apikey);
    
%Period of interest
%    startdate = datetime('01-01-2017','InputFormat','MM-dd-yyyy');
%    enddate = datetime('02-01-2019','InputFormat','MM-dd-yyyy');
    
%Identify the database
%    source = 'EOD/';        % End of Day US Stock Prices
%    periodicity = 'monthly';
    
%Create empty array for prices data
%    prices = zeros([26 28]);
%    for i=1:size(prices,2)
%        ticker = tickers(i);
%        s = strcat(source,ticker);
      
%        data = history(c,s,startdate,enddate,periodicity);  %pull data from Quandl
%       ticker_prices = data{:,'Adj_Close'};     %get adjusted close
%        prices(:,i) = ticker_prices;
%    end
    
%    dates = data{:,'Time'};


%% Create portfolios

%Set the market
    Market = Returns(:,1); %column of DJ index
    mret = mean(Market); 
    mrsk = std(Market);
    
%Set Returns
    Assets = Returns(:,2:end);
    Tickers = tickers';
    
%Compute Moments
    AssetsMean = mean(Assets); 
    AssetsCovar = cov(Assets); 
    
%Set Risk Free
    cret =  0.2; %Cash return
    crsk = 0; %cash risk 
    
%Create Portfolio Object
    p = Portfolio('AssetList',Tickers,'RiskFreeRate',cret); 
    p = setAssetMoments(p,AssetsMean,AssetsCovar);
    
%Set Equal Weight
    p = setInitPort(p,1/p.NumAssets);
    [ersk,eret]=estimatePortMoments(p,p.InitPort); %this sets up the function to estimate the mean and stdev of equal-weighted portfolio ret.
    
% Default Portfolio: fully invested, long-only portfolios (non-neg weights and must sum up to 1)
    p = setDefaultConstraints(p); %Set the portfolio constraints
    pwgt = estimateFrontier(p,20); %Estimate the efficient frontier, efficient portfolios
    [prsk,pret]=estimatePortMoments(p,pwgt); %Grab the portfolio moments
    clf; fig1 = figure(1); 
    portfolioexamples_plot('Efficient Frontier', ...
	{'line', prsk, pret}, ...
	{'scatter', [mrsk, ersk], [mret, eret], {'Market', 'Equal'}}, ...
	{'scatter', sqrt(diag(p.AssetCovar)), p.AssetMean, p.AssetList, '.r'});

%% Create Portfolio 1: long only, no leverage.

%Choose efficient portfolio #10, because it's std. dev. (0.0347) is closest to the
%market's (0.0354).

%Make weight matrix
%pwgt_matrix = repelem(pwgt(:,10)',26,1);
%PortReturns_ = sum(Assets.*pwgt_matrix, 2);

PortReturns = Assets*pwgt(:,10); %Recreate portfolio returns
PortPrice = ret2tick(PortReturns./100,100); %Recreate portfolio prices index to 100
fig2 = figure(3);
subplot(1,2,1)
plot(Dates,PortPrice(2:end,1))
ylabel('Portfolio 1 ($) Indexed at 100 on Jan`17')
subplot(1,2,2)
plot(Dates,PortReturns(1:end,1))
ylabel('% Simple Return')
%% 
%% Create Portfolio 2: long only, no leverage, 6 <= N <= 15, where N is the number of assets in the portfolio

%Settings
    p2 = p; %Reset the portfolio to the initial
    p2 = setDefaultConstraints(p2); %Set the default constraints   
    p2 = setMinMaxNumAssets(p2,6,15);   %Set minimum assets to 6 and maximum assets to 15
    p2 = setBounds(p2, .01, 'BoundType', 'conditional'); %Actual DJ30 minimum weight is 1.11%
   
%Efficient portfolios for portfolio 2
    pwgt2 = estimateFrontier(p2,20); %Estimate the efficient frontier, efficient portfolios
    [prsk2,pret2]=estimatePortMoments(p2,pwgt2); %Grab the portfolio moments
    clf; fig3 = figure(1); 
    portfolioexamples_plot('Efficient Frontier', ...
	{'line', prsk2, pret2}, ...
	{'scatter', [mrsk, ersk], [mret, eret], {'Market', 'Equal'}}, ...
	{'scatter', sqrt(diag(p2.AssetCovar)), p2.AssetMean, p2.AssetList, '.r'});
%Choose efficient portfolio #11, because it's std. dev. (0.0366) is closest to the
%market's (0.0354).

PortReturns2 = Assets*pwgt2(:,11); %Recreate portfolio returns
PortPrice2 = ret2tick(PortReturns2./100,100); %Recreate portfolio prices index to 100
fig4 = figure(3);
subplot(1,2,1)
plot(Dates,PortPrice2(2:end,1))
ylabel('Portfolio 2 ($) Indexed at 100 on Jan`17')
subplot(1,2,2)
plot(Dates,PortReturns2(1:end,1))
ylabel('% Simple Return')
%% 1. Fill in the following table. The performance metrics can be computed entirely from the in-sample estimation period.

%Compute active returns
    ActiveReturn = PortReturns - Market;
    ActiveReturn2 = PortReturns2 - Market;


%Total active return
    TotalReturn = sum(ActiveReturn);
    TotalReturn2 = sum(ActiveReturn2);

%Monthly active return
    MonthlyReturn = mean(ActiveReturn);
    MonthlyReturn2 = mean(ActiveReturn2);

%Active risk (i.e. stdev of active return)
    ActiveRisk = std(ActiveReturn);
    ActiveRisk2 = std(ActiveReturn2);

%Information ratio
    InfoRatio = MonthlyReturn / ActiveRisk;
    InfoRatio2 = MonthlyReturn2 / ActiveRisk2;

%MaxDrawdown of active return (in terms of percentage points NOT percentage)
    MaxDD = maxdrawdown(ActiveReturn, 'arithmetic');
    MaxDD2 = maxdrawdown(ActiveReturn2, 'arithmetic');

%Omega Ratio (Active @0)
    Omega = lpm(-ActiveReturn, -0, 1) / lpm(ActiveReturn, 0, 1);
    Omega2 = lpm(-ActiveReturn2, -0, 1) / lpm(ActiveReturn2, 0, 1);

%Sortino Ratio (Active @0)
    Sortino = (mean(ActiveReturn) - 0) / sqrt(lpm(ActiveReturn, 0, 2));
    Sortino2 = (mean(ActiveReturn2) - 0) / sqrt(lpm(ActiveReturn2, 0, 2));

%Upside Ratio (Active @0)
    Upside = lpm(-ActiveReturn, -0, 1) / sqrt(lpm(ActiveReturn, 0, 2));
    Upside2 = lpm(-ActiveReturn2, -0, 1) / sqrt(lpm(ActiveReturn2, 0, 2));

%Create table
    Portfolio1 = [TotalReturn; MonthlyReturn; ActiveRisk; InfoRatio; MaxDD; Omega; Sortino; Upside];
    Portfolio2 = [TotalReturn2; MonthlyReturn2; ActiveRisk2; InfoRatio2; MaxDD2; Omega2; Sortino2; Upside2];
    
    PerformanceMetrics_array = [Portfolio1 Portfolio2];
    
    PerformanceMetrics = array2table(PerformanceMetrics_array,'VariableNames',{'Portfolio1', 'Portfolio2'},'RowNames',{...
        'TotalActiveReturn_%', 'MonthlyActiveReturn_%', 'ActiveRisk', 'InformationRatio', 'MaxDrawdownOfActiveReturn',...
        'OmegaRatio', 'SortinoRatio', 'UpsideRatio'})
%% 2. Create a bar chart for Portfolio 1 that has assets on the horizontal axis and active weights on the vertical, where active weights = portfolio weight - benchmark weight. Repeat for Portfolio 2.


%Compute active weights
    ActiveWeights1 = pwgt(:,10) - BenchmarkWeights;
    ActiveWeights2 = pwgt2(:,10) - BenchmarkWeights;
    
%Create bar chart for Portfolio 1
    fig5 = figure(1);
    c = categorical(Tickers);
    bar(c, ActiveWeights1)
    title('Active Weights for Portfolio 1')

%Create bar chart for Portfolio 2
    fig6 = figure(1);
    c = categorical(Tickers);
    bar(c, ActiveWeights2)
    title('Active Weights for Portfolio 2')

    
%% 3. Discuss your findings. Specifically, interpret each of the performance metrics. Compare and contrast portfolios 1 and 2.

%Portfolios 1 and 2 contained the same number of assets (6), and they
%performed similarly, but Portfolio 2 was better in nearly every aspect.
%For each set of efficient portfolios, I picked the portfolio that had 
%the risk (i.e., standard deviation) closest to the DJ30's, which was 
%0.0354.
% 
%Portfolio 2 had higher total and active portfolio return. These metrics
%were about 8% higher than those of portfolio 1. They measure the excess
%return over the benchmark. Portfolio 2's active risk was also higher, but
%only by about 5%. The reason for this is that portfolio 2 placed more weight
%in assets which produced a higher return but were also riskier. Active risk 
%measures the standard deviation of the active returns. The information
%ratio measures the portfolio returns beyond the returns of a benchmark, 
%in this case the DJ30, compared to the volatility of those returns. A
%higher information ratio is better, and Portfolio 2 had a higher
%information ratio (by about 3%). The maximum drawdown is the largest drop 
%from a peak to a bottom in a certain investment period. I used absolute 
%percentage points (divided by 100 in table) for the maximum drawdown.
%Portfolio 2 had a maximum drawdown of active returns of 7.25 percentage 
%points, while portfolio 1 had a maximum drawdown of 6.84 percentage
%points. In this aspect, portfolio 1 performed slightly better. If you were
%unlucky enough to enter right at the peak and leave right at the bottom of 
%the trough, you would lose less active return in portfolio 1 than in
%portfolio 2. 
% 
%The omega ratio is defined as the probability weighted ratio 
%of gains versus losses for some threshold target (our threshold target for
%this and the following ratios is active return at 0). A larger ratio indicates
%that the asset provides more gains relative to losses and so would be preferred
%by an investor. Portfolio 2 has a larger omega ratio. The Sortino ratio is a 
%modification of the Sharpe ratio but penalizes only those returns falling below
%a threshold target, while the Sharpe ratio penalizes both upside and downside 
%volatility equally. Portfolio 2 has a higher Sortino ratio. The upside ratio is
%a measure of a return of an investment asset relative to the minimal acceptable
%return (0 in our case). Higher is better, and here again portfolio 2
%performs better.
%
%Even though these sets of portfolios had similar efficient frontiers, they
%differed in their amount of upside performance. The efficient frontier
%treats all risk equally, both downside and upside. The omega, Sortino, and
%upside ratios, however, do not penalize an investor for increased upside
%performance. In this respect, portfolio 2 performed solidly better than
%portfolio 1.
%
%Curiously, portfolio 1 and portfolio 2 were invested in the same assets 
%(BA, INTC, MCD, NKE, PG, and V). Where they differed was in their weights.
%I had to set a non-zero minimum weight to be able to specify the minimum
%and maximum number of potentially held assets for portfolio 2. I set this
%minimum at 1%, because that the lowest weight in the DJ30 is 1.11% (with Pfizer).
%The Portfolio object set portfolio 2's holding in INTC and PG equal to
%this minimum bound. This ended up providind portfolio 2 with higher upside
%performance than portfolio 1.
%
%P.S. If you're interested, the six assets both portfolios held were Boeing, Intel,
%McDonald's, Nike, Procter & Gamble, and Visa.