%% Homework 9 - Risk Measurement

%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon descriptiveStats.m.
%Author:
    %Malik Jabati — 8Apr2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Assumptions:
    %Database = BarChart Futures via Quandl (see sakai for API key)
    %Assets = Currency Futures: Australian Dollar Futures, British Pound Futures, Canadian Dollar
    %Futures, and Euro Futures, all with December 2018
    %Frequency = Weekly
    
%% Load in data

%HouseKeeping
    clear all; close all; clc

%Enter the API key
    apikey = 'HxTGtomxL79TZzQg_Ey4'; 

%Establish database connection
    c = quandl(apikey);
    
%Period of interest
    startdate = datetime('12-25-1999','InputFormat','MM-dd-yyyy');
    enddate = datetime('12-31-2018','InputFormat','MM-dd-yyyy');
    
%Identify the database
    source = 'BCCME/';        %Barchart CME U.S. Futures
    periodicity = 'weekly';
    tickers = {'A6Z2018','B6Z2018','D6Z2018','E6Z2018'}; 
    %Australian Dollar Futures, British Pound Futures, Canadian Dollar
    %Futures, and Euro Futures, all with December 2018

%Load in prices data
    for i=1:4
        ticker = tickers(i);
        s = strcat(source,ticker);
      
        data = history(c,s,startdate,enddate,periodicity);  %pull data from Quandl
        ticker_prices = data{:,'Close'};     %get closing price
        futures(:,i) = ticker_prices;
    end
    
    futures = flip(futures);    %flip futures prices
    dates = flip(data.Time);    %grab dates and flip
    
    ret = price2ret(futures);   %compute continuous (log) returns
    dates = dates(2:end);       %grab relevant dates for returns
    
    %Demean the returns
    DM_ret = zeros(size(ret,1),size(ret,2));
    for i=1:size(DM_ret,2)
        DM_ret(:,i) = ret(:,i) - mean(ret(:,i));
    end
    
%% 1. Conduct EDA (exploratory data analysis), including sample statistics, histogram, correlogram of returns, correlogram of squared returns, one white noise test, and a test of ARCH effects.

%Create descriptive statistics table
    
    AUS_Stats = descriptiveStats(ret(:,1));
    AUS_Stats.Properties.VariableNames = {'AUS'};
    
    UK_Stats = descriptiveStats(ret(:,2));
    UK_Stats.Properties.VariableNames = {'UK'};
    
    CAD_Stats = descriptiveStats(ret(:,3));
    CAD_Stats.Properties.VariableNames = {'CAD'};
    
    EUR_Stats = descriptiveStats(ret(:,4));
    EUR_Stats.Properties.VariableNames = {'EUR'};

    SampleStatistics = [AUS_Stats, UK_Stats, CAD_Stats, EUR_Stats]
%Plot histograms of weekly returns on currency futures from 15-Dec-2013 to 23-Dec-2018
%(all avaible data)
    Histograms = figure(1);

    subplot(2,2,1)
    histogram(ret(:,1))
    title('Returns on Australian Dollar Futures')
    
    subplot(2,2,2)
    histogram(ret(:,2))
    title('Returns on British Pound Futures')
    
    subplot(2,2,3)
    histogram(ret(:,3))
    title('Returns on Canadian Dollar Futures')
    
    subplot(2,2,4)
    histogram(ret(:,4))
    title('Returns on Euro Futures')

%Plot correlograms of weekly returns on currency futures from 15-Dec-2013 to 23-Dec-2018
%(all avaible data)
    Correlograms = figure(2);

    subplot(2,2,1)
    autocorr(DM_ret(:,1))
    title('Autocorr. on Australian Dollar Futures')
    
    subplot(2,2,2)
    autocorr(DM_ret(:,2))
    title('Autocorr. on British Pound Futures')
    
    subplot(2,2,3)
    autocorr(DM_ret(:,3))
    title('Autocorr. on Canadian Dollar Futures')
    
    subplot(2,2,4)
    autocorr(DM_ret(:,4))
    title('Autocorr. on Euro Futures')
    
%Plot correlograms of squared weekly returns on currency futures from 15-Dec-2013 to 23-Dec-2018
%(all avaible data)
    SquaredCorrelograms = figure(3);

    subplot(2,2,1)
    autocorr(DM_ret(:,1).^2)
    title('Autocorr. on Aus. Dollar Futures (^2)')
    
    subplot(2,2,2)
    autocorr(DM_ret(:,2).^2)
    title('Autocorr. on British Pound Futures (^2)')
    
    subplot(2,2,3)
    autocorr(DM_ret(:,3).^2)
    title('Autocorr. on Canadian Dollar Futures (^2)')
    
    subplot(2,2,4)
    autocorr(DM_ret(:,4).^2)
    title('Autocorr. on Euro Futures (^2)')
%White noise tests
    
%Use one-sample Kolmogorov-Smirnov test to test for white noise (normal
%distribution) - returns 1 if null hypothesis is rejected

    for i=1:size(DM_ret,2)
        [DataNotNormal(i), pValue_ks(i), stat_ks(i)] = kstest(ret(:,i));
    end
    
%1 means that the test rejects the null hypothesis at the 5% significance level
%for the null hypothesis that the data in vector x comes from a standard normal
%distribution, against the alternative that it does not come from such a distribution
    
    DataNotNormal   %None of the assets have returns from a standard normal distribution (see histogram)
    pValue_ks
    stat_ks
    

%Use Ljung-Box Q-test to test for white noise (residual autocorrelation) -
%returns 1 if null hypothesis is rejected

    for i=1:size(DM_ret,2)
        [Autocorr(i), pValue_lbq(i), stat_lbq(i)] = lbqtest(DM_ret(:,i));
    end
    
%1 means that the test rejects the null hypothesis 
%that there is no autocorrelation in the residual series

    Autocorr   %Only returns on Euro Futures reject the null hypothesis that there is no autocorrelation
                %There is significant autocorrelation for Euro Futures (see
                %correlogram for weekly returns)
    pValue_lbq
    stat_lbq
    
%Test of ARCH effects
%Engle's ARCH test tests for residual heteroscedasticity in the univariate residual series
%- returns 1 if null hypothesis that a series of residuals (rt) exhibits no conditional 
%heteroscedasticity (ARCH effects) is rejected

    for i=1:size(DM_ret,2)
        [ArchEffects(i), pValue_arch(i), stat_arch(i)] = archtest(DM_ret(:,i));
    end

%1 means that the test rejects the null hypothesis 
%that there is no conditional heteroscedasticity (ARCH effects)
    
    ArchEffects   %Only returns on British Pound and Euro Futures reject the null hypothesis that there is no heteroscedasticity
                    %The test concludes there is significant volatility clustering in the residual series
                    %for the British Pound and Euro (see correlogram for
                    %squared weekly returns)
    pValue_arch
    stat_arch
    
%% 2. Set your estimation period as beginning of data to 12/31/2017, and your holding period to the first trading week of 2018. Forecast returns and volatility of returns using an ARMA(1,1)-GARCH(1,1) with Gaussian innovations for your holding period. Estimate optimal portfolio weights using these forecasts.

%15-Dec-2013 to 31-Dec-2017 is rows 1 through 212
%07-Jan-2018 begins at row 213
    
    
    Mdl = arima('ARLags',1,'MALags',1,'Variance',garch(1,1));
    
    NumStepsAhead = 1;  %First trading week of 2018
    
    for i=1:4   %4 assets
        
        %Infer the conditional variances
        EstMdl(i) = estimate(Mdl,DM_ret(1:212,i));
        %CondVarEstimate(:,i) = infer(EstMdl(i),DM_ret(1:212,i));
        
        %Simulate a path from the model so that can forecast
        %[simCondVariance(:,i),simDM_ret(:,i)]=simulate(EstMdl(i),length(DM_ret(1:212,i))); %Simulate the de-meaned returns (i.e. epsilon=innovation process) 
        
        %Forecast from the GARCH model
        [ReturnForecast(1,i), ~, VarianceForecast(1,i)] = forecast(EstMdl(i),NumStepsAhead,'Y0',DM_ret(212,i));
    end
    %Compute covariance matrix for holding period
    Covar = cov(DM_ret(1:112,:));
    %Place forecasted variances into covariance matrix for each asset
    Covar(1,1) = VarianceForecast(1,1);
    Covar(2,2) = VarianceForecast(1,2);
    Covar(3,3) = VarianceForecast(1,3);
    Covar(4,4) = VarianceForecast(1,4);
    
    
    p = Portfolio;
    p = setAssetMoments(p,ReturnForecast(1,:),Covar);
    p = setDefaultConstraints(p); %Set the portfolio constraints
    p = setBounds(p, -1, 1,'BoundType', 'Simple');
    OptimalWeights(1,:) = estimateMaxSharpeRatio(p)';
    
%% 3. Roll your estimation and holding periods forward by one week. Repeat the forecast and optimization procedure.

%22-Dec-2013 to 07-Jan-2017 is rows 2 through 213
%14-Jan-2018 begins at row 214
    
    
%Mdl and NumStepsAhead declared above
    
    for i=1:4   %4 assets
        
        %Infer the conditional variances
        EstMdl(i) = estimate(Mdl,DM_ret(1+1:212+1,i));
        
        %Forecast from the GARCH model
        [ReturnForecast(1+1,i), ~, VarianceForecast(1+1,i)] = forecast(EstMdl(i),NumStepsAhead,'Y0',DM_ret(212+1,i));
    end
    %Compute covariance matrix for holding period
    Covar = cov(DM_ret(1+1:112+1,:));
    %Place forecasted variances into covariance matrix for each asset
    Covar(1,1) = VarianceForecast(1+1,1);
    Covar(2,2) = VarianceForecast(1+1,2);
    Covar(3,3) = VarianceForecast(1+1,3);
    Covar(4,4) = VarianceForecast(1+1,4);
    
    
    p = Portfolio;
    p = setAssetMoments(p,ReturnForecast(1+1,:),Covar);
    p = setDefaultConstraints(p); %Set the portfolio constraints
    p = setBounds(p, -1, 1,'BoundType', 'Simple');
    OptimalWeights(1+1,:) = estimateMaxSharpeRatio(p)';
    
%% 4. Repeat this rolling window until the end of 2018.

    
    %Mdl and NumStepsAhead declared above
    for j=0:50  %51 total periods
        for i=1:4   %4 assets       
            %Infer the conditional variances
            EstMdl(i) = estimate(Mdl,DM_ret(1+j:212+j,i));
            
            %Forecast from the GARCH model
            [ReturnForecast(1+j,i), ~, VarianceForecast(1+j,i)] = forecast(EstMdl(i),NumStepsAhead,'Y0',DM_ret(212+j,i));
        end
        %Compute covariance matrix for holding period
        Covar = cov(DM_ret(1+j:112+j,:));
        %Place forecasted variances into covariance matrix for each asset
        Covar(1,1) = VarianceForecast(1+j,1);
        Covar(2,2) = VarianceForecast(1+j,2);
        Covar(3,3) = VarianceForecast(1+j,3);
        Covar(4,4) = VarianceForecast(1+j,4);
        
        
        p = Portfolio;
        p = setAssetMoments(p,ReturnForecast(1+j,:),Covar);
        p = setDefaultConstraints(p); %Set the portfolio constraints
        p = setBounds(p, -1, 1,'BoundType', 'Simple');
        OptimalWeights(1+j,:) = estimateMaxSharpeRatio(p)';
    end
    
%% 5. Track the performance of your 1 week ahead rolling portfolio by i) creating a performance graphic of a notional $100 investment in this portfolio, ii) descriptive statistics, iii) VaR using "normal” method, iv) VaR using “historical” method. Provide 1-2 paragraph discussion of your results.

    
% i) create a performance graphic of a notional $100 investment in this portfolio
    PortfolioReturns = OptimalWeights .* ret(213:end,:);    %Compute weighted returns for each asset      
    PortfolioReturns = sum(PortfolioReturns,2);             %Sum weighted returns for each week
    PortfolioReturnsCum = cumsum([1; PortfolioReturns])*100;   %Add starting value of 100 for initial investment and cumulatively sum
  
    %Performance of 1-week ahead rolling investment
    PortfolioPerformance = figure(4);
    plot(dates(212:end),PortfolioReturnsCum)
    title('Portfolio Performance of $100 Investment')
    ylabel('Dollars ($)')
    grid on
    
% ii) descriptive statistics
    Portfolio_Stats = descriptiveStats(PortfolioReturns)
    
% iii) VaR using "normal" method (5%)
    Zscore   = norminv(0.05);
    Sigma = std(PortfolioReturns);
    VaR_normal = -Zscore*Sigma

% iv) VaR using "historical" method (5%)
    VaR_historical = -quantile(PortfolioReturns,0.05)
    
% Provide 1-2 paragraph discussion of your results.

%My rolling-window ARMA(1,1)-GARCH(1,1) model did not perform well in forecasting future
%returns during the holding period. My model saw total losses of ~15% at the end
%of the 51st trading week in 2018. An equally weighted portfolio without
%rebalancing would have seen total losses of ~8% by the end of 2018. Each
%asset performed negatively over the trading period, and I set my portfolio
%constraints to allow for a lower bound weight of -1 and upper bound weight
%of 1, with total portfolio weights summing to 1 (fully invested). I did this to minimize
%risks from leveraging. When determing the asset weights for each week, I
%used the forecasted return and variance of each asset. I used the rolling
%estimation period to determine covariances between assets. I maximized the
%portfolio's Sharpe ratio using these values to determine the portfolio's
%asset weights.

%My portfolio had negative mean daily returns, which were nearly double
%that of an equally weighted portfolio. The variance, a measure of 
%volatility, of my portfolio's returns was roughly equal to the average variance of the assets.
%The high volatility indicates that my portfolio did a poor job of minimizing risk,
%having the opposite effect expected from diversification. A normal
%distribution has a kurtosis of 3. My portfolio had a kurtosis of 2.74,
%which means that I had a distribution with slightly fewer outliers (or
%extreme values) than would be expected in a normal distribution. My
%skewness was slightly negative, meaning that the weekly returns spread
%out more to the left of the mean than to the right. The range between the
%highest and lowest weekly return was about 5.3%, the minimum was -3.2%,
%and the maximum was 2.1%. The VaR using both the "normal" and "historical"
%methods came out to about 2.0%, meaning that 5% of the time I could expect
%to lose at least 1/50th of my portfolio value in a week. The four assets in
%my portfolio were currency futures. Currency futures are futures contracts 
%for currencies that specify the price of exchanging one currency for another
%at a future date. The rate for currency futures contracts is derived from 
%spot rates of the currency pair. Currency futures are used to hedge the risk 
%of receiving payments in a foreign currency. If the spot rate of a currency 
%pair decreases, the futures prices have a high probability of decreasing.
%The spot rate in this case is the dollar-denominated exchange rate of Australian dollars,
%British pounds, New Zealand dollars, and Euros for U.S. dollars. The
%U.S. dollar had a strong run in 2018 (strengthening more than investors 
%expected), and the four foreign currencies weakened relative to it. This 
%caused the stable decline in prices of the December 2018 futures for the 
%four foreign currencies.