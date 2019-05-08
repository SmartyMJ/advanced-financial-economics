%%COMMENTS:
%good work, really detailed discussion
%great job mentioning all of the possible causes weights didn't move much


%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon PriceHistoryGSCI.xlsx.
%Author:
    %Malik Jabati — 31Mar2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Assumptions:
    %Use daily returns data for the assets of the S&P Dow Jones GSCI for as long a period as available, 
    %while maintaining a balanced panel of data. The data is available on Factset.
    %Use a Black-Litterman allocation strategy with a Normal reference model.
    %Benchmark: S&P Dow Jones GSCI.

%% Load in data

%All data was cleaned in Excel before importing

%HouseKeeping
    clear all; close all; clc
    
%Import data
    [num, txt, raw] = xlsread('PriceHistoryGSCI.xlsx','ASSETS');
    
    Tickers = txt(2:end);
    
%Set dates    
    Dates = x2mdate(num(:,1));
    Dates = datetime(Dates,'ConvertFrom','datenum');
    
%Set asset names
    AssetNames = {'Wheat', 'Kansas Wheat', 'Corn', 'Soybeans', 'Cotton', 'Sugar', 'Coffee', 'Cocoa',... %Agriculture
        'Crude Oil', 'Brent Crude', 'Unleaded Gasoline', 'Heating Oil', 'GasOil', 'Natural Gas',...     %Energy
        'Aluminum', 'Copper', 'Lead', 'Nickel', 'Zinc',...                                              %Industrial Metals
        'Feeder Cattle', 'Live Cattle', 'Lean Hogs',...                                                 %Livestock
        'Gold', 'Silver'};                                                                              %Precious Metals
    
%Set asset weights
    AssetWeights = [0.0277, 0.0115, 0.0436, 0.0314, 0.0141, 0.0154, 0.0072, 0.0032,...  %Agriculture
        0.2642, 0.1861, 0.0448, 0.0445, 0.0556, 0.0311,...                              %Energy
        0.0389, 0.0445, 0.0078, 0.0076, 0.0128,...                                      %Industrial Metals
        0.0127, 0.0348, 0.0191,...                                                      %Livestock
        0.0372, 0.0042];                                                                %Precious Metals
    
%Set sector names
    SectorNames = {'Agriculture', 'Energy', 'Industrial Metals', 'Livestock', 'Precious Metals'};
    
%Set sector weights
    SectorWeights = [0.1541, 0.6263, 0.1116, 0.06653, 0.04144];
    
%Compute returns
    temp = num(:,3:end); 
    N = size(temp,2); %# of assets
    T = size(temp,1); % # of months
    Returns = zeros(T-1,N);
    
%Create Returns matrix
    for i =1:N
        Returns(1:T-1,i) = temp(2:end,i)./temp(1:end-1,i)-1; 
    end
    
%Form expectations on historical avgs; ANNUALIZE them too
    ExpRet = (1+mean(Returns)).^252 - 1; 
    ExpCov = 252*cov(Returns); %Assume that daily returns are independent from each other
    
%% Portfolio 1

%BL settings
    delta1 = 3; %risk aversion calibrated outside model.
    weq = AssetWeights; %Equilibrium weights in market portfolio.
    sigma = ExpCov; %prior covariance matrix
    tau = 1/(T-N); %uncertainty in prior estimate of mean returns; common calibration. 
    P = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
         0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];%"Pick" matrix
    
    Q = [ExpRet(23);
         ExpRet(24);
         .01]; %Set the views
    Omega1 = diag(diag(P*tau*sigma*P')); %Variances of views; Common to scale with prior variances

%Optimize via BL    
    [PosteriorReturn1, PosteriorCov1, PosteriorW1, PosteriorPW1] = hlblacklitterman(delta1, weq, sigma, tau, P, Q, Omega1); 
    PriorReturn = ExpRet';
    PriorW = weq'; 

%Update Weights array
    Weights = [AssetWeights' PosteriorPW1];
  
%% Portfolio 2

%BL settings
    delta2 = 3; %risk aversion calibrated outside model.
    weq = AssetWeights; %Equilibrium weights in market portfolio.
    sigma = ExpCov; %prior covariance matrix
    tau = 1/(T-N); %uncertainty in prior estimate of mean returns; common calibration. 
    P = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
         0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];%"Pick" matrix
    
    Q = [ExpRet(23);
         ExpRet(24);
         .01]; %Set the views
    Omega2 = diag(diag(P*tau*sigma*P')); %Variances of views; Common to scale with prior variances
    Omega2(3,3) = Omega2(1,1)+Omega2(2,2); %Repeat Portfolio 1, but set the uncertainty of the second view to be twice that of the first view.
                                        %Equal to 2*((Omega(1,1)+Omega(2,2))/2))

%Optimize via BL    
    [PosteriorReturn2, PosteriorCov2, PosteriorW2, PosteriorPW2] = hlblacklitterman(delta2, weq, sigma, tau, P, Q, Omega2); 
    PriorReturn = ExpRet';
    PriorW = weq'; 

%Update Weights array
    Weights = [Weights PosteriorPW2];
    
%% Portfolio 3

%BL settings
    delta3 = 2*delta1; %risk aversion calibrated outside model; Repeat Portfolio 1, but double risk aversion.
    weq = AssetWeights; %Equilibrium weights in market portfolio.
    sigma = ExpCov; %prior covariance matrix
    tau = 1/(T-N); %uncertainty in prior estimate of mean returns; common calibration. 
    P = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
         0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];%"Pick" matrix
    
    Q = [ExpRet(23);
         ExpRet(24);
         .01]; %Set the views
    Omega3 = diag(diag(P*tau*sigma*P')); %Variances of views; Common to scale with prior variances

%Optimize via BL    
    [PosteriorReturn3, PosteriorCov3, PosteriorW3, PosteriorPW3] = hlblacklitterman(delta3, weq, sigma, tau, P, Q, Omega3); 
    PriorReturn = ExpRet';
    PriorW = weq'; 

%Update Weights array
    Weights = [Weights PosteriorPW3];
    
%% Portfolio 4

%BL settings
    delta4 = 3; %risk aversion calibrated outside model
    weq = AssetWeights; %Equilibrium weights in market portfolio.
    sigma = ExpCov; %prior covariance matrix
    tau = 1/(T-N); %uncertainty in prior estimate of mean returns; common calibration. 
    P = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;...
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
         0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];%"Pick" matrix
    
    Q = [ExpRet(23);
         ExpRet(24);
         .01]; %Set the views
    Omega4 = diag(diag(P*tau*sigma*P')); %Variances of views; Common to scale with prior variances
    Omega4(1,1) = .00000001;    %Repeat Portfolio 1, but with virtually no uncertainty in either view.
    Omega4(2,2) = .00000001;
    Omega4(3,3) = .00000001;
    
%Optimize via BL    
    [PosteriorReturn4, PosteriorCov4, PosteriorW4, PosteriorPW4] = hlblacklitterman(delta4, weq, sigma, tau, P, Q, Omega4); 
    PriorReturn = ExpRet';
    PriorW = weq'; 

%Update Weights array
    Weights = [Weights PosteriorPW4];
    
%% 1. Fill in the following table


    WeightsTable = array2table(Weights,'VariableNames',{'Benchmark', 'Port_1', 'Port_2', 'Port_3', 'Port_4'});
    IndexTable = array2table(AssetNames','VariableNames',{'Index_Name'});
    
    WeightsTable = [IndexTable WeightsTable];
    
    WeightsTable.Properties.RowNames = Tickers(2:end)';
    
    WeightsTable
%% 2. Within 1-3 paragraphs, please compare and contrast each of the portfolios. Make sure to largely focus on the differences in weights between portfolios and the reasons for these differences.

ChangeTable = WeightsTable;

ChangeTable.Port_1 = ((ChangeTable.Port_1 ./ ChangeTable.Benchmark) - 1)*100;
ChangeTable.Port_2 = ((ChangeTable.Port_2 ./ ChangeTable.Benchmark) - 1)*100;
ChangeTable.Port_3 = ((ChangeTable.Port_3 ./ ChangeTable.Benchmark) - 1)*100;
ChangeTable.Port_4 = ((ChangeTable.Port_4 ./ ChangeTable.Benchmark) - 1)*100;

ChangeTable

SectorTable = [array2table(SectorNames',"VariableNames",{'Sector_Name'}) array2table(SectorWeights',"VariableNames",{'Sector_Weights'})]
mean(ExpRet)

%In all 4 portfolios, the weights decreased by 0.0338% for every asset
%except for the four assets on which we had views -- WTI Crude Oil, Brent
%Crude, and the Precious Metals (i.e., Gold and Silver). In Portfolio 1, 
%View 1a was an absolute view that set Gold returns equal to their historical
%average, View 1b was an absolute view that set Silver returns equal to
%their historical average, and View 2 was a relative view that set Brent
%Crude outperforming WTI Crude Oil by 1 percentage point. Uncertainty was
%set to the variance of assets as per equation 39 of Walters (2014), with
%tau equal to 1/(T-k) as the best quadratic unbiased estimator. Portfolio 2
%was a repeat of Portfolio 1, but set the uncertainty of the second view to be 
%twice that of the first view. Portfolio 3 was a repeat of Portfolio 1, but with
%double the risk aversion. Portfolio 1 was a repeat of Portfolio 1, but with 
%virtually no uncertainty in either view, which I set to 0.00000001.
%Assets in the Agriculture sector made up ~15% of the portfolios, Energy
%assets made up ~63%, Industrial Metals assets made up ~11%, Livestock
%assets made up ~7%, and Precious Metals assets made up ~4%. These sector
%weights stayed roughly constant between the portfolios and benchmark.

%For all portfolios, the change in asset weights was relatively small.
%Looking at percentage change is potentially misleading for intra-portfolio
%analysis, but is useful in comparing differences between portfolios.
%The biggest percentage change in weight occured in Portfolio 4, where
%Silver's weight increased by 2.7%. The smallest percentage change in
%weight occured in Portfolios 1 and 3, where WTI Crude Oil's weight
%decreased by only 0.0072%. The specification of the variance, or 
%uncertainty, of the views in Portfolio 1 essentially equally weights
%the investor's views and the market equilibrium weights. By including
%tau in the expression, the posterior estimate of the returns becomes
%independent of tau as well. This seems to be the most common method 
%used in the literature (Walters 2014, p. 14). In Portfolio 1, the weight
%for Silver increased by 0.94%, Gold's weight increased by 0.056%, Brent's
%weight decreased by 0.072%, and WTI's weight decreased by .0072%. Our analysis
%has only partial views, that is views on a subset of the assets, so by using
%a posterior estimate of the variance we will tilt the posterior weights towards
%assets with lower variance (higher precision of the estimated mean) and away from
%assets with higher variance (lower precision of the estimated mean). Thus the 
%existence of the views and the updated covariance will tilt the optimizer towards
%using or not using those assets. This tilt will not be very large if we are working
%with a small value of tau, but it will be measurable. Our views on Brent,
%WTI, and precious metals tilted our optimizer towards and away from those
%assets. In Portfolio 2, the weight for Silver increased by 0.96%, Gold's
%weight increased by 0.058%, Brent's weight decreased by 0.040%, and WTI's
%weight decreased by .030%. The weight of the precious metals increased by
%more in Portfolio 2 than in Portfolio 1, and the weights of Brent and WTI
%decreased by more in Portfolio 1 than in Portfolio 2. This is because the
%uncertainty for View 2 was greater in Portfolio 2 than in Portfolio 1.
%That tilted our optimizer towards the precious metals and away from WTI
%and Brent. In Portfolio 3, the weight for Silver increased by 0.94%, Gold's
%weight increased by 0.056%, Brent's weight decreased by 0.072%, and WTI's
%weight decreased by .0072%. The percentage change in weights for Portfolio
%3 is equal to that of Portfolio 1. The doubling of risk aversion then
%affected every view equally and did not affect our final weighing of
%assets. Traditionally, depending on their risk aversion an investor will 
%hold arbitrary fractions of their wealth in the risk free asset and/or 
%the CAPM market portfolio. All of our assets were risky assets, so the
%risk aversion level did not have any impact (or a non-applicable impact)
%on our holdings in risk free assets. In Portfolio 4, the weight for Silver
%increased by 2.7%, Gold's weight decreased by 0.14%, Brent's weight
%decreased by 0.10%, and WTI's weight increased by .016%. Portfolio had
%virtually no uncertainty in either view. I also set the uncertainty equal
%amongst all views. This affected the tilt of the optimizer and changed our
%weightings from their values in Portfolio 1.
% 
%There are three additional considerations when constructing our Black-Litterman model
%that affects our output portfolio weights.
%[1] The tau used in our calculations was extremely small at 0.00039 (about two orders
%of magnitude smaller than the 0.025 to 0.050 recommended by Black and
%Litterman). Since we often build the known covariance matrix of returns, 
%sigma, from historical data we can use methods from basic statistics to compute
%tau, as tau*sigma is analogous to the standard error. We can also estimate 
%tau based on our confidence in the prior distribution. Note that both of 
%these techniques provide some intuition for selecting a value of tau which
%is closer to 0 than to 1. Black and Litterman (1992), He and Litterman
%(1999) and Idzorek (2005) all indicate that in their calculations
%they used small values of tau, on the order of 0.025 – 0.050
%(Walters 2014, p. 20). Our small tau affected the optimization of our
%Black-Litterman model and decreased the difference in our change in weights because tau
%is a constant of proportionality. A slightly larger tau would likely make for a better model.
%[2] Weights within the portfolios are based on relative views, but there was high
%correlation between the assets because they are all commodities. We would prefer
%less correlation between the assets in our model. Some views may actually be pulling
%the posterior towards the prior, and the investor could strengthen these views, or 
%weaken views which pull the posterior away from the prior. This point may seem 
%non-intuitive. Given that the views are indirectly coupled by the covariance matrix,
%one would expect that the views only push the posterior distribution away from the prior.
%However, because the views can be conflicting, either directly or via the correlations, 
%any individual view can have a net impact pushing the posterior closer to the prior,
%or pushing it further away.
%[3] Variances in omega are unbounded and are low. This affects our confidence and certainty
%of the views. Idzorek's extension would help here, because it provides a method for specifying
%the confidence in the view in terms of a percentage move of the weights on the interval
%from 0% confidence to 100% confidence.