%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon tickerData_clean.xlsx and MacroFactors.xlsx files.
    %Interpretations are included in Malik_Jabati_HW1_interpretations.txt.
    %I pulled price data from WRDS and removed tickers with any NaN values.
%Author:
    %Malik Jabati — 24Jan2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
    
%Housekeeping
    clear all; close all; clc;
    
%Open the data file
    [data,label,raw] = xlsread('tickerData_clean.xlsx');
    [data_f,label_f,raw_f] = xlsread('MacroFactors.xlsx');
    
%Cleaning data
    Tickers_temp = cell2table(label(2:end,3));
    Tickers_temp.Properties.VariableNames = {'Tickers'};
    
    PricesTable_temp = cell2table(raw(2:end,[2 3 5])); %Convert prices array to table    
    PricesTable_temp.Properties.VariableNames = {'Date', 'Ticker', 'Price'};
    
    PricesTable = unstack(PricesTable_temp,'Price','Ticker'); %Variables are ordered by Permno from lowest to highest
    
    PricesArray = table2array(PricesTable);
    
    VarNames = PricesTable.Properties.VariableNames; %List of variables
    
    NaNCols = any(isnan(PricesArray));
    PricesArray = PricesArray(:,~NaNCols); %Clear out cols with ANY NaNs
    
    VarNames = VarNames(:,~NaNCols); %Clear out bad variable names
    
    %Recreate cleaned table
    PricesTable = array2table(PricesArray,'VariableNames',VarNames);
    
    Dates_temp = PricesArray(2:end,1);
    
    %Line up factors
    Factors = data_f(1:end-1,2:end);
    
    %Generating returns
    Returns = real(diff(log(PricesArray(:,2:end))));
    
    %Convert dates format
    year = floor(Dates_temp/10000);
    month = floor((Dates_temp-year*10000)/100);
    day = Dates_temp - (year*10000 + month*100);
    Dates = datenum(year,month,day);

% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%First Pass: Estimate the betas
    beta = zeros(size(Returns,2),4); 
    for i = 1:size(Returns,2) %loop over each asset
        y = Returns(:,i); %Set the "y"
        X = [ones(size(Factors,1),1) Factors]; %set the "X"; Add ones vector for intercept. 
        [b,bint,r,rint,stats] = regress(y,X); %Conduct regression. 
        beta(i,1:end) = b(2:end,1); %grab the betas and stack, b is a px1 vector where p is the # of predictors, so it's only each factor and NFP
    end
    
    AvgReturns = mean(Returns);
    OneA = array2table([AvgReturns' beta], 'VariableNames', {'AvgReturn', 'B_10yr', 'B_Spread', 'B_SP500Yields', 'B_IndProd'});
    OneA.Properties.RowNames = VarNames(1,2:end)';
    
    %OneB: See attached doc for analysis
    OneB_var = var(OneA{:,{'B_10yr', 'B_Spread', 'B_SP500Yields', 'B_IndProd'}});
    OneB_corr = corr(OneA{:,{'B_10yr', 'B_Spread', 'B_SP500Yields', 'B_IndProd'}});
    summary(OneA)
    
%Second Pass: Regress avg returns on estimated betas (assume no risk free)
    clear b y X
    y = AvgReturns'; %Set the y 
    X = [beta]; % Set the X. Note that the fitlm function includes intercept by default. 
    mdl = fitlm(X,y); %this is another way to run regressions in Matlab.  Nice formatting for display. 
    lambda = mdl.Coefficients(:,[1 3]);
    
    OneC = lambda;
    OneC.Properties.VariableNames = {'Lambdas', 'T_Stat'};
    OneC.Properties.RowNames = {'Intercept';'L_10yr';'L_Spread';'L_SP500Yields';'L_IndProd'};
    
    %OneD: See attached doc for analysis
    OneD_stats = mdl.Coefficients
    
% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%First Pass: Estimate the betas
    beta = zeros(size(Returns,2),4); 
    for i = 1:size(Returns,2) %loop over each asset
        y = Returns(:,i); %Set the "y"
        X = [ones(size(Factors,1),1) Factors]; %set the "X"; Add ones vector for intercept. 
        [b,bint,r,rint,stats] = regress(y,X); %Conduct regression. 
        beta(i,1:end) = b(2:end,1); %grab the betas and stack, b is a px1 vector where p is the # of predictors, so it's only each factor and NFP
    end
    
    AvgReturns = mean(Returns);
    TwoA = array2table([AvgReturns' beta], 'VariableNames', {'AvgReturn', 'B_10yr', 'B_Spread', 'B_SP500Yields', 'B_IndProd'});
    TwoA.Properties.RowNames = VarNames(1,2:end)';
    
    %TwoB: See attached doc for analysis
    TwoB_var = var(TwoA{:,{'B_10yr', 'B_Spread', 'B_SP500Yields', 'B_IndProd'}});
    TwoB_corr = corr(TwoA{:,{'B_10yr', 'B_Spread', 'B_SP500Yields', 'B_IndProd'}});
    summary(TwoA)

%Second Pass: Regress returns on estimated betas for each t, then estimate lambdas (assume no risk free)
    T = size(Returns,1); %# of months, since the data is given by month
    for t = 1:T %Loop over each month
    y = Returns(t,:); %Set the y 
    X = [ones(size(beta,1),1) beta]; % Set the X; add ones for intercept, beta now but still loop over T's
    [lambda,bint,r,rint,stats] = regress(y',X); %Regress returns on estimated betas
    %Grab the Lambda for that month T
        Lambda1(t,:) = lambda(:,1)'; 
    end
    
%TwoC: Create five time series graphs for lambda
    
    Dates_plot = datetime(Dates,'ConvertFrom','datenum');
    
    %L0
    plot(Dates_plot,Lambda1(:,1))
    title('Factor premium of intercept over time')
    xlabel('Date') 
    ylabel('L0 premium')

    %L1
    figure
    plot(Dates_plot,Lambda1(:,2))
    title('Factor premium of US 10yr Treasury yield over time')
    xlabel('Date') 
    ylabel('L1 premium')
    
    %L2
    figure
    plot(Dates_plot,Lambda1(:,3))
    title('Factor premium of 10yr-2yr US Treasury yield spread over time')
    xlabel('Date') 
    ylabel('L2 premium')
    
    %L3
    figure
    plot(Dates_plot,Lambda1(:,4))
    title('Factor premium of SP500 dividend yield over time')
    xlabel('Date') 
    ylabel('L3 premium')
    
    %L4
    figure
    plot(Dates_plot,Lambda1(:,5))
    title('Factor premium of % change in industrial production over time')
    xlabel('Date') 
    ylabel('L4 premium')

    %TwoD: See attached doc for analysis

    %Summarize by avg the lambdas over time
    AvgL = mean(Lambda1)';
    StdL = std(Lambda1)';
    %Create a t-statistic for this avg lambda
    tstat = AvgL ./ ((1/sqrt(T)).*StdL); 
    
    %Clean up for readability
    out = [AvgL tstat]; %Group for the table
    TwoE = array2table(out,'VariableNames',{'Lambdas', 'T_Stat'},'RowNames',{'Intercept';'L_10yr';'L_Spread';'L_SP500Yields';'L_IndProd'}); %create a table
    
    %TwoF: See attached doc for analysis
    
    