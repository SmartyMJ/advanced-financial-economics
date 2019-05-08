%%COMMENTS
%did you clean your data here at all? there should have been some negative prices?
%Good job with your interpretations, but definitely need more
%you are lacking detail, though initially you did try to add theory
%you spend a lot of time comparing the sample stats but do not give a reason for it
%you need to be more specific-- why would one factor out perform the other
%what is the specific cause? economies of scale? cost control? etc....








%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon tutorial.db and descriptiveStats.m.
%Author:
    %Malik Jabati — 31Jan2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Table Names:
    %Ticker, Year, prices, NumShare, BookV
%Assumptions:
    %I use simple returns, assuming each portfolio is value weighted.
    %I create 8 sets of the 6 Fama-French portfolios for the years 2010-2017.
    %Big portfolios are defined as those with size bigger than than median value.
    %Value portfolios are those with BE/ME greater than the 70th percentile.
    %Growth portfolios are those with growth less than or equal to the 30th percentile.
    %The 6 Fama-French portfolios are constructed from the intersection of the 2 size portfolios and 3 BE/ME portfolios.
    
%HouseKeeping
    clear all; close all; clc
    
%Create a read-only SQLite connection to database
    conn = sqlite('tutorial.db','readonly');
      
%Grab and standardize fields from database
    sql1 = 'Select * from Year'; %Set up the sql statement
    Year = fetch(conn,sql1); %fetch the data
    Year = double(cell2mat(Year)); %convert to array

    sql2 = 'Select * from Ticker'; %Set up the sql statement
    Ticker = fetch(conn,sql2); %fetch the data
    Ticker = double(cell2mat(Ticker)); %convert to array
    
    sql3 = 'Select * from prices'; %Set up the sql statement
    Prices = fetch(conn,sql3); %fetch the data
    Prices = cell2mat(Prices); %convert to array
    
    sql4 = 'Select * from NumShare'; %Set up the sql statement
    NumShare = fetch(conn,sql4); %fetch the data
    NumShare = cell2mat(NumShare); %convert to array
    
    sql5 = 'Select * from BookV'; %Set up the sql statement
    BookV = fetch(conn,sql5); %fetch the data
    BookV = cell2mat(BookV); %convert to array
    
%Calculate size and BE/ME
    Size = Prices .* NumShare;
    
    BookValue = BookV .* NumShare;
    
    BookToMarket = BookValue ./ Size;
      
%Create Big and Small matrices
    Big = zeros(size(Size,1),size(Size,2));
    Small = zeros(size(Size,1),size(Size,2));

    %Give Big and Small matrices dummy variables
    for i = 1:size(Size,2) %loop over each year
        for j = 1:size(Size,1) %loop over each asset
            if Size(j,i) > median(Size(:,i)) %if current asset's size is greater than median size for year
                Big(j,i) = 1; %mark 1 in Big table
            else              %if current asset's size is less than or equal to median size for year
                Small(j,i) = 1; %mark 1 in Small table
            end
        end
    end
    
%Create Value, Neutral, and Growth matrices
    Value = zeros(size(BookToMarket,1),size(BookToMarket,2));
    Neutral = zeros(size(BookToMarket,1),size(BookToMarket,2));
    Growth = zeros(size(BookToMarket,1),size(BookToMarket,2));
    
    %Give Value, Neutral, and Growth matrices dummy variables
    for i = 1:size(BookToMarket,2) %loop over each year
        for j = 1:size(BookToMarket,1) %loop over each asset
            if BookToMarket(j,i) > prctile(BookToMarket(:,i),70) %if current asset's BE/ME is greater than 70th percentile BE/ME for year
                Value(j,i) = 1; %mark 1 in Value table
            elseif BookToMarket(j,i) <= prctile(BookToMarket(:,i),30) %if current asset's BE/ME is less than or equal to 30th percentile BE/ME for year
                Growth(j,i) = 1; %mark 1 in Growth table
            else %if between 30th and 70th percentile BE/ME
                Neutral(j,i) = 1;
            end
        end
    end  
    
%Create market equity matrices for 6 Fama-French portfolios

    SmallValueEquity = Size .* Small .* Value;
    SmallNeutralEquity = Size .* Small .* Neutral;
    SmallGrowthEquity = Size .* Small .* Growth;
    
    BigValueEquity = Size .* Big .* Value;
    BigNeutralEquity = Size .* Big .* Neutral;
    BigGrowthEquity = Size .* Big .* Growth;
    
%Create portfolio price matrices for 6 Fama-French portfolios

    SmallValuePrices = sum(SmallValueEquity .* Prices) ./ sum(SmallValueEquity);
    SmallNeutralPrices = sum(SmallNeutralEquity .* Prices) ./ sum(SmallNeutralEquity);
    SmallGrowthPrices = sum(SmallGrowthEquity .* Prices) ./ sum(SmallGrowthEquity);
    
    BigValuePrices = sum(BigValueEquity .* Prices) ./ sum(BigValueEquity);
    BigNeutralPrices = sum(BigNeutralEquity .* Prices) ./ sum(BigNeutralEquity);
    BigGrowthPrices = sum(BigGrowthEquity .* Prices) ./ sum(BigGrowthEquity);
    
%Create portfolio simple returns
    SmallValueReturns = tick2ret(SmallValuePrices');
    SmallNeutralReturns = tick2ret(SmallNeutralPrices');
    SmallGrowthReturns = tick2ret(SmallGrowthPrices');
    
    BigValueReturns = tick2ret(BigValuePrices');
    BigNeutralReturns = tick2ret(BigNeutralPrices');
    BigGrowthReturns = tick2ret(BigGrowthPrices');
    
%Create SMB and HML factors
    SMB = (1/3)*(SmallValueReturns + SmallNeutralReturns + SmallGrowthReturns) - (1/3)*(BigValueReturns + BigNeutralReturns + BigGrowthReturns);
    HML = (1/2)*(SmallValueReturns + BigValueReturns) - (1/2)*(SmallGrowthReturns+BigGrowthReturns);
    
%% A: Create a time series chart that overlays the Small Value, Small Neutral, and Small Growth portfolios.
    figure('Name','Small Portfolios')

    SV_ts = timeseries(SmallValueReturns(12:end),1:8);
    SV_ts.TimeInfo.Units = 'years';
    SV_ts.TimeInfo.StartDate = '2010';     % Set start date.
    SV_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    SV_ts.Time = SV_ts.Time - SV_ts.Time(1);        % Express time relative to the start date.

    plot(SV_ts)
    grid on
    hold on

    SN_ts = timeseries(SmallNeutralReturns(12:end),1:8);
    SN_ts.TimeInfo.Units = 'years';
    SN_ts.TimeInfo.StartDate = '2010';     % Set start date.
    SN_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    SN_ts.Time = SN_ts.Time - SN_ts.Time(1);        % Express time relative to the start date.
    
    plot(SN_ts)
    
    SG_ts = timeseries(SmallGrowthReturns(12:end),1:8);
    SG_ts.TimeInfo.Units = 'years';
    SG_ts.TimeInfo.StartDate = '2010';     % Set start date.
    SG_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    SG_ts.Time = SG_ts.Time - SG_ts.Time(1);        % Express time relative to the start date.
    
    plot(SG_ts)
    hold off
    
    title('Time Series: Small Portfolios')
    xlabel('Year')
    ylabel('Returns')
    legend('Small Value','Small Neutral','Small Growth','Location','northwest')

%% B: Provide descriptive statistics of the Small Value, Small Neutral, and Small Growth portfolios. Present the descriptive statistics in a table. Interpret your results 1-2 paragraphs.
    % i.e., mean, std error, median, mode, std. deviation, sample variance, kurtosis, skewness, range, min, max, sum, count
    StatRows = {'Mean', 'StdError', 'Median', 'Mode', 'StdDeviation', 'SampleVariance', 'Kurtosis', 'Skewness', 'Range', 'Min', 'Max', 'Sum', 'Count'};
    SmallVars = {'SmallValue', 'SmallNeutral','SmallGrowth'};
    
    SV_stats = descriptiveStats(SmallValueReturns(12:end));
    SN_stats = descriptiveStats(SmallNeutralReturns(12:end));
    SG_stats = descriptiveStats(SmallGrowthReturns(12:end));
    
    SmallStats = array2table([SV_stats SN_stats SG_stats],'RowNames',StatRows,'VariableNames',SmallVars)

%                      SmallValue    SmallNeutral    SmallGrowth
%                      __________    ____________    ___________
%
%    Mean               0.0011038     -0.006516      -0.00070752
%    StdError             0.01373      0.012688         0.011054
%    Median            -0.0035534      -0.01678        0.0016056
%    Mode                -0.05052      -0.04895        -0.039848
%    StdDeviation        0.038835      0.035886         0.031265
%    SampleVariance     0.0015081     0.0012878        0.0009775
%    Kurtosis               3.155        1.8997           2.1047
%    Skewness             0.75632       0.55119          0.31806
%    Range                0.12862      0.098717         0.092612
%    Min                 -0.05052      -0.04895        -0.039848
%    Max                 0.078096      0.049768         0.052764
%    Sum                0.0088302     -0.052128       -0.0056601
%    Count                      8             8                8

%Interpretation
%   The small value portfolio was the only portfolio with positive mean returns.
%   On the other hand, the small growth portolio was the only portfolio with positive median returns.
%   All the portfolios had average returns close to 0, and it is hard to draw conclusions given the small medians/means and relatively high standard deviations.
%   The small value portfolio's standard deviation, sample variance, kurtosis, and range suggested that it had the highest variability and most extreme annual returns.
%   The small sample size, though, further impedes our ability to draw conclusions from the data.
%   Each portfolio is skewed to the right (meaning higher values), but the magnitudes of the skew are less than one, meaning that the returns are not highly asymmetrical.
%
%   Looking at mean and max returns, the small value portfolio had the highest returns over our estimation period.
%   This could be explained by the combination of smaller companies' tendency to grow much faster than larger ones
%   and value companies' history of outgrowing growth companies (which is surprisingly counter to expectations).
%   Small companies have higher risk that is compensated with higher returns, and value companies are commonly thought to be purchased at a discount as compared to other companies.

%% C: Create a separate time series chart that overlays the Big Value, Big Neutral, and Big Growth portfolios.
    figure('Name','Big Portfolios')
    
    BV_ts = timeseries(BigValueReturns(12:end),1:8);
    BV_ts.TimeInfo.Units = 'years';
    BV_ts.TimeInfo.StartDate = '2010';     % Set start date.
    BV_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    BV_ts.Time = BV_ts.Time - BV_ts.Time(1);        % Express time relative to the start date.

    plot(BV_ts)
    grid on
    hold on
    
    BN_ts = timeseries(BigNeutralReturns(12:end),1:8);
    BN_ts.TimeInfo.Units = 'years';
    BN_ts.TimeInfo.StartDate = '2010';     % Set start date.
    BN_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    BN_ts.Time = BN_ts.Time - BN_ts.Time(1);        % Express time relative to the start date.
    
    plot(BN_ts)
    
    BG_ts = timeseries(BigGrowthReturns(12:end),1:8);
    BG_ts.TimeInfo.Units = 'years';
    BG_ts.TimeInfo.StartDate = '2010';     % Set start date.
    BG_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    BG_ts.Time = BG_ts.Time - BG_ts.Time(1);        % Express time relative to the start date.
    
    plot(BG_ts)
    hold off 
    
    title('Time Series: Big Portfolios')
    xlabel('Year')
    ylabel('Returns')
    legend('Big Value','Big Neutral','Big Growth','Location','northwest')    

%% D: Provide descriptive statistics of the Big Value, Big Neutral, and Big Growth portfolios. Present the descriptive statistics in a table. Interpret your results 1-2 paragraphs.
    BigVars = {'BigValue', 'BigNeutral','BigGrowth'};
    
    BV_stats = descriptiveStats(BigValueReturns(12:end));
    BN_stats = descriptiveStats(BigNeutralReturns(12:end));
    BG_stats = descriptiveStats(BigGrowthReturns(12:end));
    
    BigStats = array2table([BV_stats BN_stats BG_stats],'RowNames',StatRows,'VariableNames',BigVars)
    
%                      BigValue     BigNeutral    BigGrowth
%                      _________    __________    _________
%
%    Mean              0.0017501     0.013707      0.12294 
%    StdError           0.032493     0.012086      0.11431 
%    Median            -0.015433     0.006968      0.16111 
%    Mode                -0.1184    -0.025206     -0.57009 
%    StdDeviation       0.091903     0.034184      0.32331 
%    SampleVariance    0.0084461    0.0011685      0.10453 
%    Kurtosis             4.4981       2.3269       3.7558 
%    Skewness             1.3269      0.64681      -1.2119 
%    Range               0.32479      0.10135        1.048 
%    Min                 -0.1184    -0.025206     -0.57009 
%    Max                 0.20639     0.076141      0.47788 
%    Sum                0.014001      0.10966      0.98356 
%    Count                     8            8            8

%Interpretation
%   The big growth portfolio had both the highest mean and median returns.
%   Big value had the lowest mean and median returns.
%   The variance, standard error, standard deviation, and range for the big growth portfolio is highest.
%   This suggests high variability in the big growth returns.
%   Big value had the highest kurtosis, suggesting relatively many extreme values.
%   Big value is right-skewed (higher values) and big growth is left-skewed (lower values).
%   Big neutral was relatively symmetric.
%   Big growth had wild swings in returns, with both the most negative and most positive annual returns.
%
%   The findings for the big growth portfolio are surprising, given the assuption that large companies are supposedly more stable than smaller ones.
%   It is not surprising that the growth portfolio gives higher returns than the neutral and value portfolios, though.
%   The small sample size may have affected the accuracy and precision of our statistics.

%% E: Create a chart that overlays the SMB and HML factors.
    figure('Name','SMB and HML Factors')
    
    SMB_ts = timeseries(SMB(12:end),1:8);
    SMB_ts.TimeInfo.Units = 'years';
    SMB_ts.TimeInfo.StartDate = '2010';     % Set start date.
    SMB_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    SMB_ts.Time = SMB_ts.Time - SMB_ts.Time(1);        % Express time relative to the start date.

    plot(SMB_ts)
    grid on
    hold on
    
    HML_ts = timeseries(HML(12:end),1:8);
    HML_ts.TimeInfo.Units = 'years';
    HML_ts.TimeInfo.StartDate = '2010';     % Set start date.
    HML_ts.TimeInfo.Format = 'yyyy';       % Set format for display on x-axis.
    
    HML_ts.Time = HML_ts.Time - HML_ts.Time(1);        % Express time relative to the start date.
    
    plot(HML_ts)
    hold off 
    
    title('Time Series: SMB and HML Factors')
    xlabel('Year')
    ylabel('Returns')
    legend('SMB','HML','Location','northwest')


%% F: Provide descriptive statistics of both SMB and HML factors. Present the descriptive statistics in a table. Interpret your results 1-2 paragraphs.
    FactorVars = {'SMB', 'HML'};

    SMB_stats = descriptiveStats(SMB(12:end));
    HML_stats = descriptiveStats(HML(12:end));
    
    FactorStats = array2table([SMB_stats HML_stats],'RowNames',StatRows,'VariableNames',FactorVars)

%                          SMB          HML   
%                      _________    _________
%
%    Mean              -0.048174    -0.059692
%    StdError           0.042228     0.049438
%    Median            -0.050607    -0.075552
%    Mode               -0.17486     -0.26479
%    StdDeviation        0.11944      0.13983
%    SampleVariance     0.014265     0.019553
%    Kurtosis             3.4055       3.6002
%    Skewness             1.0178      0.76223
%    Range               0.37401      0.48981
%    Min                -0.17486     -0.26479
%    Max                 0.19915      0.22502
%    Sum                -0.38539     -0.47753
%    Count                     8            8

%Interpretation
%   Both the mean and median for the SMB and HML factors are negative.
%   This suggests that during our estimation period big portfolios performed better than small portfolios
%   and low-value (i.e., growth) performed better than high-value (i.e., value) portfolios.
%   These results are counter to what one woud expect.
%   The small sample size, high standard deviation, and relatively large standard error means that it is hard to draw accurate and precise conclusions from our data.
%   
%   From the mean values, one could interpret a stock with 100% exposure to the SMB factor to lose 4.8% of annual returns.
%   The same process could lead one to the conclusion that a stock with 100% exposure to the HML factor would lose 6% of annual returns.
%   Again, the small sample size hinders the accuracy and precision of our data.