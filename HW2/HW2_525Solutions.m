%  Solutions for EC525 HW2 Constructing your own Fama French Factors
%  This code is just one way to do this, I will note other acceptable functions or methods
%  it is not exhaustive, there may be other acceptable methods.
%  There are two different discussions provided for each question 
%  Each one is an example of what a good answer could have been
%  The discussions are not exhaustive, there could obviously be different discussions
%  that are still very much so correct. However, the two different answers below
%  illustrate a discussion that successfully touches on everything.


% HouseKeeping & Connect to SQLite
    clear all; close all; clc
 
%% Create a read-only SQLite connection to database
    conn = sqlite('tutorial.db','readonly');% create the connection with SQL
    Year1 = fetch(conn,'SELECT * FROM Year'); % select the dataset needed from SQL
    Ticker = fetch(conn,'SELECT * FROM Ticker');
    Prices = fetch(conn,'SELECT * FROM prices');
    NumShare = fetch(conn,'SELECT * FROM NumShare');
    BookValue = fetch(conn,'SELECT * FROM BookV');
    
%% Turn Cell to number for calculation, select the year needed
% Pulling specifically for our data period history, could have used all years not a big deal if you did

    Year=cell2mat(Year1(13:20));     
    Temp_Price=cell2mat(Prices(:,12:20)); % This selected 2009-2017 for the purpose of return calculation
    %remember, when you are calculating returns from prices you need to pull price from t-1 otherwise
    %you will have one observation less or your time will not match
    
    Prices=cell2mat(Prices(:,13:20));
    NumShare=cell2mat(NumShare(:,13:20));
    BookValue=cell2mat(BookValue(:,13:20));
    
%% Check Negative Prices
% I found that there are some negative values in prices and BV, so to clean the data, I use a checker
% to perform a logical statement and select the data needed

% many of you just made the prices positive-- this is okay because I simulated the data here
% I said you can pretend the prices are from WRDS and the negative prices are the average of bid ask because
% closing prices are missing, however you cannot do this all the time
% you also might get -999 from WRDS of NaN, so you need to take a look at your data before you move on

% here I run a quick check on my data, looking for anything weird
% I could also plot a time series of prices and do a visual check
% the min and max will tell me if I need to flag anything or if
% I have negative data points
% I am using the mean function intentionally because it will return
% NaN in any column that I might have NaN in
% there is a function nanmean that will take a mean and ignore NaNs

checkpr=[min(Prices) max(Prices) mean(Prices)];
checkbook=[min(BookValue) max(BookValue) mean(BookValue)];
checkshare=[min(NumShare) max(NumShare) mean(NumShare)];
checktemcal=[min(Temp_Price) max(Temp_Price) mean(Temp_Price)];

    check_price=ones(size(Ticker));
    for i=1:size(Ticker)
        for j=1:8
            if Prices(i,j)<=0 | BookValue(i,j)<=0 %this is important because French throws out any firm with negative
            check_price(i)=0;                    %book value so I am looking for either or  
            end
        end
    end
    
    check_price=logical(check_price); % This is used to turn the (zeros, ones) array to a logical statement
    
%% Select the needed Data
% This process is used to select the clean data according to the checker
% If there is a negative price or book value, I will throw everything out

   Prices=Prices(check_price,:);
   NumShare=NumShare(check_price,:);
   BookValue=BookValue(check_price,:);
   Ticker=Ticker(check_price,:);
   Temp_Price=Temp_Price(check_price,:);
   
%% Calculate MarketCap and BookToMarket

   TotalBookVal=BookValue.*NumShare; %Total Book Value is the pershare BV* shares
   MarketCap = Prices.*NumShare;
   BookToMar=TotalBookVal./MarketCap; % the BTM is calculated by total BV/Market Capitalization
   %or could've done BTM=bookvalue/price as book value was per share here
   
%% Get the returns

   Return = 100*tick2ret(Temp_Price);%I multiply by 100 to avoid scaling issues and to have as a percent for
   %ease of readability and interpretation
   %don't work with decimals convert everything to percent, unless you are working with simple returns of (1+r)

%% I write a loop to find the portfolios each year
% The following is the process of the loop so that you can follow the logic
% many of you used intersect functions, this is find too..
% 1. for each time t, I calculate the percentiles for market cap and btm
% 2. for each asset I check to see where they belong with logical statements
% 3. if the statement is true then I throw the return of that asset and the market cap into arrays
% 4. then calculate the weights of each asset based on the market cap to get value weighted 
%    as is the most common way. 

   S_Value=zeros(size(Year));
   S_Neut=zeros(size(Year));
   S_Growth=zeros(size(Year));
   B_Value=zeros(size(Year));
   B_Neut=zeros(size(Year));
   B_Growth=zeros(size(Year));
   
for i=1:8 % each year
        %could use quantile or percentile in matlab
        
       break_mc=quantile(MarketCap(:,i),0.5);
       break_bmtop=quantile(BookToMar(:,i),0.7);
       break_bmbot=quantile(BookToMar(:,i),0.3);
       
       %I clear the temporary stuff I am storing each loop to avoid issues
       SVal=[];
       SNeu=[];
       SGro=[];
       BVal=[];
       BNeu=[];
       BGro=[];
       
       WeSVal=[];
       WeSNeu=[];
       WeSGro=[];
       WeBVal=[];
       WeBNeu=[];
       WeBGro=[];
       
                   for j=1:size(Prices)
                       if MarketCap(j,i)>=break_mc && BookToMar(j,i)>=break_bmtop
                           BVal=[BVal Return(j,i)];
                           WeBVal=[WeBVal MarketCap(j,i)];
                       end
                       if MarketCap(j,i)>=break_mc && break_bmtop>BookToMar(j,i)>=break_bmbot
                           BNeu=[BNeu Return(j,i)];
                           WeBNeu=[WeBNeu MarketCap(j,i)];
                       end
                       if MarketCap(j,i)>=break_mc && break_bmbot>BookToMar(j,i)
                           BGro=[BGro Return(j,i)];
                           WeBGro=[WeBGro MarketCap(j,i)];
                       end
                       if MarketCap(j,i)<break_mc && BookToMar(j,i)>=break_bmtop
                           SVal=[SVal Return(j,i)];
                           WeSVal=[WeSVal MarketCap(j,i)];
                       end
                       if MarketCap(j,i)<break_mc && break_bmtop>BookToMar(j,i)>=break_bmbot
                           SNeu=[SNeu Return(j,i)];
                           WeSNeu=[WeSNeu MarketCap(j,i)];
                       end
                       if MarketCap(j,i)<break_mc && break_bmbot>BookToMar(j,i)
                           SGro=[SGro Return(j,i)];
                           WeSGro=[WeSGro MarketCap(j,i)];
                       end
                   end
       
           %I create my weights here
           WeBVal=WeBVal./sum(WeBVal);
           WeBNeu=WeBNeu./sum(WeBNeu);
           WeBGro=WeBGro./sum(WeBGro);
           WeSVal=WeSVal./sum(WeSVal);
           WeSNeu=WeSNeu./sum(WeSNeu);
           WeSGro=WeSGro./sum(WeSGro);
           
  
           SVal=sum(WeSVal.*SVal);
           SNeu=sum(WeSNeu.*SNeu);
           SGro=sum(WeSGro.*SGro);
           BVal=sum(WeBVal.*BVal);
           BNeu=sum(WeBNeu.*BNeu);
           BGro=sum(WeBGro.*BGro);
           
           %Finally I store to matrices I initialized outside of my loop
           S_Value(i)=SVal;
           S_Neut(i)=SNeu;
           S_Growth(i)=SGro;
           B_Value(i)=BVal;
           B_Neut(i)=BNeu;
           B_Growth(i)=BGro;
end
    
%% Time series overlay SVal, SNeu, SGro

figure
hold on
plot(Year,S_Value);
plot(Year,S_Neut);
plot(Year,S_Growth);
legend('Small Value','Small Neutrual','Small Growth');
title('Small Value/Neutral/Growth Returns')
xlabel('Year');
ylabel('Value-weighted Return (%)')
hold off

%% descriptive statistics

% I choose to use average, std, max, min, skewness and Kurtosis to explain their pattern.
% these are the most common ones chosen
% I stack them first
% then I put them all in a table together

SmallStat=[S_Value S_Neut S_Growth];
AvgS=mean(SmallStat);
StdS=std(SmallStat);
MaxS=max(SmallStat);
MinS=min(SmallStat);
SkewS=skewness(SmallStat);
KS=kurtosis(SmallStat);
SmallStat=array2table(vertcat(AvgS,StdS,MaxS,MinS,SkewS,KS));
SmallStat.Properties.VariableNames={'SmallValue' 'SmallNeutral' 'SmallGrowth'};
SmallStat.Properties.RowNames={'Average' 'Std' 'Max' 'Min' 'Skewness' 'KS'};
SmallStat

%% Interpretation for 1(b)

%%NOTE: BOTH DISCUSSIONS BELOW ARE EXAMPLES OF WHAT CONSTITUTES A CORRECT AND COMPLETE ANSWER TO THE QUESTION.

%%%%%%%%First discussion of 1(b) 

% For similarities, I found that most of the Small size companies have negative returns on average
% and their minimum is significantly low , which indicates that small firms
% are in general underperform than the big firms and they are more risky. 
% For the pattern across their BTM values, on average, small growth outperform than the small value
% because they generally have higher expectation of profits from investors and they are favorred by
% the market. However, they have higher standard deviation, which means they are more risky. This can
% be explain by the fact that investors might need to be compensated by the risk-return trade-off and 
% with high expectation, if their growth plans don't materialize, the price would plummet.

% In the maximum of growth happens during 2011 and 2017: 2011 is the recovery period from 2008 financial crisis
% where a lot of small growth firms enjoys the benefits, while in 2017, the huge tax cut lead to high pofits
% from small growth firms and the stock market is performed well during the time. 
% In time series, skewness measures the probability of amplitude, the negative skewness for Small Neutral means
% they have more higher value in general and relatively stable. The kurtosis for all of them are small due to 
% the limit of data size.


%%%%%%%%Second discussion of 1(b) 
%Interpret the results: 

% The table above provides basic descriptive statistics for the three small portfolios. It is important to
% note that the study only calculated return for eight years.
% Overall the statistics follow what we studied in class about the difference between Value, Neutral , 
% and Growth. Value companies are categorized as ones with large amounts of real capital and low amounts of
% new investment, which is why they have high book to market ratios. They normally encompass established less risky 
% blue chip companies. They generally exhibit consistent, yet low, earnings and stock returns.
% This definitions is partially supported by the data, as the small value portfolios did, in fact, 
% have the lowest average return yet they also had the highest standard deviation and largest spread. 
% Here spread is defined as the difference between the max and min returns. 

% On the other hand, growth companies have more investments (debt) and they are expected to grow at a 
% higher pace and therefore have greater return. They also do not hold as much real capital, so they
% have a low book to market ratio. However, these generally 'tech' and or 'start up' companies, are more
% risky so investors should receive a premium for their investment. Again, the data partially support
% this definition. Growth assets enjoyed the highest returns and lowest standard deviation and spread.
% Here, growth companies also had the highest skew, meaning that they had more highly successful outliers. 
% Finally, the the neutral portfolio encompasses stocks that are not definitively measured as value or growth.
% Therefore, we expected them to exhibit medium returns and medium risk. However, the data showed that the neutral
% portfolio had the lowest average return, but medium spread and standard deviation.


       
%% 1(c) BVal BNeu BGro Plot

figure % The plot process is same as Small portfolios
hold on
plot(Year,B_Value)
plot(Year,B_Neut);
plot(Year,B_Growth);
legend('Big Value','Big Neutrual','Big Growth');
title('Big Value/Neutral/Growth Returns')
xlabel('Year');
ylabel('Value-weighted Return (%)')
hold off

%% 1(d) Descriptive Statistics

BigStat=[B_Value B_Neut B_Growth]; %The stat process is same as Small portfolios
AvgS=mean(BigStat);
StdS=std(BigStat);
MaxS=max(BigStat);
MinS=min(BigStat);
SkewS=skewness(BigStat);
KB=kurtosis(BigStat);
%I create a table and store everything
BigStat=array2table(vertcat(AvgS,StdS,MaxS,MinS,SkewS,KB));
BigStat.Properties.VariableNames={'BigValue' 'BigNeutral' 'BigGrowth'};
BigStat.Properties.RowNames={'Average' 'Std' 'Max' 'Min' 'Skewness' 'Kurtosis'};
BigStat

%% Interpretation for 1(d) 

%%NOTE: BOTH DISCUSSIONS BELOW ARE EXAMPLES OF WHAT CONSTITUTES A CORRECT AND COMPLETE ANSWER TO THE QUESTION.

%%%%%%%%First discussion of 1(b) 

% For similarities, all BIG portfolios have positive average return since they have higher market expectations. 
% Here the Big Growth portfolios beat the Big Value portfolios again and even by higher amount comparing to small
% portfolios. Besides the reason stated in 1(b), one of the main reason is the superstar firm effect that change
% the market foundamentals in recent years. They growing fastly and have the ability to disrups other industries.
% This phenomenon leads to the larger spread between Big Growth and Big Value.
% The superstar firm effect can also be used to explain the abnormal greater std for value since the big growing
% firms has dominance in the market.
% The min of Big value and max of Big growth happens simultaneously in 2014 where the superstar effect is growing. And 
% growth's success is on the expense of value assets. Since the skewness and kurtosis are both small and similar
% from one another, they are not significant for the discussion.


%%%%%%%%Second discussion of 1(b) 
%Interpret the results: 

% Discussion of the differences of value, neutral, and growth portfolios of large companies is similar
% to the discussion for the three small portfolios in terms of the macro economic theory. However, here 
% the portfolios control for only large companies as opposed to small companies above. 

% Again, value companies should have lower return for lower risk (as measured here by the spread, standard
% deviation and skew). Large value companies do exhibit the lowest mean and median and average spread and 
% standard deviation. However, they also have the highest skew, meaning that this specific portfolio also
% has several assets with high returns. Next, we expect the growth portfolio to have high returns for 
% high risk. We do see this, in fact, the average returns on the growth portfolio are 61 times the average
% returns on the value portfolio. The growth portfolio also exhibits the highest standard deviation and 
% lowest skew of all 6 portfolios. This is largely driven by a handful of companies whose market caps 
% decreased drastically around 2016, only to increase back up the following year. Assuming this is not a data
% entry error, this phenomenon showcases weaknesses in a market value weighted portfolio. Lastly, as mentioned
% in the discussion above, we expect the neutral portfolio to exhibit average returns for average risk. Overall, 
% that is exactly what we see here with the neutral portfolio exhibiting average returns (average meaning less
% than growth but more than value portfolios) and average skew.


%% Constructing the factors SMB and HML

Small_Big=zeros(size(Year)); % Calculate the SMB and HML based on the formula given by French's website
%doing this for each year as per French's method
for i=1:8
    SMB=1/3*(S_Value(i)+S_Neut(i)+S_Growth(i))-1/3*(B_Value(i)+B_Neut(i)+B_Growth(i));
    Small_Big(i)=SMB;
end

High_Low=zeros(size(Year));
for i=1:8
    HML=1/2*(S_Value(i)+B_Value(i))-1/2*(S_Growth(i)+B_Growth(i));
    High_Low(i)=HML;
end

%% 1(e) Time series chart for SMB and HML

figure 
hold on
plot(Year,Small_Big)
plot(Year,High_Low);
legend('SMB','HML');
title('SMB and HML Returns')
xlabel('Year');
ylabel('Value-weighted Return (%)')
hold off

%% 1(f) Descriptive Data for SMB and HML
%Again I take a bunch of sample statistics and then store them in a table together

SMBHMLStat=[Small_Big High_Low]; 
AvgS=mean(SMBHMLStat);
StdS=std(SMBHMLStat);
MaxS=max(SMBHMLStat);
MinS=min(SMBHMLStat);
SkewS=skewness(SMBHMLStat);
KSum=kurtosis(SMBHMLStat);
SMBHMLStat=array2table(vertcat(AvgS,StdS,MaxS,MinS,SkewS,KSum));
SMBHMLStat.Properties.VariableNames={'SMB' 'HML'};
SMBHMLStat.Properties.RowNames={'Average' 'Std' 'Max' 'Min' 'Skewness' 'Kurtosis'};
SMBHMLStat

%% Interpretation of 1(f)

%%NOTE: BOTH DISCUSSIONS BELOW ARE EXAMPLES OF WHAT CONSTITUTES A CORRECT AND COMPLETE ANSWER TO THE QUESTION.

%%%%%%%%First discussion of 1(b) 

% The table above compares the High-minus-Low and Small-minus-Big factors.  
% It is interesting to note, that the correlation between SMB and HML is low. This is to be
% expected as the two factors should explain different aspects of the market and therefore have
% very little or no correlation. (In theory no correlation) 
% For SMB we expect the mean and median to be positive, meaning that smaller companies have greater 
% return than larger ones. However here we see that SMB has a negative mean and median meaning that 
% large companies outperform smaller ones. This could potentially be explained by the fact that the 
% represented years were a time of post financial crisis recovery, where larger companies
% that were hit hard during the recession regained their footing. The min, max, and spread simply measure
% the magnitude at which the two sizes differ and the second and third moments measure how their difference
% tends to move. In the data the skew for SMB is highly positive when compared to the HML factor meaning that this
% factor had a greater amount of highly positive values, which we would expect due to several high return 
% small assets. 

% For the HML factor we expect growth companies to outperform value since growth companies must pay 
% investors a premium for the additional risk, meaning that the mean and median of the
% factor should be negative. This phenomenon occurs in the data. However we also see a larger spread and 
% standard deviation in HML when compared to SMB, suggesting that the difference of value and growth 
% companies fluctuates and at times value outperforms growth. This can, again, be explained by the fact that
% the data reflects a period of post recession recovery. 



%%%%%%%%Second discussion of 1(b) 
%Interpret the results: 

% For SMB, it is on average negative, which means the Small portfolio usually has lower
% value-weighted average return than the large firms when they are nutrual to the value effects. 
% This can be explain by the superstar effect since the largest drop is during 2014. 
% For HML, it is also on average negative, which means the value portfolio outperform than the Gorwth
% portfolio when they are size neutral. It can be explained by the compensation of risky
% growth portfolio. However, in recent years, the HML is growing due to the superstar effect stated above.
% Furthermore, the large standard deviation for both of them means the volatility of SMB HML factor is 
% large and hard to predict through time series, so the main purpose is to explain the 
% returns rather than predict. Same as the Skewness and kurtosis. 



