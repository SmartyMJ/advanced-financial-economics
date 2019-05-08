%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon Estimize_Expectations.xlsx, tickhistory.mat, EventStudyNFP_Data.mat, and Malik_Jabati_HW3.docx.
%Author:
    %Malik Jabati — 11Feb2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Table Names:
    %
%Assumptions:
    %I used absolute in all of computations and visualizations.

%% Load in data
%HouseKeeping
    clear all; close all; clc
    
    %I went through all the data, and it was clean. There were no zeros and no negaitves.
    load EventStudyNFP_Data.mat %first column is dates second column is NFP release numbers
    load tickhistory.mat
    [num_est,txt_est,raw_est] = xlsread('Estimize_Expectations.xlsx');
    
    %Extract asset data from cells
    txt_asset=All_data{2,1}; %Text values for assets
    raw_asset=All_data{3,1}; %Raw (all) data for assets
    num_asset=All_data{1,1}; %Number values for assets
    
    %Determine total number of periods of interest
    total_dates = unique(num_asset(:,6)); %identify all unique dates
    daily_minutes = unique(txt_asset(2:end,6)); %identify all recorded minutes in day
    total_assets = unique(txt_asset(2:end,1)); %identify all unique dates
    
    temp_dates=datestr(NFP_dat(:,1),'mmmddyyyy');%convert dates to characters
    dates = cellstr(temp_dates);%convert to cell so that can be a variable name
    
    X = [NFP_dat(:,1) (NFP_dat(:,2)-num_est(:,2))]; %Create independent variable, i.e. NFP surprise
    
%% Create clean table with assets and prices    

    %Create and fill date column
    date_col = zeros([length(total_dates)*length(daily_minutes),1]);
    for i=1:length(total_dates)     
        for j=(i-1)*length(daily_minutes)+1:i*length(daily_minutes)
            date_col(j,1) = NFP_dat(i,1);
        end
    end

    %Create and fill time column  
    time_col = txt_asset(2:6336,6);
    
    %Create and fill price columns
    prices = reshape(num_asset(:,5),[length(total_dates)*length(daily_minutes),10]);
    
    %Create table with concatenated colums 
    AssetPrices = [array2table(date_col) array2table(time_col) array2table(prices)];
    AssetPrices.Properties.VariableNames = {'Date','Time','AUD','CAD','CHF', 'EUR', 'GBP', 'JPY', 'NZD', 'RBc1', 'XAG', 'XAU'};

%% Create table with average returns
    
    avg_returns = [NFP_dat(:,1) zeros(length(total_dates),10)]; %Get average minute returns for each day
    for i=1:length(total_dates)
        start = (i-1)*length(daily_minutes)+1;
        finish = start+length(daily_minutes)-1;
        avg_returns(i,2:end) = 100*mean(tick2ret(prices(start:finish,:))); %Get returns as percentage
        %avg_returns(i,2:end) = 100*mean(diff(log(prices(start:finish,:))));
    end

%% Estimate alpha and beta using first 32 periods

    coefficients = zeros([2,length(total_assets)]); %allocate space for coefficient array. Alphas are row 1, Betas row 2
    
    for i=1:length(total_assets)
        coefficients(:,i) = [ones(32,1) X(1:32,2)] \ avg_returns(1:32,1+i); % X \ y finds the coefficients for the linear regression that fits X to y
    end
    %Equivalent to fitlm(X(1:32,2),avg_returns(1:32,1+i))
    
%% 1. Conduct an event study for each asset
    
%Do this three times for each event: day 33, day 34, day 35
%Then combine table to find daily average
    
    event_time = (-5:30)';    

    %Create table of event 1 (day 33) abnormal returns for each asset
    %Time is -5 to +30 minutes in 1-minute increments. 8:25 = t(-5), 8:30 = t(0), 9:00 = t(+30)
    %8:25 is 85 minutes after the start of the day (7:00), so use 85 as the offset
    
    %Create an event window for days 33, 34, and 35
    
    %Day 33: event 1
    event_window1 = [event_time zeros([length(event_time),10])];
    day_start1 = (33-1)*length(daily_minutes)+1; % 7:00 price on day 33
    wstart1 = day_start1+85; % event window start is 8:25 on that day
    event_window1(:,2:end) = prices(wstart1:wstart1+35,:);
    
    %Simple test to ensure price values align for first asset
        %event_window1(1,2)
        %prices(5878,1) %equals (33-1)*181+1+85 which is the index of 8:25 on day 33
    
    %Day 34: event 2
    event_window2 = [event_time zeros([length(event_time),10])];
    day_start2 = (34-1)*length(daily_minutes)+1; % 7:00 price on day 34
    wstart2 = day_start2+85; % event window start is 8:25 on that day
    event_window2(:,2:end) = prices(wstart2:wstart2+35,:);   

    %Simple test to ensure price values align for first asset
        %event_window2(1,2)
        %prices(6059,1) %equals (34-1)*181+1+85 which is the index of 8:25 on the day 34
    
    %Day 35: event 3
    event_window3 = [event_time zeros([length(event_time),10])];
    day_start3 = (35-1)*length(daily_minutes)+1; % 7:00 price on day 35
    wstart3 = day_start3+85; % event window start is 8:25 on that day
    event_window3(:,2:end) = prices(wstart3:wstart3+35,:);        
        
    %Simple test to ensure price values align for first asset
        %event_window3(1,2)
        %prices(6240,1) %equals (35-1)*181+1+85 which is the index of 8:25 on the day 35
   
    %Create abnormal return window for each event
    
    %Event 1 (day 33) returns
    expected_return1 = coefficients(1,:) + coefficients(2,:)*X(33,2); % E(r) = alpha + beta*factor
    
    ret1 = zeros([length(event_time),10]); %array for regular returns
    abret1 = zeros([length(event_time),10]); %array for abnormal returns
    
    ret1(1,:)=(event_window1(1,2:end)-event_window1(1,2:end))./(event_window1(1,2:end)); %clearly start at 0
    for i=2:length(event_window1) %for each minute starting from t(-4)
        ret1(i,:)=(event_window1(i,2:end)-event_window1((i-1),2:end))./(event_window1((i-1),2:end)); %taking mt/mt returns percentage returns
        abret1(i,:)=ret1(i,:)-expected_return1;
    end  
        
    %Event 2 (day 34) returns
    expected_return2 = coefficients(1,:) + coefficients(2,:)*X(34,2); % E(r) = alpha + beta*factor
    
    ret2 = zeros([length(event_time),10]); %array for regular returns
    abret2 = zeros([length(event_time),10]); %array for abnormal returns
    
    ret2(1,:)=(event_window2(1,2:end)-event_window2(1,2:end))./(event_window2(1,2:end)); %clearly start at 0
    for i=2:length(event_window2) %for each minute starting from t(-4)
        ret2(i,:)=(event_window2(i,2:end)-event_window2((i-1),2:end))./(event_window2((i-1),2:end)); %taking mt/mt returns percentage returns
        abret2(i,:)=ret2(i,:)-expected_return2;
    end          
        
    %Event 3 (day 35) returns 
    expected_return3 = coefficients(1,:) + coefficients(2,:)*X(35,2); % E(r) = alpha + beta*factor
    
    ret3 = zeros([length(event_time),10]); %array for regular returns
    abret3 = zeros([length(event_time),10]); %array for abnormal returns

    ret3(1,:)=(event_window3(1,2:end)-event_window3(1,2:end))./(event_window3(1,2:end)); %clearly start at 0
    for i=2:length(event_window3) %for each minute starting from t(-4)
        ret3(i,:)=100*((event_window3(i,2:end)-event_window3((i-1),2:end))./(event_window3((i-1),2:end))); %taking mt/mt returns percentage returns
        abret3(i,:)=ret3(i,:)-expected_return3;
    end         
        
    %Compute average abnormal returns and cumulative avg. abnormal returns
    avg_abreturns = (abret1+abret2+abret3)/3;
    CAR = cumsum(avg_abreturns);
    
    %Combine avg. abnormal returns in CAR, interleaving according to asset
    event_study_array = zeros(length(CAR),20); %20 columns because there are 2 columns for each of the 10 assets
    for i=1:10 %for each asset
        event_study_array(:,2*(i-1)+1) = avg_abreturns(:,i); %first column is avg. abnormal returns
        event_study_array(:,2*(i-1)+2) = CAR(:,i); %second column is CAR
    end
    
    %Create cell array for variable names
    AssetNames = {'AUD','CAD','CHF', 'EUR', 'GBP', 'JPY', 'NZD', 'RBc1', 'XAG', 'XAU'};
    AssetNames2x = repelem(AssetNames,2);
    for i=1:10 %for each asset
        AssetNames2x(2*(i-1)+1) = strcat(AssetNames2x(2*(i-1)+1),{'_AR'});
        AssetNames2x(2*(i-1)+2) = strcat(AssetNames2x(2*(i-1)+2),{'_CAR'});
    end
        
    %Create cell array for row names
    event_time_rows = cellstr(strcat(num2str(event_time),'mts'));
    event_time_rows(7:end) = strcat({'+'},event_time_rows(7:end));
    
    %Create abnormal returns table
    EventStudy = array2table(event_study_array,'RowNames',event_time_rows,'VariableNames',AssetNames2x)
    
%% 2. Create a graph with the asset CARs on Y and event time on X

    plot(event_time, CAR(:,1:7)) %plot forex
    hold on
    plot(event_time, CAR(:,8:end),':') %plot commodities
    xlabel('Event time (minutes)')
    ylabel('Cumulative average abnormal return (%)')
    title('CAR over Time')
    legend(AssetNames,'Location','northwest')
    hold off
    
%% 3. Compute SCAR test statistic and associated p-value
    
    %Follows methodology from in-class  breakout
    scar_stat = zeros([1,10]);
    pvalue = zeros([1,10]);
    for i=1:10 %for each asset
        car_hat = CAR(end,i);
        var_ar = var([abret1(:,i) abret2(:,i) abret3(:,i)]);
        sigma_e=sqrt((1/(length(var_ar).^2))*sum(var_ar));
        scar_stat(i)=(car_hat)/sigma_e; %construct SCAR statistic
        pvalue(i) = 1-normcdf(abs(scar_stat(i)));
    end
    
    %Make a table with scar stat and p-value
    StatsTable = array2table([scar_stat; pvalue],'RowNames',{'SCAR','P-Value'},'VariableNames',AssetNames)
    
%% 4. Comment on results from the tables above
    % See Malik_Jabati_HW3.docx for answer to this question