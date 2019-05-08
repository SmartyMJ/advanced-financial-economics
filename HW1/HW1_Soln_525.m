%%EC525 HW1 Solutions
%% Add in SP500 stocks

%%%%%for those of you who had complex numbers for returns!1%%%%%%%%
%ret=log(abs(returns1(1:end-1,:))./abs(returns1(2:end,:)));
%%%%%you did not check your data to see if you had negatives which caused this
%%%%%do the above to ensure no negatives when calculating returns from WRDS
%%%%%as WRDS gives negative prices sometimes!



%Housekeeping
    clear all; close all; clc; 
%Get the list of S&P500 as of Jan 2018
    %[~,List] = xlsread();
    %This is where you would put in your ticker list from wikipedia so you can run a loop over them

% Quandl had #67 BF.B as BF_B and #75 BRK.B as BRK_B
% This happens- sometimes the databases we are working with name stocks slightly differently
    List{67} = 'BF_B';
    List{75} = 'BRK_B';
    
%Provide the API Key     
    Quandl.api_key('xaFxr9SP6Wd5sKFHdEax'); %This is the QFE key. DO NOT DISTRIBUTE

%Settings for data pull 
    database = 'WIKI'; %Could also use EOD
    startdate = '2011-12-30'; 
    enddate = '2017-12-31'; 
    type='data';

    %ticker loop
    for ii = 1:505
        ticker = List{ii}; 
        code = [database,'/',ticker];
        WikiData2=Quandl.get(code,'start_date',startdate,'end_date',enddate,'collapse','monthly','type',type);
        x{ii}=WikiData2(:,12);%store in cells because some might be incomplete
    end
    
%% take out the stocks that pose issues with data i.e. incomplete data or empty data
    y = [];
    for ii = 505:-1:1
        if size(x{ii},1)~=73 %get those stocks
            List(ii) = []; % delete the names of those stocks
        else
            y = [x{ii} y];% y will be the matrix containing all the data
        end
    end
%% need to flip the data
    y = flip(y);
%% calculate the returns
ret_y = (log(y([2:end],:))-log(y([1:end-1],:)));% monthly log return for stocks
    
    %% load the data for 10yr treasury yield, 10yr-2yr spread, sp500 dividend yield, and %change industrial production
    % convert dates
    %this will be different for everyones code depending on what you call it....
    %I call mine--
    %data_10yr
    %data_10Y2Y 
    %data_DY
    %data_Ind
    
%% calculate the average returns and betas
    avg_ret = mean(ret_y);
    
    %%first pass over columns of returns
    betas=zeros(4,size(ret_y,2));
    for ii = 1:469
        [B,TSTAT,S2,VCV,VCV_WHITE,R2,RBAR,YHAT] = ols(ret_y(:,ii),[data_10yr,data_10Y2Y,data_DY,data_Ind],1); %conduct the multivariate regression
        betas(:,ii) = B(2:5,1); %store the betas
    end
%% create table for 1a and 2a (the answers for them are the same)
    Results = array2table([avg_ret' betas'],'VariableNames',{'Average_return','Beta_10yr','Beta_10Y2Y','Beta_DY','Beta_Ind'}); %create a table
    Results.Properties.RowNames = List % this is the table for both 1a and 2a
    % find statistics of those variables
    avg_data = mean([avg_ret' betas'])
    std_data = std([avg_ret' betas'])
    tstat_data = avg_data./std_data*sqrt(size(ret_y,1))
    min_data = min([avg_ret' betas'])
    max_data = max([avg_ret' betas'])
    
%             Average_return     Beta_10yr     Beta_10Y2Y       Beta_DY      Beta_Ind 
%             ______________    ___________    ___________    ___________    _________
%
%   mean          0.0125          -0.6806        0.8225        -2.9940       -0.2106
%   std           0.0075           3.0138        1.8423        10.9134        1.4609 
%   tstat        14.1192          -1.9161        3.7883        -2.3279       -1.2233
%   min          -0.0216         -15.8257       -6.2187       -53.7922       -5.4717
%   max           0.0412          10.1857       10.0086        28.6721        4.7738



%% 2nd step of two-pass model
% regress average return for all stocks calling avg ret from above

    [B,TSTAT,S2,VCV,VCV_WHITE,R2,RBAR,YHAT] = ols(avg_ret',betas',1);
    Lamda2step = [B TSTAT];
    Results_2step = array2table(Lamda2step,'VariableNames',{'Lamda','TSTAT'},'RowNames',{'Lamda0','10yr','10Y2Y','DY','Ind'}) %create a table
    
%                 Lamda        TSTAT 
%              ___________    _______
%
%    Lamda0       0.012206      31.898
%    10yr        0.0006041      2.2806
%    10Y2Y     -0.00017146    -0.47954
%    DY        -0.00026561     -5.0269
%    Ind       -0.00018997     -0.5963


%% Fama-Macbeth Method

    LamdaF = [];
    STAT = [];
    for ii = 1:size(ret_y,1)
        ret2 = ret_y(ii,:)';
        [B,TSTAT,S2,VCV,VCV_WHITE,R2,RBAR,YHAT] = ols(ret2,betas',1);
        LamdaF = [LamdaF B];
    end
    LamdaFama = mean(LamdaF')';
    STDTFama = std(LamdaF')';
    TSTATFama = LamdaFama./STDTFama*sqrt(size(ret_y,1));
    Results_Fama = array2table([LamdaFama TSTATFama],'VariableNames',{'Lamda','TSTAT'},'RowNames',{'Lamda0','10yr','10Y2Y','DY','Ind'}) %create a table

%                 Lamda        TSTAT  
%              ___________    ________
%
%    Lamda0       0.012206      3.7116
%    10yr        0.0006041      0.9081
%    10Y2Y     -0.00017146    -0.22808
%    DY        -0.00026561     -1.6378
%    Ind       -0.00018997    -0.28699


%% Plot monthly Lamda 2b

% plot Lamda0
    figure
    hold on
    plot(date,LamdaF(1,:));
    datetick('x','mmmyy')
    xlabel('Year');
    ylabel('Lamda');
    title('Lamda0');
    hold off;

% plot Lamda for 10yr treasury yield
    figure
    hold on
    plot(date,LamdaF(2,:));
    datetick('x','mmmyy')
    xlabel('Year');
    ylabel('Lamda');
    title('Lamda for 10yr treasury yield');
    hold off;
    
% plot Lamda for 10yr-2yr spread
    figure
    hold on
    plot(date,LamdaF(3,:));
    datetick('x','mmmyy')
    xlabel('Year');
    ylabel('Lamda');
    title('Lamda for 10yr-2yr spread');
    hold off;

% plot Lamda for sp500 dividend yield
    figure
    hold on
    plot(date,LamdaF(4,:));
    datetick('x','mmmyy')
    xlabel('Year');
    ylabel('Lamda');
    title('Lamda for sp500 dividend yield');
    hold off;

% plot Lamda for %change industrial production
    figure
    hold on
    plot(date,LamdaF(5,:));
    datetick('x','mmmyy')
    xlabel('Year');
    ylabel('Lamda');
    title('Lamda for %change industrial production');
    hold off;
