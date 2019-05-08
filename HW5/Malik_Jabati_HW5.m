%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon GDP.xlsx, 10_Industry_Portfolios.xlsx, MarketReturns.xlsx, and the MIDASv2.3 folder.
    %Find interpretations for 1e in Malik_Jabati_HW5_1e.docx.
%Author:
    %Malik Jabati — 28Feb2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Assumptions:
    %Used value-weighted returns for Fama-French 10 industry portfolios
    %Used value-weighted returns for S&P 500 market index returns

%% Load in data
%%

%HouseKeeping
    clear all; close all; clc
    
    addpath(genpath('MATLAB Drive/HW5/MIDASv2.3')); %Add MIDAS toolbox to path
    
    %I went through all the data in Excel, and it was clean. There were no NaNs, extreme values, or recurring zeros.
    [GDP_num,GDP_txt,GDP_raw] = xlsread('GDP.xlsx','Sheet1');
    
    [Portfolios_num,Portfolios_txt,Portfolios_raw] = xlsread('10_Industry_Portfolios.xlsx','Sheet1');
    
    [Index_num,Index_txt,Index_raw] = xlsread('MarketReturns.xlsx','Sheet1');
    Index_num(:,2) = Index_num(:,2)*100; %convert to percentage to standardize dataset
    
    %Convert dates to text
    GDP_dates = datetime(GDP_num(:,1),'ConvertFrom','datenum') + calyears(1900);
    GDP_dates = cellstr(datestr(GDP_dates,'mm/dd/yyyy'));
    
    Portfolios_dates = datetime(Portfolios_num(:,1),'ConvertFrom','datenum') + calyears(1900);
    Portfolios_dates = cellstr(datestr(Portfolios_dates, 'mm/dd/yyyy'));
    
    Index_dates = datetime(Index_num(:,1),'ConvertFrom','datenum') + calyears(1900);
    Index_dates = cellstr(datestr(Index_dates,'mm/dd/yyyy'));   
    
    %Standardize variable names
    GDP_Y = GDP_num(:,2);
    Portfolios_X = Portfolios_num(:,2:end);
    Index_X = Index_num(:,2);
    
    %Quick test to see if data is clean
    dirtyGDP = isnan(GDP_Y);
    dirtyPortfolios = isnan(Portfolios_X);
    dirtyIndex = isnan(Index_X);
    
    if(any(dirtyGDP,'all') == 1 || any(dirtyPortfolios,'all') == 1 || any(dirtyIndex,'all') == 1)   %check for any NaN values
        disp('WARNING: Missing values in data!')
    elseif(min(GDP_Y,[],'all')<=-50 || min(Portfolios_X,[],'all')<=-50 || min(Index_X,[],'all')<=-50) %check for any losses of 50+%
        disp('WARNING: Check data for unexpectedly small values!')
    elseif(max(GDP_Y,[],'all')>=50 || max(Portfolios_X,[],'all')>=50 || max(Index_X,[],'all')>=50) %check for any gains of 50+%
        disp('WARNING: Check data for unexpectedly large values!')
    else
        disp('Data appears to be clean.')
    end
    
    % Specify lag structure and sample size
    Xlag = 3;
    Ylag = 1;
    EstStart = '04/01/1947';
    EstEnd = '10/01/2016';
    Method = 'rollingwindow';
    
%% 1.a Run the model for each portfolio and each horizon. Create a table that has industry names as rows, and horizons {3, 2, 1} as column headings, where table entries are the RMSE. See table below.
%%

    %For a horizon of 3   
    Horizon = 3;
    
    Portfolios_forecasts3 = cell(10,1);
    
    for i=1:10 %10 industry portfolios
        Portfolios_forecasts3{i} = MIDAS_ADL(GDP_Y,GDP_dates,Portfolios_X(:,i),Portfolios_dates,'Xlag',Xlag,'Ylag',Ylag,...
            'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
    end
    
    Portfolios_RMSE3 = zeros([10 1]);
    
    for i=1:10 %10 industry portfolios
        Portfolios_RMSE3(i) = Portfolios_forecasts3{i,1}.RMSE;
    end
    
    %For a horizon of 2   
    Horizon = 2;
    
    Portfolios_forecasts2 = cell(10,1);
    
    for i=1:10 %10 industry portfolios
        Portfolios_forecasts2{i} = MIDAS_ADL(GDP_Y,GDP_dates,Portfolios_X(:,i),Portfolios_dates,'Xlag',Xlag,'Ylag',Ylag,...
            'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
    end
    
    Portfolios_RMSE2 = zeros([10 1]);
    
    for i=1:10 %10 industry portfolios
        Portfolios_RMSE2(i) = Portfolios_forecasts2{i,1}.RMSE;
    end
 
    %For a horizon of 1   
    Horizon = 1;
    
    Portfolios_forecasts1 = cell(10,1);
    
    for i=1:10 %10 industry portfolios
        Portfolios_forecasts1{i} = MIDAS_ADL(GDP_Y,GDP_dates,Portfolios_X(:,i),Portfolios_dates,'Xlag',Xlag,'Ylag',Ylag,...
            'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
    end
    
    Portfolios_RMSE1 = zeros([10 1]);
    
    for i=1:10 %10 industry portfolios
        Portfolios_RMSE1(i) = Portfolios_forecasts1{i,1}.RMSE;
    end
    
    RmseForecastResults = array2table([Portfolios_RMSE3 Portfolios_RMSE2 Portfolios_RMSE1],...
        'RowNames',Portfolios_txt(1,2:end),'VariableNames',{'Hor_3','Hor_2', 'Hor_1'})
    
        
%% 1.b  Construct a forecast combination of the results above (using the default combination scheme).
%%

    %Horizon of 3
    Combo3 = ForecastCombine(Portfolios_forecasts3{1},Portfolios_forecasts3{2},Portfolios_forecasts3{3},Portfolios_forecasts3{4},...
        Portfolios_forecasts3{5},Portfolios_forecasts3{6},Portfolios_forecasts3{7},Portfolios_forecasts3{8},...
        Portfolios_forecasts3{9},Portfolios_forecasts3{10});
    
    Combo_RMSE3 = sqrt(mean((Combo3-GDP_Y(280:end)).^2));
    
    %Horizon of 2
    Combo2 = ForecastCombine(Portfolios_forecasts2{1},Portfolios_forecasts2{2},Portfolios_forecasts2{3},Portfolios_forecasts2{4},...
        Portfolios_forecasts2{5},Portfolios_forecasts2{6},Portfolios_forecasts2{7},Portfolios_forecasts2{8},...
        Portfolios_forecasts2{9},Portfolios_forecasts2{10});
    
    Combo_RMSE2 = sqrt(mean((Combo2-GDP_Y(280:end)).^2));
    
    %Horizon of 1
    Combo1 = ForecastCombine(Portfolios_forecasts1{1},Portfolios_forecasts1{2},Portfolios_forecasts1{3},Portfolios_forecasts1{4},...
        Portfolios_forecasts1{5},Portfolios_forecasts1{6},Portfolios_forecasts1{7},Portfolios_forecasts1{8},...
        Portfolios_forecasts1{9},Portfolios_forecasts1{10});
    
    Combo_RMSE1 = sqrt(mean((Combo1-GDP_Y(280:end)).^2));
    
%% 1.c Now run the same model from part a with market index returns instead of industry returns.
%%
    %For a horizon of 3   
    Horizon = 3;
    
    Index_forecast3 = MIDAS_ADL(GDP_Y,GDP_dates,Index_X,Index_dates,'Xlag',Xlag,'Ylag',Ylag,...
        'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
    
    Index_RMSE3 = Index_forecast3.RMSE;
    
    %For a horizon of 2   
    Horizon = 2;
    
    Index_forecast2 = MIDAS_ADL(GDP_Y,GDP_dates,Index_X,Index_dates,'Xlag',Xlag,'Ylag',Ylag,...
        'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
    
    Index_RMSE2 = Index_forecast2.RMSE;
    
    %For a horizon of 1   
    Horizon = 1;
    
    Index_forecast1 = MIDAS_ADL(GDP_Y,GDP_dates,Index_X,Index_dates,'Xlag',Xlag,'Ylag',Ylag,...
        'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','beta','Method',Method);
    
    Index_RMSE1 = Index_forecast1.RMSE;
    
%% 1.d Append two rows to the bottom of the table created in part a. Input the MSE of the default combination scheme and the MSE from the results in part c.
%%
    
    ComboRMSEs = array2table([Combo_RMSE3 Combo_RMSE2 Combo_RMSE1],...
        'RowNames',{'combo default'},'VariableNames',{'Hor_3','Hor_2', 'Hor_1'});
    
    IndexRMSEs = array2table([Index_RMSE3 Index_RMSE2 Index_RMSE1],...
        'RowNames',{'market index'},'VariableNames',{'Hor_3','Hor_2', 'Hor_1'});
    
    RmseForecastResultsFinal = [RmseForecastResults;ComboRMSEs;IndexRMSEs]
  
%% 1.e Interpret your findings 1-2 paragraphs.
%%

%See Malik_Jabati_HW5_1e.docx