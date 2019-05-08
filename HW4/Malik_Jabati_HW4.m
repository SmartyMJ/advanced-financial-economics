%Purpose:
    %Econ 525-Spring2019
%Note:
    %This m-file is dependent upon S&P_Top5.xlsx.
%Author:
    %Malik Jabati — 21Feb2019
    %UNC Honor Pledge: I certify that no unauthorized assistance has been received or given in the completion of this work.
%Assumptions:
    % Used daily HPY log returns
    % Used 50 epochs to train neural network, in the interest of speed and
    % flexibility for the grader
    
%% Load in data
%%
%HouseKeeping
    clear all; close all; clc
    
    %I went through all the data, and it was clean. There were no zeros and no negatives.
    [num,txt,raw] = xlsread('S&P_Top5.xlsx','Sheet1');
    
    %Prepare the data
    assets = txt(1,2:end);
    returns = diff(log(num))*100; %convert to percentage returns
    dates = txt(3:end,1);
    T = length(returns); %Count observations
    
    %Partition the data into training and evaluation
    numTimeStepsTrain = floor(0.9*T);
    dataTrain = returns(1:numTimeStepsTrain,:);
    dataTest = returns(numTimeStepsTrain+1:end,:);
    
    %Standardize the data for LTSM model
    mu = mean(dataTrain);
    sig = std(dataTrain);
    dataTrainStandardized = (dataTrain - mu) ./ sig;
    
%% 1.a Conduct a one day ahead forecast with an ARMA(1,1) model for each individual asset.
%%

    model1 = arima(1,0,1);

    MSFT_mdl = estimate(model1,dataTrain(:,1));
    AAPL_mdl = estimate(model1,dataTrain(:,2));
    AMZN_mdl = estimate(model1,dataTrain(:,3));
    BRKB_mdl = estimate(model1,dataTrain(:,4));
    FB_mdl = estimate(model1,dataTrain(:,5));  
   
    %Conduct one-day ahead forecast
    MSFT_fcst1 = forecast(MSFT_mdl,1,'Y0',dataTrain(:,1));
    AAPL_fcst1 = forecast(AAPL_mdl,1,'Y0',dataTrain(:,2));
    AMZN_fcst1 = forecast(AMZN_mdl,1,'Y0',dataTrain(:,3));
    BRKB_fcst1 = forecast(BRKB_mdl,1,'Y0',dataTrain(:,4));
    FB_fcst1 = forecast(FB_mdl,1,'Y0',dataTrain(:,5));
    
    
    
%% 1.b Conduct a one day ahead forecast with an LSTM model for each individual asset.
%%
    
    %Use the time t to predict t+1
    XTrain = dataTrainStandardized(1:end-1,:); %Lagged values of the data for predictors
    YTrain = dataTrainStandardized(2:end,:);
    
    %Generate the regression network
    numFeatures = 1;
    numResponses = 1;
    numHiddenUnits = 200;
    %Define the options for the LSTM network
    layers = [ ...
        sequenceInputLayer(numFeatures)
        lstmLayer(numHiddenUnits)
        fullyConnectedLayer(numResponses)
        regressionLayer];
    
    options = trainingOptions('adam', ...
        'MaxEpochs',50, ...
        'GradientThreshold',1, ...
        'InitialLearnRate',0.005, ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropPeriod',125, ...
        'LearnRateDropFactor',0.2, ...
        'Verbose',0, ...
        'Plots','training-progress');
    
    %Train the network
    
    MSFT_net = trainNetwork(XTrain(:,1)',YTrain(:,1)',layers,options);
    AAPL_net = trainNetwork(XTrain(:,2)',YTrain(:,2)',layers,options);
    AMZN_net = trainNetwork(XTrain(:,3)',YTrain(:,3)',layers,options);
    BRKB_net = trainNetwork(XTrain(:,4)',YTrain(:,4)',layers,options);
    FB_net = trainNetwork(XTrain(:,5)',YTrain(:,5)',layers,options);
    
    
%%
    %Forecast one day ahead using LSTM model

    %Standardize
    dataTestStandardized = (dataTest - mu) ./ sig;
    XTest = dataTestStandardized(1:end-1,:);
    
    %To initialize the network state, make the first prediction using the last time step of the training response YTrain(end). 
    %Loop over the remaining predictions and input the previous prediction to predictAndUpdateState.
    
    MSFT_YPred = predict(MSFT_net,YTrain(end,1));
    AAPL_YPred = predict(AAPL_net,YTrain(end,2));
    AMZN_YPred = predict(AMZN_net,YTrain(end,3));
    BRKB_YPred = predict(BRKB_net,YTrain(end,4));
    FB_YPred = predict(FB_net,YTrain(end,5));
    
    %[MSFT_net,MSFT_YPred] = predictAndUpdateState(MSFT_net,YTrain(end,1));
    %[AAPL_net,AAPL_YPred] = predictAndUpdateState(AAPL_net,YTrain(end,2));
    %[AMZN_net,AMZN_YPred] = predictAndUpdateState(AMZN_net,YTrain(end,3));
    %[BRKB_net,BRKB_YPred] = predictAndUpdateState(BRKB_net,YTrain(end,4));
    %[FB_net,FB_YPred] = predictAndUpdateState(FB_net,YTrain(end,5));
    
    
%% 1.c Repeat part a for 2 days ahead, 3 days ahead, 4 days ahead and 5 days ahead.
%%
    
    %Conduct five-day ahead forecast using ARMA model
    
    [MSFT_fcst5, MSFT_MSE] = forecast(MSFT_mdl,5,'Y0',dataTrain(:,1));
    [AAPL_fcst5, AAPL_MSE] = forecast(AAPL_mdl,5,'Y0',dataTrain(:,2));
    [AMZN_fcst5, AMZN_MSE] = forecast(AMZN_mdl,5,'Y0',dataTrain(:,3));
    [BRKB_fcst5, BRKB_MSE] = forecast(BRKB_mdl,5,'Y0',dataTrain(:,4));
    [FB_fcst5, FB_MSE] = forecast(FB_mdl,5,'Y0',dataTrain(:,5));
    
    %MSFT_fcst5 = zeros([1,5]);
    %AAPL_fcst5 = zeros([1,5]);
    %AMZN_fcst5 = zeros([1,5]);
    %BRKB_fcst5 = zeros([1,5]);
    %FB_fcst5 = zeros([1,5]);
    
    %for i=1:5
    %    [MSFT_fcst5(i), MSFT_MSE] = forecast(MSFT_mdl,1,'Y0',dataTrain(end-(i-1),1));
    %    [AAPL_fcst5(i), AAPL_MSE] = forecast(AAPL_mdl,1,'Y0',dataTrain(end-(i-1),2));
    %    [AMZN_fcst5(i), AMZN_MSE] = forecast(AMZN_mdl,1,'Y0',dataTrain(end-(i-1),3));
    %    [BRKB_fcst5(i), BRKB_MSE] = forecast(BRKB_mdl,1,'Y0',dataTrain(end-(i-1),4));
    %    [FB_fcst5(i), FB_MSE] = forecast(FB_mdl,1,'Y0',dataTrain(end-(i-1),5));
    %end
    
%% 1.d Repeat part b for 2 days ahead, 3 days ahead, 4 days ahead and 5 days ahead.
%%
    
    %Predict five days ahead using static LSTM model
    
    %Create vectors to store predicted values
    MSFT_YPred5 = zeros([1,5]);
    AAPL_YPred5 = zeros([1,5]);
    AMZN_YPred5 = zeros([1,5]);
    BRKB_YPred5 = zeros([1,5]);
    FB_YPred5 = zeros([1,5]);
    
    for i=1:5 %five steps ahead
        MSFT_YPred5(i) = predict(MSFT_net,YTrain(end-(i-1),1));
        AAPL_YPred5(i) = predict(AAPL_net,YTrain(end-(i-1),2));
        AMZN_YPred5(i) = predict(AMZN_net,YTrain(end-(i-1),3));
        BRKB_YPred5(i) = predict(BRKB_net,YTrain(end-(i-1),4));
        FB_YPred5(i) = predict(FB_net,YTrain(end-(i-1),5));
    end
    
    %for i=1:5 %five steps ahead
    %    [MSFT_net,MSFT_YPred5(i)] = predictAndUpdateState(MSFT_net,YTrain(end-(i-1),1));
    %    [AAPL_net,AAPL_YPred5(i)] = predictAndUpdateState(AAPL_net,YTrain(end-(i-1),2));
    %    [AMZN_net,AMZN_YPred5(i)] = predictAndUpdateState(AMZN_net,YTrain(end-(i-1),3));
    %    [BRKB_net,BRKB_YPred5(i)] = predictAndUpdateState(BRKB_net,YTrain(end-(i-1),4));
    %    [FB_net,FB_YPred5(i)] = predictAndUpdateState(FB_net,YTrain(end-(i-1),5));
    % end
    
    
    
%% 1.e Repeat part c, but now with a dynamic forecasting structure.
%%
    
    %Conduct five-day ahead dynamic forecast using ARMA model
    
    %Create vectors for dynamic forecasts, putting final training value
    %into first position
    
    MSFT_fcst5d = zeros([1,5]);
    AAPL_fcst5d = zeros([1,5]);
    AMZN_fcst5d = zeros([1,5]);
    BRKB_fcst5d = zeros([1,5]);
    FB_fcst5d = zeros([1,5]);
    
    %Initialize first values
    
    MSFT_fcst5d(1) = MSFT_fcst1;
    AAPL_fcst5d(1) = AAPL_fcst1;
    AMZN_fcst5d(1) = AMZN_fcst1;
    BRKB_fcst5d(1) = BRKB_fcst1;
    FB_fcst5d(1) = FB_fcst1;
    
    dataTrainMSFT = [dataTrain(:,1); MSFT_fcst1];
    dataTrainAAPL = [dataTrain(:,2); AAPL_fcst1];
    dataTrainAMZN = [dataTrain(:,3); AMZN_fcst1];
    dataTrainBRKB = [dataTrain(:,4); BRKB_fcst1];
    dataTrainFB = [dataTrain(:,5); FB_fcst1];
    
    MSFT_mdl = estimate(model1,dataTrainMSFT);
    AAPL_mdl = estimate(model1,dataTrainAAPL);
    AMZN_mdl = estimate(model1,dataTrainAMZN);    
    BRKB_mdl = estimate(model1,dataTrainBRKB);
    FB_mdl = estimate(model1,dataTrainFB);      

    %Forecast i-step ahead, append this result to training data, and
    %re-estimate model. Repeat until you reach five steps ahead
    for i=2:5 %five steps ahead
              
        [MSFT_fcst5d(i), MSFT_MSEd] = forecast(MSFT_mdl,1,'Y0',dataTrainMSFT);
        [AAPL_fcst5d(i), AAPL_MSEd] = forecast(AAPL_mdl,1,'Y0',dataTrainAAPL);
        [AMZN_fcst5d(i), AMZN_MSEd] = forecast(AMZN_mdl,1,'Y0',dataTrainAMZN);
        [BRKB_fcst5d(i), BRKB_MSEd] = forecast(BRKB_mdl,1,'Y0',dataTrainBRKB);
        [FB_fcst5d(i), FB_MSEd] = forecast(FB_mdl,1,'Y0',dataTrainFB);
        
        dataTrainMSFT = [dataTrainMSFT; MSFT_fcst5d(i)];
        dataTrainAAPL = [dataTrainAAPL; AAPL_fcst5d(i)];
        dataTrainAMZN = [dataTrainAMZN; AMZN_fcst5d(i)];
        dataTrainBRKB = [dataTrainBRKB; BRKB_fcst5d(i)];
        dataTrainFB = [dataTrainFB; FB_fcst5d(i)];
        
        MSFT_mdl = estimate(model1,dataTrainMSFT);
        AAPL_mdl = estimate(model1,dataTrainAAPL);
        AMZN_mdl = estimate(model1,dataTrainAMZN);    
        BRKB_mdl = estimate(model1,dataTrainBRKB);
        FB_mdl = estimate(model1,dataTrainFB);  
    end
%% 1.f Repeat part d, but now with a dynamic forecasting structure.
%%
%Conduct five-day ahead dynamic forecast using LSTM model

    %Create vectors to store predicted values
    MSFT_YPred5d = zeros([1,5]);
    AAPL_YPred5d = zeros([1,5]);
    AMZN_YPred5d = zeros([1,5]);
    BRKB_YPred5d = zeros([1,5]);
    FB_YPred5d = zeros([1,5]);
    
    %Initialize first values
    [MSFT_net,MSFT_YPred5d(1)] = predictAndUpdateState(MSFT_net,YTrain(end,1));
    [AAPL_net,AAPL_YPred5d(1)] = predictAndUpdateState(AAPL_net,YTrain(end,2));
    [AMZN_net,AMZN_YPred5d(1)] = predictAndUpdateState(AMZN_net,YTrain(end,3));
    [BRKB_net,BRKB_YPred5d(1)] = predictAndUpdateState(BRKB_net,YTrain(end,4));
    [FB_net,FB_YPred5d(1)] = predictAndUpdateState(FB_net,YTrain(end,5));
    
    
    for i=2:5 %five steps ahead
        [MSFT_net,MSFT_YPred5d(i)] = predictAndUpdateState(MSFT_net,MSFT_YPred5d(i-1));
        [AAPL_net,AAPL_YPred5d(i)] = predictAndUpdateState(AAPL_net,AAPL_YPred5d(i-1));
        [AMZN_net,AMZN_YPred5d(i)] = predictAndUpdateState(AMZN_net,AMZN_YPred5d(i-1));
        [BRKB_net,BRKB_YPred5d(i)] = predictAndUpdateState(BRKB_net,BRKB_YPred5d(i-1));
        [FB_net,FB_YPred5d(i)] = predictAndUpdateState(FB_net,FB_YPred5d(i-1));
    end
    
%% 2.a Create two tables
%%
%Table 1: Static Forecast Results

    StaticVariables = {'static_1day_ARMA','static_1day_LSTM','static_2day_ARMA','static_2day_LSTM', ...
        'static_3day_ARMA','static_3day_LSTM','static_4day_ARMA','static_4day_LSTM','static_5day_ARMA','static_5day_LSTM'};

    StaticForecasts_array = zeros([5,10]);
    
    for i=1:5      
        %ARMA static forecasts
        StaticForecasts_array(1,2*(i-1)+1) = MSFT_fcst5(i);
        StaticForecasts_array(2,2*(i-1)+1) = AAPL_fcst5(i);
        StaticForecasts_array(3,2*(i-1)+1) = AMZN_fcst5(i);
        StaticForecasts_array(4,2*(i-1)+1) = BRKB_fcst5(i);
        StaticForecasts_array(5,2*(i-1)+1) = FB_fcst5(i);
        
        %LSTM static forecasts (need to unstandardize before placing in table)
        StaticForecasts_array(1,2*i) = MSFT_YPred5(i)*sig(1)+mu(1);
        StaticForecasts_array(2,2*i) = AAPL_YPred5(i)*sig(2)+mu(2);
        StaticForecasts_array(3,2*i) = AMZN_YPred5(i)*sig(3)+mu(3);
        StaticForecasts_array(4,2*i) = BRKB_YPred5(i)*sig(4)+mu(4);
        StaticForecasts_array(5,2*i) = FB_YPred5(i)*sig(5)+mu(5);
    end
    
    StaticForecasts = array2table(StaticForecasts_array,'RowNames',assets,'VariableNames',StaticVariables)
    
%Table 2: Dynamic Forecast Results
    
    DynamicVariables = {'dyn_1day_ARMA','dyn_1day_LSTM','dyn_2day_ARMA','dyn_2day_LSTM', ...
        'dyn_3day_ARMA','dyn_3day_LSTM','dyn_4day_ARMA','dyn_4day_LSTM','dyn_5day_ARMA','dyn_5day_LSTM'};
    
    
    DynamicForecasts_array = zeros([5,10]);
    
    for i=1:5      
        %ARMA dynamic forecasts
        DynamicForecasts_array(1,2*(i-1)+1) = MSFT_fcst5d(i);
        DynamicForecasts_array(2,2*(i-1)+1) = AAPL_fcst5d(i);
        DynamicForecasts_array(3,2*(i-1)+1) = AMZN_fcst5d(i);
        DynamicForecasts_array(4,2*(i-1)+1) = BRKB_fcst5d(i);
        DynamicForecasts_array(5,2*(i-1)+1) = FB_fcst5d(i);
        
        %LSTM dynamic forecasts (need to unstandardize before placing in table)
        DynamicForecasts_array(1,2*i) = MSFT_YPred5d(i)*sig(1)+mu(1);
        DynamicForecasts_array(2,2*i) = AAPL_YPred5d(i)*sig(2)+mu(2);
        DynamicForecasts_array(3,2*i) = AMZN_YPred5d(i)*sig(3)+mu(3);
        DynamicForecasts_array(4,2*i) = BRKB_YPred5d(i)*sig(4)+mu(4);
        DynamicForecasts_array(5,2*i) = FB_YPred5d(i)*sig(5)+mu(5);
    end    
        
        
    DynamicForecasts = array2table(DynamicForecasts_array,'RowNames',assets,'VariableNames',DynamicVariables)
    
%% 2.b Comment on the two tables. Be sure to discuss if the static and dynamic forecasts look very different for either of the models.
%%

%For the ARMA models, the static and dynamic forecasts appear to be the
%same. For the LSTM models, the static and dynamic models are different but
%the discrepancy is small. For a five-day horizon, the choice between a
%static or dynamic model does not appear to make much of a difference. This
%is likely due to the magnitude of the compounding errors still being small in the dynamic model.

%% 2.c Calculate the RMSE for the static forecasts for each asset and model
%%

    %RMSE for static ARMA models
    MSFT_RMSEa = sqrt(MSFT_MSE(5));
    AAPL_RMSEa = sqrt(AAPL_MSE(5));
    AMZN_RMSEa = sqrt(AMZN_MSE(5));
    BRKB_RMSEa = sqrt(BRKB_MSE(5));
    FB_RMSEa = sqrt(FB_MSE(5));
    
    %Compute RMSE for static LSTM models
    YTest = dataTest(1:5,:);
    
    MSFT_RMSEl = sqrt(mean((MSFT_YPred5'-YTest(:,1)).^2));
    AAPL_RMSEl = sqrt(mean((AAPL_YPred5'-YTest(:,2)).^2));
    AMZN_RMSEl = sqrt(mean((AMZN_YPred5'-YTest(:,3)).^2));
    BRKB_RMSEl = sqrt(mean((BRKB_YPred5'-YTest(:,4)).^2));
    FB_RMSEl = sqrt(mean((FB_YPred5'-YTest(:,5)).^2));
    
    
%% 2.d Calculate the RMSE for the dynamic forecasts for each asset and model.
%%

    %RMSE for dynamic ARMA models
    MSFT_RMSEad = sqrt(MSFT_MSEd);
    AAPL_RMSEad = sqrt(AAPL_MSEd);
    AMZN_RMSEad = sqrt(AMZN_MSEd);
    BRKB_RMSEad = sqrt(BRKB_MSEd);
    FB_RMSEad = sqrt(FB_MSEd);
    
    %Compute RMSE for dynamic LSTM models
    YTest = dataTest(1:5,:);
    
    MSFT_RMSEld = sqrt(mean((MSFT_YPred5d'-YTest(:,1)).^2));
    AAPL_RMSEld = sqrt(mean((AAPL_YPred5d'-YTest(:,2)).^2));
    AMZN_RMSEld = sqrt(mean((AMZN_YPred5d'-YTest(:,3)).^2));
    BRKB_RMSEld = sqrt(mean((BRKB_YPred5d'-YTest(:,4)).^2));
    FB_RMSEld = sqrt(mean((FB_YPred5d'-YTest(:,5)).^2));
    
%% 2.e Create a table where you present the RMSE’s calculated above with asset names as the row labels.
%%

    RMSEForecasts_array = zeros([5,4]);
    
    %Static RMSE ARMA
    RMSEForecasts_array(1,1) = MSFT_RMSEa;
    RMSEForecasts_array(2,1) = AAPL_RMSEa;
    RMSEForecasts_array(3,1) = AMZN_RMSEa;
    RMSEForecasts_array(4,1) = BRKB_RMSEa;
    RMSEForecasts_array(5,1) = FB_RMSEa;
    
    %Static RMSE LSTM
    RMSEForecasts_array(1,2) = MSFT_RMSEl;
    RMSEForecasts_array(2,2) = AAPL_RMSEl;
    RMSEForecasts_array(3,2) = AMZN_RMSEl;
    RMSEForecasts_array(4,2) = BRKB_RMSEl;
    RMSEForecasts_array(5,2) = FB_RMSEl;    
    
    %Dynamic RMSE ARMA
    RMSEForecasts_array(1,3) = MSFT_RMSEad;
    RMSEForecasts_array(2,3) = AAPL_RMSEad;
    RMSEForecasts_array(3,3) = AMZN_RMSEad;
    RMSEForecasts_array(4,3) = BRKB_RMSEad;
    RMSEForecasts_array(5,3) = FB_RMSEad;
    
    %Dynamic RMSE LSTM
    RMSEForecasts_array(1,4) = MSFT_RMSEld;
    RMSEForecasts_array(2,4) = AAPL_RMSEld;
    RMSEForecasts_array(3,4) = AMZN_RMSEld;
    RMSEForecasts_array(4,4) = BRKB_RMSEld;
    RMSEForecasts_array(5,4) = FB_RMSEld;    
    
    RMSEForecasts = array2table(RMSEForecasts_array,'RowNames',assets,'VariableNames',{'static_RMSE_ARMA','static_RMSE_LSTM','dyn_RMSE_ARMA','dyn_RMSE_LSTM'})
    
    RMSE_means = mean(table2array(RMSEForecasts))
%% 2.f Comment on the above table. Be sure to discuss which model had the lowest RMSE on average, static vs dynamic accuracy, assets with high accuracy, etc...
%%
%The static LSTM model had the lowest RMSE on average, although its RMSE was not
%much smaller than that of the dynamic LSTM model.
% 
%There did not appear to be much difference between the RMSE for the
%static and dynamic models. For ARMA, the dynamic model was slightly more
%accurate. The Berkshire Hathaway Class B (BRK.B) and Microsoft (MSFT) assets generally had the
%highest accuracy. BRK.B predictions were the most accurate among the ARMA
%models, and MSFT predictions were the most accurate among the LSTM models.
%This is likely due to the low variability of returns near the end of the
%training period and beginning of the test period.