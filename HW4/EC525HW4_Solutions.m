%THIS IS THE SOLUTIONS FOR HOMEWORK 4 SPRING 2019
%Housekeeping
    clear all; close all; clc
    
%% Establish Data Connection
%I am going to pull my data from Quandl
%Enter the API key
    apikey = 'MJ5soq3BhmRgiUyjF7aM'; 
%Establish database connection 
    c = quandl(apikey);
    
%reads excel sheet with S&P500 tickers and places them into array called tickers
tickers = ["MSFT" "AAPL" "AMZN" "BRK_B" "FB"];
%creates a vector that is the same size as the number of tickers with the string Wiki/ (quandl key)
tick(1:5) = "Wiki/";

tickers = tickers';
tick = tick';
%concatenates tick array to tickers array for each ticker to make the full quandl key (wiki/"ticker")
cattick = strcat(tick,tickers);   


%% Pull Raw Returns

%Create array for returns
%I pulled the first ticker separately because it was the easierst way to keep the dates in the first column
data1 = history(c,cattick(1),'12-31-2014','12-31-2017','daily',"column_index",11);
%creates a table to hold the adjusted close prices for each ticker over the time period
Ret = timetable2table(data1(:,'Adj0x2EClose'));
%changes the name of the adjusted close variable to the ticker name to avoid a repeating variable name error
Ret.Properties.VariableNames{'Adj0x2EClose'} = 'MSFT';    
    
%Pull the data for each ticker
    %for loop runs through each ticker
    for i = 2:5
        %pulls data from quandl, "column_index, 11" grabs the adjusted close price
        data1 = history(c,cattick(i,1),'12-31-2014','12-31-2017','daily',"column_index",11);
        %creates a temporary varibale to store the adjusted close price for each day and makes it a table 
        temp = timetable2table(data1(:,'Adj0x2EClose'));
        %creates another temporary variable to store a character vector of the ticker name that will be used to change the name of the variable
        tempName = char(tickers(i,1));
        %changes the name of the adjusted close variable to the ticker name to avoid a repeating variable name error 
        temp.Properties.VariableNames{'Adj0x2EClose'} = tempName;
            %if statement that checks if the data has values for every day, if not then it is not used (cleans data)
            if height(data1) == 755
                %If data for all days are available, it is concatenated to the full table of ticker prices (Ret)
                Ret = cat(2,Ret,temp(:,2));
            end
    end
%creates a copy of the prices that are now oldest to newest which is more desirable
Ret = flip(Ret);  


%% Clean Data (specifically AAPL and AMZN)
%%they had issues with data length

%Fix Apple's data
apple = timetable2table((history(c,"Wiki/AAPL",'12-31-2014','12-31-2017','daily',"column_index",11))); %Pull apples history 
apple2 = table2array(apple(:,2));
apple101 = array2table(["07-Aug-2017" 157.596]); %the average between the day before and day after was calculate to fill in the missing day
apple101.Properties.VariableNames = {'Time' 'Adj0x2EClose'};
applefinal = cat(1, apple(1:100,:), apple101);
apple = cat(1, applefinal, apple(101:end,:)); %Produces a final vector for apple that is the correct size with all the days filled in
apple = flip(apple);
apple.Properties.VariableNames{'Adj0x2EClose'} = 'AAPL';
Ret2 = cat(2,Ret(:,1:2), apple(:,2)); %puts apple into final returns array
 

%Fix amazon's data
amazon = timetable2table((history(c,"Wiki/AMZN",'12-31-2014','12-31-2017','daily',"column_index",11))); %pulls amazons hisotry
amazon2 = table2array(amazon(:,2));
amazon101 = array2table(["07-Aug-2017" 988.71]); %the average between the day before and day after was calculate to fill in the missing day
amazon101.Properties.VariableNames = {'Time' 'Adj0x2EClose'};
amazonfinal = cat(1,amazon(1:100,:), amazon101); 
amazon = cat(1,amazonfinal,amazon(101:end,:)); %Produces a final vector for apple that is the correct size with all the days filled in 
amazon = flip(amazon);
amazon.Properties.VariableNames{'Adj0x2EClose'} = 'AMZN';
%Place amazon into returns array
Ret3 = cat(2, Ret2(:,1:3),amazon(:,2));
Ret4 = cat(2,Ret3,Ret(:,'BRK_B'));
Ret = cat(2,Ret4,Ret(:,'FB')); %puts amazon into final returns array



%% Calculate Log Returns
Returns_all = table2array(Ret(:,2:end));
Returns_all = str2double(Returns_all);
Returns_all = diff(log(Returns_all));

%Create 90% training_part and testing_part period data
j = ceil(.9*size(Returns_all));
training_part = Returns_all(1:j,:);
testing_part = Returns_all(j+1:end,:);



%% 1A) ARMA(1,1) One Day Ahead Forecast Static
mdl = arima(1,0,1); %setup ARMA(1,1) model

%% 1C) ARMA(1,1) Static forecasts..

five_Static= zeros(5,5);
for i = 1:5
    est = estimate(mdl, training_part(:,i));
    five_Static(i,:) = forecast(est,5,'Y0',training_part(:,i));   
end



%% 1E) Dynamic ARMA(1,1) 

dyn_frcst = zeros(5,5);

for i = 1:5 %loop through each asset
    
   dyn_data = training_part(:,i); %take the asset's data from the training_part data
   estDyn = estimate(mdl, dyn_data); %estimate based on asset's data
   
   for k = 1:5 %loop through each day
      temp = forecast(estDyn,k,'Y0',dyn_data); %temp variable to store forecast data
      dyn_frcst(i,k) = temp(k,1); %add data to forecast array
      dyn_data = cat(1, dyn_data, temp(k,1)); %add forecasted data into array to use in the next forecast
      estDyn = estimate(mdl,dyn_data); %recursively forecast the next day with the new array
   end    
end    



%% 1B&D) LSTM
finalYPred = [];


for i = 1:5
%Standardize the training data
mu = mean(training_part(:,i));
sig = std(training_part(:,i));
dataTrainStandardized(:,i) = (training_part(:,i) - mu) ./ sig;


%Lag x for prediction
XTrain = dataTrainStandardized(1:end-1,i); 
YTrain = dataTrainStandardized(2:end,i);


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
  'MaxEpochs',2, ...
  'GradientThreshold',1, ...
  'InitialLearnRate',0.005, ...
  'LearnRateSchedule','piecewise', ...
  'LearnRateDropPeriod',125, ...
  'LearnRateDropFactor',0.2, ...
  'Verbose',0);
  %'Plots','training-progress');
  
  
%Train the network
net = trainNetwork(XTrain',YTrain',layers,options);


%Forecast
%Standardize
dat_test_std = (testing_part(:,i) - mu) ./ sig;

%forecast using predict, could have used other function
one_LSTM = predict(net,dat_test_std(end,:));%one day ahead
two_LSTM = predict(net,dat_test_std(end-1,:));
three_LSTM = predict(net,dat_test_std(end-2,:));
four_LSTM = predict(net,dat_test_std(end-3,:));
five_LSTM = predict(net,dat_test_std(end-4,:));


%unstanderdize forecasts
one = one_LSTM*sig + mu;
two = two_LSTM*sig + mu;
three = three_LSTM*sig + mu;
four = four_LSTM*sig + mu;
five = five_LSTM*sig + mu;
all = [one;two;three;four;five];
finalYPred = cat(2, finalYPred, all);
end


finalYPred = finalYPred';




%% 1F) LSTM Dynamic

dynamicPredict = zeros(5,5);
dynamicPredict(:,1) = finalYPred(:,1);


for i = 1:5
    %Standardize the training data
    mu = mean(training_part(:,i));
    sig = std(training_part(:,i));
    dataTrainStandardized(:,i) = (training_part(:,i) - mu) ./ sig;
    %Use the time t to predict t+1
    XTrain = dataTrainStandardized(1:end-1,i); %training of the training
    YTrain = dataTrainStandardized(2:end,i);
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
      'MaxEpochs',2, ...
      'GradientThreshold',1, ...
      'InitialLearnRate',0.005, ...
      'LearnRateSchedule','piecewise', ...
      'LearnRateDropPeriod',125, ...
      'LearnRateDropFactor',0.2, ...
      'Verbose',0);
      %'Plots','training-progress');
    %Train the network
    net = trainNetwork(XTrain',YTrain',layers,options);
    %Forecast
    %Standardize
    
    net = predictAndUpdateState(net,XTrain');
    [net,YPred] = predictAndUpdateState(net,YTrain(end));

    for i = 2:5
        [net,YPred(i,1)] = predictAndUpdateState(net,YPred(i-1,1),'ExecutionEnvironment','cpu');
    end
    
     
       %Unstandardize the prediction for y
    dynamicPredict(:,j) = sig(j)*YPred + mu(j);
       
     
end
   

%% 2C&D) Calculate the RMSE

%grab testing data....
testing_partY = testing_part(1:5,:);
testing_partY = testing_partY.';


LSTM_statrmse = sqrt(mean((testing_partY - finalYPred).^2,2));
dyn_rmse = sqrt(mean((testing_partY - dynamicPredict).^2,2));
five_rmse = sqrt(mean((testing_partY - five_Static).^2,2));
dyn_rmsefcst = sqrt(mean((testing_partY - dyn_frcst).^2,2));



%% 2A&E) Create final output tables

%create names of columns and rows
rowNames = ["Microsoft", "Apple", "Amazon", "Berkshire", "Facebook"];
colNames = ["Static1dayARMA","Static1dayLSTM","Static2dayARMA","Static2dayLSTM","Static3dayARMA","Static3dayLSTM","Static4dayARMA","Static4dayLSTM","Static5dayARMA","Static5dayLSTM",];
colNames2 = ["Dyn1dayARMA","Dyn1dayLSTM","Dyn2dayARMA","Dyn2dayLSTM","Dyn3dayARMA","Dyn3dayLSTM","Dyn4dayARMA","Dyn4dayLSTM","Dyn5dayARMA","Dyn5dayLSTM",];
RMSEnames = ["StaticRMSE_ARMA","StaticRMSE_LSTM","DynamicRMSE_ARMA","DynamicRMSE_LSTM"];

StaticFrcst_final = table(five_Static(:,1),finalYPred(:,1),five_Static(:,2),finalYPred(:,2),five_Static(:,3),finalYPred(:,3),five_Static(:,4),finalYPred(:,4),five_Static(:,5),finalYPred(:,5), 'VariableNames', colNames,'RowNames',rowNames);
dyn_frcstResults = table(dyn_frcst(:,1),dynamicPredict(:,1),dyn_frcst(:,2),dynamicPredict(:,2),dyn_frcst(:,3),dynamicPredict(:,3),dyn_frcst(:,4),dynamicPredict(:,4),dyn_frcst(:,5),dynamicPredict(:,5),'VariableNames', colNames2,'RowNames',rowNames);
RMSEForecastResults = table(five_rmse,LSTM_statrmse,dyn_rmsefcst,dyn_rmse,'VariableNames', RMSEnames,'RowNames',rowNames);

