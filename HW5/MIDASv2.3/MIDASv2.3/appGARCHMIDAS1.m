clear
clc

% Purpose:
%
% This is a demo of GARCH-MIDAS model
%
% Reference:
%
% Engle,R.F., Ghysels, E. and Sohn, B. (2013), Stock Market Volatility 
% and Macroeconomic Fundamentals. The Review of Economics and Statistics,
% 95(3), 776-797.
%
% Written by Hang Qian
% Contact: matlabist@gmail.com

% NASDAQ Composite Index, daily percentage change 1971 - 2015
% Data Source: FRED Economic Data
% https://research.stlouisfed.org/fred2/series/NASDAQCOM
y = xlsread('NASDAQCOM.xls','B22:B11669') ./ 100;

% Estimate the GARCH-MIDAS model, and extract the volatilities
period = 22;
numLags = 24;
[estParams,EstParamCov,Variance,LongRunVar] = GarchMidas(y,'Period',period,'NumLags',numLags); %#ok<*ASGLU>

% Plot the conditional volatility and its long-run component
figure(1)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.1,2015.8,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1974,2016])
ylim([0,0.5])
hold off

% Estimate the rolling window version of the GARCH-MIDAS model
[estParams,EstParamCov,Variance,LongRunVar]...
    = GarchMidas(y,'Period',period,'NumLags',numLags,'RollWindow',1);

% Plot the conditional volatility and its long-run component
figure(2)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.1,2015.8,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1974,2016])
ylim([0,0.5])
hold off

% Industrial Production Index growth rate, 1971-2015
% Data Source: FRED database 
% https://research.stlouisfed.org/fred2/series/INDPRO
xMonth = xlsread('INDPRO.xls','B42:B576') ./ 100;

% Repeat the monthly value throughout the days in that month
[~,yDate] = xlsread('NASDAQCOM.xls','A22:A11669');
[~,yDateMonth] = datevec(yDate);
xDay = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)    
        count = count + 1;
        if count > length(xMonth)
            break
        end
    end
    xDay(t) = xMonth(count);
end

% Estimate the GARCH-MIDAS model with the exogenous regressor
[estParams,EstParamCov,Variance,LongRunVar] = GarchMidas(y,'Period',period,'NumLags',32,'X',xDay,'ThetaM',1);

% Plot the conditional volatility and its long-run component
figure(3)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.1,2015.8,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1974,2016])
ylim([0,0.5])
hold off

% In-sample forecast validation
GarchMidas(y,'Period',period,'NumLags',numLags,'estSample',8000);

% Out-of-sample forecast
estParams = GarchMidas(y,'Period',period,'NumLags',numLags);
nForecast = 5;
yBig = [y;0];
for t = 1:nForecast    
    [~,~,Variance,LongRunVar] = GarchMidas(yBig,'Period',period,'NumLags',numLags,'Params',estParams);
    yPseudo = estParams(1) + sqrt(Variance(end));
    yBig = [yBig(1:end-1);yPseudo;0];
end
VarianceForecast = Variance(nobs+1:nobs+nForecast);
LongRunVarForecast = LongRunVar(nobs+1:nobs+nForecast);


