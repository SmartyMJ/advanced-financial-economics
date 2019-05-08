function Yf = ForecastCombine(varargin)

%ForecastCombine Forecast combination of mixed frequency models
%
% Syntax:
%
%  Yf = ForecastCombine(OutputForecastAIO)
%  Yf = ForecastCombine(OutputForecastAIO, weightFlag)
%  Yf = ForecastCombine(OutputForecast1, OutputForecast2, ...)
%  Yf = ForecastCombine(OutputForecast1, OutputForecast2, ...,weightFlag)
%
% Description:
%
%   Forecast of different models/data are combined by discounted or plain
%   mean squared error of forecast.
%
% Input Arguments:
%
%   OutputForecastAIO - a cell array of OutputForecast struct using
%                       different models/data, as created by MIDAS_ADL
%
%   weightFlag        - a string that specifies the combination criterion
%      o 'MSFE':  weights by mean squared error of forecast (default)
%      o 'DMSFE': weights by discounted mean squared error of forecast
%      o 'aic':   weights by Akaike information criteria of the regression
%      o 'bic':   weights by Bayesian information criteria of the regression
%      o 'Flat':  equal weights
%
% Output Arguments:
%
%   Yf - weighted average forecast
%
% Notes:
%
% o The forecasting periods must be identical for all models.
%
% Contact: matlabist@gmail.com

weightFlag = 'MSFE';

% Parse inputs
OutputForecastAIO = varargin{1};
if isstruct(OutputForecastAIO)
    % syntax: Yf = ForecastCombine(OutputForecast1, OutputForecast2, ...,weightFlag)
    nmodel = length(varargin);
    OutputForecastAIO = cell(nmodel,1);
    for m = 1:nmodel
        if isstruct(varargin{m})
            OutputForecastAIO{m} = varargin{m};
        end
        if ischar(varargin{m})
            weightFlag = varargin{m};
            nmodel = m - 1;
            break
        end
    end
elseif iscell(OutputForecastAIO)
    % syntax: Yf = ForecastCombine(OutputForecastAIO, weightFlag)
    nmodel = length(OutputForecastAIO);
    if length(varargin) == 2
        weightFlag = varargin{2};
    end
else
    error('Unrecognized input arugments.')    
end

% Fetch forecast values and weights for models
nperiod = length(OutputForecastAIO{1}.Yf);
YfAIO = zeros(nperiod,nmodel);
weight = zeros(1,nmodel);
for m = 1:nmodel
    
    if length(OutputForecastAIO{m}.Yf) ~= nperiod
        error('Forecasting periods must be identical for all models')
    end
    YfAIO(:,m) = OutputForecastAIO{m}.Yf;
    
    switch upper(weightFlag)
        case 'MSFE'
            weight(m) = 1 / OutputForecastAIO{m}.MSFE;
        case 'DMSFE'
            weight(m) = 1 / OutputForecastAIO{m}.DMSFE; 
        case 'AIC'
            weight(m) = exp(-OutputForecastAIO{m}.aic);    
        case 'BIC'
            weight(m) = exp(-OutputForecastAIO{m}.bic);
        case 'FLAT'
            weight(m) = 1;
    end    
end

% Forecast combination
weight = weight ./ sum(weight);
YfAIOWight = bsxfun(@times,YfAIO,weight);
Yf = sum(YfAIOWight,2);        
            
            
    
 