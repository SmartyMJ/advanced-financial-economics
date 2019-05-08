function [estParams,condQuantile,yLowFreq,xHighFreq,yDates] = MidasQuantile(y,varargin)
%MidasQuantile: MIDAS quantile regression
%
% Syntax:
%
%   [estParams,condQuantile] = MidasQuantile(y)
%   [estParams,condQuantile,yLowFreq,xHighFreq,yDates] = MidasQuantile(y, name,value,...)
%
% Description:
%
%  Suppose that we have a return series y(t), t = 1,...,T.
%  MIDAS quantile regression estimates the conditional quantile of n-period
%  returns (obtained by aggregating y(t),...,y(t+n)). The conditioning
%  variable (predictor) is sampled at high frequency with MIDAS weights.
%  The default predictor is |y(t)|, as "absolute returns successfully 
%  capture time variation in the conditional distribution of returns".
%
%  Refer to Ghysels, Plazzi and Valkanov (2016) for model specification.
%
% Input Arguments:
%
%   y           T-by-1 observation data
%
% Optional Input Name/Value Pairs:
%
%  'Quantile'  A scalar between zero and one that specifies the level
%              (i.e., alpha) of quantile. The default is 0.05.
%
%  'X'         T-by-1 high frequency conditioning variable (predictor). 
%              This variable must have the same length as y, and only one
%              predictor is supported. By default, it is |y|, as "absolute
%              returns successfully capture time variation in the 
%              conditional distribution of returns".
%
%  'Period'    A scalar integer that specifies the aggregation periodicity.
%              y will be aggregated so as to formulate n-period returns.
%              How many days in a week/month/quarter/year?
%              The default is 22 (as in a day-month aggregation)
%
%  'NumLags'   A scalar integer that specifies the number of lags for the
%              high frequency predictor, to which MIDAS weights is
%              assigned. The default is 250.
%
%  'Dates'     T-by-1 vector or cell array for the dates of y. This is
%              for book-keeping purpose and does not affect estimation.
%              The default is 1:length(y).
%
%  'Smoother'  A non-negative scalar that specifies how to smooth the
%              non-differentiable objective function. If it is zero, there
%              is no smoothing. The default is average absolute residuals. 
%              This is the starting smoother. The software will run a
%              series of optimizations; each time the smoother will be
%              reduced by an half.
%
%  'Search'    A logical value that indicates numerical minimization
%              via pattern search in the Global Optimization Toolbox.
%              If not available, it resorts to fminsearch in base MATLAB.
%              The default is false (and will use gradient-based methods
%              under smoothed objective functions)
%
%  'Options'   The options for numerical optimization. 
%              The default is the FMINCON default choice.
%
%  'Gradient'  A logical value that indicates analytic gradients.
%              The default is false.
%
%  'Bootstrap' A character vector that specifies bootstrap standard error
%              method: 'Residual' (Default) or 'XY'.
%
%  'Params'    Parameter values for [intercept;slope;k].
%              In that case, the program skips estimation, and just infers
%              conditional quantiles based on the specified parameters.
%              The default is empty (need parameter estimation).
%
%  'Params0'   Starting parameter values of [intercept;slope;k] for numeric
%              optimization. This will overload the software default
%              choice, which starts from the OLS estimator.
%
% Output Arguments:
%
%   estParams   Estimated parameters for [intercept;slope;k],
%               where intercept and slope are the coefficients of the
%               quantile regression, and k is the parameter in the
%               MIDAS Beta polynomial
%
%   condQuantile R-by-1 conditional quantile. This is the fitted value of
%                the right-hand-side of the quantile regression.
%                R = T - Period - NumLags + 1.
%
%   yLowFreq     R-by-1 n-period returns obtained from overlapping
%                aggregation of y. This is the left-hand-side of the
%                quantile regression.
%
%   xHighFreq    R-by-NumLags high-frequency predictor data used by the
%                quantile regression
%
%   yDates       R-by-1 serial dates for the output variables.
%
%
% Notes:
%
% o If 'Smoother' = 0, the software will numerically optimize a
%   non-differentiable objective function. The codes will run faster. 
%   Sometimes it is necessary to fine-tune 'Smoother' for successful
%   optimization.
%
% o If numerical optimization works poorly, try 'Gradient' = 1.
%   The codes will be slower for MATLAB versions earlier than R2015b.
%
% Reference:
%
% Ghysels, E., Plazzi, A. and Valkanov, R. (2016) Why Invest in Emerging
% Market? The Role of Conditional Return Asymmetry. Journal of Finance,
% forthcoming.


% Parse inputs and set defaults
quantileDefault = 0.05;
periodDefault = 22;
nlagDefault = 250;
callerName = 'MidasQuantile';
parseObj = inputParser;
addParameter(parseObj,'Quantile',quantileDefault,@(x)validateattributes(x,{'numeric'},{'scalar','>',0,'<',1},callerName));
addParameter(parseObj,'X',[],@(x)validateattributes(x,{'numeric'},{'2d'},callerName));
addParameter(parseObj,'Period',periodDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'NumLags',nlagDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'Dates',[],@(x)validateattributes(x,{'numeric','cell'},{},callerName));
addParameter(parseObj,'Smoother',[],@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},callerName));
addParameter(parseObj,'Search',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Options',[],@(x)validateattributes(x,{},{},callerName));
addParameter(parseObj,'Gradient',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Bootstrap','Residual',@(x)validateattributes(x,{'char'},{},callerName));
addParameter(parseObj,'Params',[],@(x)validateattributes(x,{'numeric'},{'column'},callerName));
addParameter(parseObj,'Params0',[],@(x)validateattributes(x,{'numeric'},{'column'},callerName));
parse(parseObj,varargin{:});
q = parseObj.Results.Quantile;
Regressor = parseObj.Results.X;
period = parseObj.Results.Period;
nlag = parseObj.Results.NumLags;
yDates = parseObj.Results.Dates;
smoother = parseObj.Results.Smoother;
searchFlag = parseObj.Results.Search;
options = parseObj.Results.Options;
autodiffFlag = parseObj.Results.Gradient;
Bootstrap = parseObj.Results.Bootstrap;
estParams = parseObj.Results.Params;
params0 = parseObj.Results.Params0;

% Replace missing values by the sample average
y = y(:);
y(isnan(y)) = nanmean(y);
nobs = length(y);

% Load the conditioning variable (predictor)
if isempty(Regressor)
    Regressor = abs(y);    
else
    if numel(Regressor) ~= numel(y)
        error('Conditioning variable (predictor) must be a vector of the same length as y.')
    end
    Regressor = Regressor(:);
    Regressor(isnan(Regressor)) = nanmean(Regressor);
end

% Load dates
if isempty(yDates)
    yDates = (1:nobs)';
elseif iscell(yDates)
    yDates = datenum(yDates);    
end

if numel(yDates) ~= nobs
    error('Length of Dates must equal the number of observations.')
end

% Prepare data for the LHS and RHS of the quantile regression
% LHS: n-period returns by aggregating y(t),...,y(t+period-1)
% RHS: many lagged returns by extracting y(t-1),...,y(t-nlag)
nobsShort = nobs-nlag-period+1;
yLowFreq = zeros(nobsShort,1);
xHighFreq = zeros(nobsShort,nlag);
for t = nlag+1 : nobs-period+1
    yLowFreq(t-nlag,1) = sum(y(t:t+period-1));  
    xHighFreq(t-nlag,:) = Regressor(t-1:-1:t-nlag)';
end
yDates = yDates(nlag+1:nobs-period+1);

% In case of known parameters, just compute conditional quantile and exit.
if ~isempty(estParams)
    if numel(estParams) ~= 3
        error('The length of Params must be 3.')
    end
    [~,condQuantile] = objFun(estParams,yLowFreq,xHighFreq,q,0);
    return
end

% Initial parameters by OLS
if isempty(params0)
    k0 = 5;
    X0 = [ones(nobsShort,1),xHighFreq * midasBetaWeights(nlag,1,k0)'];
    OLS = X0 \ yLowFreq;
    params0 = [OLS;k0];    
    resid = yLowFreq - X0 * OLS;
else
    k0 = params0(3);
    X0 = [ones(nobsShort,1),xHighFreq * midasBetaWeights(nlag,1,k0)'];
    resid = yLowFreq - X0 * params0(1:2);    
end

if isempty(smoother)
    smoother = mean(abs(resid));
end

% Bounds for numerical optimization
lb = [-Inf;-Inf;0];
ub = [Inf;Inf;100];

% Optimization options
if isempty(options)
    if autodiffFlag
        % options = optimoptions('fmincon','GradObj','on','DerivativeCheck','on','FinDiffType','central','Algorithm','interior-point','Display','notify-detailed');
        options = optimoptions('fmincon','GradObj','on','Algorithm','interior-point','Display','off');
    else
        options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
    end        
end

% Numeric minimization
if smoother == 0 && ~searchFlag && ~autodiffFlag
    % Minimize non-differentiable function by Newton method, finite-difference
    estParams = fmincon(@(params) objFun(params,yLowFreq,xHighFreq,q,0),params0,[],[],[],[],lb,ub,[],options);    
elseif smoother == 0 && ~searchFlag && autodiffFlag
    % Minimize non-differentiable function by Newton method, analytic grad
    estParams = fmincon(@(params) objFunGrad(params,yLowFreq,xHighFreq,q,0),params0,[],[],[],[],lb,ub,[],options);    
elseif searchFlag
    if exist('patternsearch','file') ~= 0
        % Minimize non-differentiable function by pattern search
        estParams = patternsearch(@(params) objFun(params,yLowFreq,xHighFreq,q,0),params0,[],[],[],[],lb,ub);        
    else
        % Minimize non-differentiable function by Nelder-Mead search
        estParams = fminsearch(@(params) objFun(params,yLowFreq,xHighFreq,q,0),params0);        
    end
else
    % Minimize a series of smooth functions with decreasing intensities
    % In each round, smoother is reduced by an half
    maxIter = 10;
    estParamsCand = zeros(3,maxIter);
    fvalCand = Inf(1,maxIter);    
    smootherUse = smoother;
    for r = 1:maxIter
       estParamsCand(:,r) = fmincon(@(params) objFunGrad(params,yLowFreq,xHighFreq,q,smootherUse),params0,[],[],[],[],lb,ub,[],options);
       fvalCand(r) = objFun(estParamsCand(:,r),yLowFreq,xHighFreq,q,0); % Evaluate the original objective function
       if r>1 && abs((fvalCand(r)-fvalCand(r-1))/fvalCand(r)) < 1e-4
           break
       end
       smootherUse = smootherUse / 2;
       params0 = estParamsCand(:,r);
    end    
    [~,minInd] = min(fvalCand);
    estParams = estParamsCand(:,minInd);    
end
[fval,condQuantile] = objFun(estParams,yLowFreq,xHighFreq,q,0);

% Bootstrap standard errors
currentSeed = rng;
rng(12345);
nsim = 100;
resid = yLowFreq - condQuantile;
ResidualBootstrap = strncmpi(Bootstrap,'Residual',1);
paramSim = zeros(3,nsim);
for r = 1:nsim
    ind = randi(nobsShort,[nobsShort,1]);
    if ResidualBootstrap
        yLowFreqSim = condQuantile + resid(ind);
        paramSim(:,r) = fminsearch(@(params) objFun(params,yLowFreqSim,xHighFreq,q,0),estParams);
    else
        paramSim(:,r) = fminsearch(@(params) objFun(params,yLowFreq(ind),xHighFreq(ind,:),q,0),estParams);
    end
end
se = std(paramSim,0,2);
zstat = estParams ./ se;
pval = 0.5 * erfc(0.7071 * abs(zstat)) * 2;
pval(pval<1e-6) = 0;
rng(currentSeed);

% Display the estimation results
if smoother == 0 && ~searchFlag && ~autodiffFlag
    fprintf('Method: Asymmetric loss function minimization\n');
elseif smoother == 0 && ~searchFlag && autodiffFlag
    fprintf('Method: Asymmetric loss function minimization, Analytic gradient Newton iterations\n');
elseif searchFlag
    if exist('patternsearch','file') ~= 0        
        fprintf('Method: Asymmetric loss function minimization, Pattern search\n');
    else        
        fprintf('Method: Asymmetric loss function minimization, Nelder-Mead search\n');
    end
else
    if autodiffFlag
        fprintf('Method: Smoothed asymmetric loss function minimization, Analytic gradient Newton iterations\n');
    else
        fprintf('Method: Smoothed asymmetric loss function minimization\n');
    end    
end
fprintf('Sample size:                 %d\n',nobs);
fprintf('Adjusted sample size:        %d\n',nobsShort);
fprintf('Minimized function value: %10.6g\n',fval);
columnNames = {'Coeff','StdErr','tStat','Prob'};
rowNames = {'Intercept';'Slope';'k2'};
try
    Table = table(estParams,se,zstat,pval,'RowNames',rowNames,'VariableNames',columnNames);
catch
    Table = [estParams,se,zstat,pval];
end
disp(Table)

%{
% plot quantiles
figure(3)
plot(yDates, condQuantile)
dateaxis
%}

end


%-------------------------------------------------------------------------
% Local function: the objective function
function [fval,condQuantile] = objFun(params,y,X,q,smoother)

% Allocate parameters
intercept = params(1);
slope = params(2);
k2 = params(3);

% Compute MIDAS weights
nlag = size(X,2);
k1 = 1;
weights = midasBetaWeights(nlag,k1,k2)';

% Conditional quantile
condQuantile = intercept + slope .* (X * weights);

% Asymmetric loss function
loss = y - condQuantile;
if smoother == 0
    % Non-differentiable loss function    
    fval = loss' * (q - (loss<0));
else
    % Piecewise linear-quadratic smoothing loss function
    maskSmall = loss < (q-1)*smoother;
    maskBig = loss > q*smoother;
    maskMedian = ~maskSmall & ~maskBig;
    fvalSmall = -0.5*(q-1)^2*smoother + (q-1)*loss(maskSmall);
    fvalBig = -0.5*q^2*smoother + q*loss(maskBig);
    fvalMedian = 0.5/smoother .* loss(maskMedian).^2;
    fval = sum(fvalSmall) + sum(fvalBig) + sum(fvalMedian);    
end

end


%-------------------------------------------------------------------------
% Local function: the objective function
function [fval,condQuantile] = objFunCapital(params,y,X,q,smoother)

% Allocate parameters
intercept = params(1);
slope = params(2);
k2 = params(3);

% Compute MIDAS weights
nlag = size(X,2);
k1 = 1;
weights = midasBetaWeightsCapital(nlag,k1,k2).';

% Conditional quantile
condQuantile = intercept + TIMES(slope, MTIMES(X, weights));

% Asymmetric loss function
loss = y - condQuantile;
if smoother == 0
    % Non-differentiable loss function    
    fval = MTIMES(loss.', (q - (loss<0)));
else
    % Piecewise linear-quadratic smoothing loss function
    maskSmall = loss < (q-1)*smoother;
    maskBig = loss > q*smoother;
    maskMedian = ~maskSmall & ~maskBig;
    fvalSmall = -0.5*(q-1)^2*smoother + (q-1)*loss(maskSmall);
    fvalBig = -0.5*q^2*smoother + q*loss(maskBig);
    fvalMedian = 0.5/smoother .* TIMES(loss(maskMedian),loss(maskMedian));
    fval = sum(fvalSmall) + sum(fvalBig) + sum(fvalMedian);    
end

end


%-------------------------------------------------------------------------
% Local function: Compute objective function value and gradient
function [fval,Gradient] = objFunGrad(params,y,X,q,smoother)

fval = objFun(params,y,X,q,smoother);
if nargout == 1    
    return
end

fun = @(params) objFunCapital(params,y,X,q,smoother);
x = params;

nfunVal = numel(fval);
nparams = numel(x);
fvalVec = zeros(nfunVal,nparams);
Gradient = zeros(nfunVal,nparams);

for m = 1:nparams
    
    % One element of x carries an imaginary number
    xComplex = x;
    xComplex(m) = xComplex(m) + 1i;
    
    % The real component is the function value 
    % The imaginary component is the derivative
    complexValue = fun(xComplex);
    fvalVec(:,m) = real(complexValue);
    Gradient(:,m) = imag(complexValue);
end

if any(norm(bsxfun(@minus,fvalVec, fval))>1e-6)
    warning('Gradients might be falsely computed. Capitalize arithmetic operator and utility function names.')
end

end

%-------------------------------------------------------------------------
% Local function: MIDAS beta polynomial weights
function weights = midasBetaWeights(nlag,param1,param2)
seq = linspace(eps,1-eps,nlag);
if param1 == 1    
    weights = (1-seq).^(param2-1);    
else
    weights = (1-seq).^(param2-1) .* seq.^(param1-1);    
end
weights = weights ./ nansum(weights);
end


%-------------------------------------------------------------------------
% Local function: MIDAS beta polynomial weights
function weights = midasBetaWeightsCapital(nlag,param1,param2)
seq = linspace(eps,1-eps,nlag);
if param1 == 1
    weights = POWER(1-seq, param2-1);
else
    weights = TIMES( POWER(1-seq, param2-1), POWER(seq, param1-1)) ;
end
weights = RDIVIDE(weights,nansum(weights));
end


%-------------------------------------------------------------------------
function B = EXP(A)
Areal = real(A);
Aimag = imag(A);
Breal = exp(Areal);
Bimag = Breal .* Aimag;
B = complex(Breal, Bimag);
end


%-------------------------------------------------------------------------
function B = LOG(A)
Areal = real(A);
Aimag = imag(A);
Breal = log(Areal);
Bimag = Aimag ./ Areal;
B = complex(Breal, Bimag);
end


%-------------------------------------------------------------------------
function C = MTIMES(A,B)
C = A * B + imag(A) * imag(B);
end


%-------------------------------------------------------------------------
function C = POWER(A,B)
if isscalar(B) && (B == 3)
    C = TIMES(TIMES(A, A), A);
elseif isscalar(B) && (B == 2)
    C = TIMES(A, A);
elseif isscalar(B) && (B == 1)
    C = A;
elseif isscalar(B) && (B == 0.5)
    C = SQRT(A);
elseif isscalar(B) && (B == 0)
    C = complex(ones(size(A)),0);
elseif isscalar(B) && (B == -0.5)
    C = RDIVIDE(1, SQRT(A));
elseif isscalar(B) && (B == -1)
    C = RDIVIDE(1, A);
elseif isscalar(B) && (B == -2)
    C = RDIVIDE(1, TIMES(A,A));
elseif all(A(:)>0)
    C = EXP(TIMES(B, LOG(A)));
else
    % This is not necessarily correct
    Creal = real(A).^real(B);
    Cimag = imag(exp(B.*log(abs(A))));
    C = complex(Creal,Cimag);
end
end


%-------------------------------------------------------------------------
function C = RDIVIDE(A,B)
Areal = real(A);
Aimag = imag(A);
Breal = real(B);
Bimag = imag(B);
Creal = Areal ./ Breal;
Cimag = (Aimag .* Breal - Areal .* Bimag) ./ (Breal.*Breal);
C = complex(Creal, Cimag);
end


%-------------------------------------------------------------------------
function C = TIMES(A,B)
C = A .* B + imag(A) .* imag(B);
end
