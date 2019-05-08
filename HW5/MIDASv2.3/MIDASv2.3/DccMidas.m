function [estParamsStep1,EstParamCovStep1,estParamsStep2,EstParamCovStep2,Variance,LongRunVar,CorrMatrix,LongRunCorrMatrix,logL] = DccMidas(Data,varargin)
%DccMidas: Maximum likelihood estimation of DCC-MIDAS
%
% Syntax:
%
%   [...] = DccMidas(Data)
%   [...] = DccMidas(Data, name,value,...)
%
% Description:
%
%   The DCC-MIDAS model decomposes the conditional covariance matrix into
%   the variances and the correlation matrix, with a two-step model
%   specification and estimation strategy. In the first step, conditional
%   variances are obtained from the univariate GARCH-MIDAS models with  the
%   long-run component being the realized volatility. In the second step,
%   observations are deflated by the estimated mean and conditional
%   variances, and the standardized residuals are thus constructed. The
%   standardized residuals have a correlation matrix with GARCH-MIDAS-like
%   dynamics. The long-run component is determined by the recent history of
%   sample autocorrelations under MIDAS weights.
%
%   Refer to Colacito, Engle and Ghysels (2011) for model specification.
%
% Input Arguments:
%
%   Data           T-by-n multivariate observation data
%
% Optional Input Name/Value Pairs:
%
%   'Period'    A scalar integer that specifies the aggregation periodicity
%               How many days in a week/month/quarter/year?
%               How long is the secular component fixed?
%               This value is for both realized volatility aggregation
%               (first step estimation) and historical sample correlation
%               matrix aggregation (second step estimation)
%               The default is 22 (as in a day-month aggregation)
%
%   'NumLagsVar' A scalar integer that specifies the number of lags in
%               filtering the secular component by MIDAS weights.
%               This is for the first step GARCH-MIDAS model.
%               The default is 10 (say a history of 10 weeks/months/quarters/years)
%
%   'NumLagsCorr' A scalar integer that specifies the number of lags in
%               filtering the secular component by MIDAS weights.
%               This is for the second step estimation of correlation matrix
%               The default is 10 (say a history of 10 weeks/months/quarters/years)
%
%   'EstSample' A scalar integer that specifies estimation sample size.
%               That is, Data(1:EstSample,:) are used for parameter
%               estimation. The remaining sample is for conditional               
%               variance forecast.
%               The default is length(y), no forecast
%
%   'RollWindow' A logical value that indicates rolling window estimation
%               on the long-run component. If true, the long-run component
%               varies every period. If false, the long-run component will
%               be fixed for a week/month/quarter/year.
%               This is for both two steps of estimation
%               The default is false            
%
%   'LogTau'    A logical value that indicates logarithmic long-run
%               volatility component. 
%               This is for the first step GARCH-MIDAS model
%               The default is false
%
%   'Beta2Para' A logical value that indicates two-parameter Beta MIDAS polynomial
%               This is for the both steps of estimation of DCC-MIDAS model
%               The default is false (one-parameter Beta polynomial)
%
%   'Options'   The FMINCON options for numerical optimization. For example,
%               Display iterations: optimoptions('fmincon','Display','Iter');
%               Change solver: optimoptions('fmincon','Algorithm','active-set');
%               The default is the FMINCON default choice
%
%   'Mu0'       MLE starting value for the location-parameter mu
%               This is for the first step GARCH-MIDAS model
%               The default is the sample average of observations
%
%   'Alpha0'    MLE starting value for alpha in the short-run GARCH(1,1) component
%               This is for the first step GARCH-MIDAS model
%               The default is 0.05 * ones(n,1)
%
%   'Beta0'     MLE starting value for beta in the short-run GARCH(1,1) component
%               This is for the first step GARCH-MIDAS model
%               The default is 0.9 * ones(n,1)
%
%   'Theta0'    MLE starting value for the MIDAS coefficient sqrt(theta) in the long-run component
%               This is for the first step GARCH-MIDAS model
%               The default is 0.1 * ones(n,1)
%
%   'W0'        MLE starting value for the MIDAS parameter w in the long-run component
%               This is for the first step GARCH-MIDAS model
%               The default is 5 * ones(n,1)
%
%   'M0'        MLE starting value for the location-parameter sqrt(m) in the long-run component
%               This is for the first step GARCH-MIDAS model
%               The default is 0.01 * ones(n,1)
%
%   'CorrA0'    MLE starting value for a in the GARCH(1,1) component
%               It is either a scalar (if all variables share it) or a
%               column vector (if each variable has its own parameter).
%               This is for the second step correlation matrix estimation
%               The default is 0.05 (or a vector expansion)
%
%   'CorrB0'    MLE starting value for b in the GARCH(1,1) component
%               It is either a scalar (if all variables share it) or a
%               column vector (if each variable has its own parameter).
%               This is for the second step correlation matrix estimation
%               The default is 0.9 (or a vector expansion)
%
%   'CorrW0'    MLE starting value for the MIDAS parameter w in the long-run component
%               It is either a scalar (if all variables share it) or a
%               column vector (if each variable has its own parameter).
%               This is for the second step correlation matrix estimation
%               The default is 5 (or a vector expansion)
%
%   'MorePara'  A logical value that indicates multivariate series have
%               different a,b. However, the program only supports a single w.
%               This is for the second step correlation matrix estimation
%               The default is false (parameters a,b,w are shared by all variables)
%
%   'Gradient'  A logical value that indicates analytic gradients in MLE
%               The default is false
%
%   'AdjustLag' A logical value that indicates MIDAS lag adjustments for
%               initial observations due to missing presample values
%               The default is false
%
%   'ThetaM'    A logical value that indicates not taking squares for the
%               parameter theta and m in the long-run volatility component
%               This is for the first step GARCH-MIDAS model
%               The default is false (they are squared)
%
%   'ZeroLogLStep1' A vector of indices between 1 and T, which select some dates
%              and forcefully reset likelihood values of those dates to zero.
%              For example, use ZeroLogL to ignore initial likelihood values. 
%              This is for the first step GARCH-MIDAS model
%              The default is empty (no reset)
%
%   'ZeroLogL' A vector of indices between 1 and T, which select some dates
%              and forcefully reset likelihood values of those dates to zero.
%              For example, use ZeroLogL to ignore initial likelihood values. 
%              This is for the second step DCC-MIDAS model
%              The default is empty (no reset)
%
% Output Argument:
%
%   estParamsStep1   6-by-n estimated parameters for [mu;alpha;beta;theta;w;m]
%                    obtained from univariate GARCH-MIDAS models
%
%   EstParamCovStep1 6-by-6-by-n estimated parameter covariance matrix,
%                    obtained from univariate GARCH-MIDAS models
%
%   estParamsStep2   3-by-1 or (2n+1)-by-1 estimated parameters,
%                    obtained from the second-step autocorrelation matrix estimation
%
%   EstParamCovStep2 3-by-3 or (2n+1)-by-(2n+1) estimated parameter covariance matrix
%
%   Variance         T-by-n conditional variances obtained from univariate
%                    GARCH-MIDAS models
%
%   LongRunVar       T-by-n long-run component of the conditional variances 
%                    obtained from univariate GARCH-MIDAS models
%
%   CorrMatrix       n-by-n-by-T conditional correlation matrices
%
%   LongRunCorrMatrix n-by-n-by-T long-run component of the correlation matrices
%
%   logL             T-by-1 log likelihood. Initial observations may be
%                    assigned a flag of zero.
%
%
% Reference:
%
% Colacito,R., Engle,R.F. and Ghysels, E. (2011), A Component Model for
% Dynamic Correlations, Journal of Econometrics, 164, 45-59.
%
% Written by Hang Qian
% Contact: matlabist@gmail.com


% Data dimension
dim = size(Data,2);

% Parse inputs and set defaults.
periodDefault = 22;
nlagVarDefault = 10;
nlagCorrDefault = 10;
callerName = 'DccMidas';
parseObj = inputParser;
addParameter(parseObj,'NumLagsVar',nlagVarDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'NumLagsCorr',nlagCorrDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'Period',periodDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'EstSample',[],@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'RollWindow',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'LogTau',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Beta2Para',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Options',[],@(x)validateattributes(x,{'struct','optim.options.Fmincon'},{},callerName));
addParameter(parseObj,'Mu0',[],@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'Alpha0',0.05*ones(dim,1),@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'Beta0',0.9*ones(dim,1),@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'Theta0',0.1*ones(dim,1),@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'W0',5*ones(dim,1),@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'M0',0.01*ones(dim,1),@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'CorrA0',0.05,@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'CorrB0',0.9,@(x)validateattributes(x,{'numeric'},{'vector'},callerName));
addParameter(parseObj,'CorrW0',5,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'MorePara',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Gradient',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'AdjustLag',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'ThetaM',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'ZeroLogL',[],@(x)validateattributes(x,{'numeric'},{'2d','integer'},callerName));
addParameter(parseObj,'ZeroLogLStep1',[],@(x)validateattributes(x,{'numeric'},{'2d','integer'},callerName));
parse(parseObj,varargin{:});
nlagVar = parseObj.Results.NumLagsVar;
nlagCorr = parseObj.Results.NumLagsCorr;
period = parseObj.Results.Period;
estSample = parseObj.Results.EstSample;
rollWindow = parseObj.Results.RollWindow;
logTau = parseObj.Results.LogTau;
beta2Para = parseObj.Results.Beta2Para;
options = parseObj.Results.Options;
mu0 = parseObj.Results.Mu0;
alpha0 = parseObj.Results.Alpha0;
beta0 = parseObj.Results.Beta0;
theta0 = parseObj.Results.Theta0;
w0 = parseObj.Results.W0;
m0 = parseObj.Results.M0;
corrA0 = parseObj.Results.CorrA0;
corrB0 = parseObj.Results.CorrB0;
corrW0 = parseObj.Results.CorrW0;
morePara = parseObj.Results.MorePara;
autodiffFlag = parseObj.Results.Gradient;
adjustLag = parseObj.Results.AdjustLag;
thetaM = parseObj.Results.ThetaM;
zeroLogLStep1 = parseObj.Results.ZeroLogLStep1;
zeroLogL = parseObj.Results.ZeroLogL;

% Replace missing values by the sample average
for m = 1:dim
    mask = isnan(Data(:,m));
    if any(mask)        
        Data(mask,m) = nanmean(Data(:,m));
    end
end

% Divide data as the estimation and forecast samples
if isempty(estSample)
    estSample = size(Data,1);
end
DataFull = Data;
Data = Data(1:estSample,:);
nobs = size(Data,1);

% Starting values
if isempty(mu0); mu0 = mean(Data); end
if isscalar(mu0); mu0 = repmat(mu0,dim,1); end
if isscalar(alpha0); alpha0 = repmat(alpha0,dim,1); end
if isscalar(beta0); beta0 = repmat(beta0,dim,1); end
if isscalar(theta0); theta0 = repmat(theta0,dim,1); end
if isscalar(w0); w0 = repmat(w0,dim,1); end
if isscalar(m0); m0 = repmat(m0,dim,1); end

% First step estimation: univariate GARCH-MIDAS for variances
% Refer to Eq (2.3) - (2.6) in Colacito et al. (2011)
if beta2Para
    estParamsStep1 = zeros(7,dim);
    EstParamCovStep1 = zeros(7,7,dim);
else
    estParamsStep1 = zeros(6,dim);
    EstParamCovStep1 = zeros(6,6,dim);
end
Variance = zeros(nobs,dim);
for m = 1:dim
    [estParamsStep1(:,m),EstParamCovStep1(:,:,m),Variance(:,m)] = GarchMidas(Data(:,m),...
        'Period',period,'NumLags',nlagVar,'RollWindow',rollWindow,...
        'LogTau',logTau,'Beta2Para',beta2Para,'Options',options,...
        'Mu0',mu0(m),'Alpha0',alpha0(m),'Beta0',beta0(m),'Theta0',theta0(m),'W0',w0(m),'M0',m0(m),...
        'Gradient',autodiffFlag,'AdjustLag',adjustLag,'ThetaM',thetaM,'ZeroLogL',zeroLogLStep1);    
end

% Construct standardized residuals (T-by-n matrix)
Resid = bsxfun(@minus,Data,estParamsStep1(1,:)) ./ sqrt(Variance);

% Construct cross-product of residuals
CrossResid = zeros(dim,dim,nobs);
for t = 1:nobs
    CrossResid(:,:,t) = Resid(t,:)' * Resid(t,:);
end

% Construct residual correlation matrix
% ResidCorr is reshaped as a T-by-n^2 matrix for efficient computation
% Essentially, it is corr(Resid(t-period+1:t,:))
% Refer to the c(i,j,t) term in Eq (2.7) of Colacito et al. (2011)
ResidCorr = zeros(nobs,dim^2);
if rollWindow
    
    % correlation matrix varies every period due to rolling windows
    for t = period:nobs
        partialSum = nansum(CrossResid(:,:,t-period+1:t),3);
        partialSd = sqrt(diag(partialSum));
        CorrTemp = partialSum ./ (partialSd * partialSd');
        ResidCorr(t,:) = CorrTemp(:);
    end
    
    % Presample correlation matrix
    ResidCorr(1:period-1,:) = repmat(ResidCorr(period,:),period-1,1);
    
else
    
    % correlation matrix is fixed for a week/month/quarter/year
    for t = period:period:nobs
        partialSum = nansum(CrossResid(:,:,t-period+1:t),3);
        partialSd = sqrt(diag(partialSum));
        CorrTemp = partialSum ./ (partialSd * partialSd');
        ResidCorr(t-period+1:t,:) = repmat(CorrTemp(:)',period,1);
    end
    
    % The last week/month/quarter/year may have odd days
    oddDay = mod(nobs,period);
    if oddDay > 0
        partialSum = nansum(CrossResid(:,:,end-period+1:end),3);
        partialSd = sqrt(diag(partialSum));
        CorrTemp = partialSum ./ (partialSd * partialSd');
        ResidCorr(end-oddDay+1:end,:) = repmat(CorrTemp(:)',oddDay,1);
    end
    
end

% Scalar expansion for the second step estimation
if morePara
    if isscalar(corrA0)
        corrA0 = repmat(corrA0,dim,1);
    end
    if isscalar(corrB0)
        corrB0 = repmat(corrB0,dim,1);
    end    
end

% MLE specification
if beta2Para
    params0 = [corrA0(:);corrB0(:);corrW0;corrW0];
else
    params0 = [corrA0(:);corrB0(:);corrW0];
end
if morePara
    if beta2Para
        lbMat = repmat([0 0],dim,1);
        ubMat = repmat([1 1],dim,1);
        lb = [lbMat(:);1.001;1.001];
        ub = [ubMat(:);50;50]; 
    else
        lbMat = repmat([0 0],dim,1);
        ubMat = repmat([1 1],dim,1);
        lb = [lbMat(:);1.001];
        ub = [ubMat(:);50];        
    end        
else
    if beta2Para
        lb = [0 0 1.001 1.001]';
        ub = [1 1 50 50]';
    else
        lb = [0 0 1.001]';
        ub = [1 1 50]';
    end
end

% Maximum likelihood estimation
if autodiffFlag
    if isempty(options)
        options = optimoptions('fmincon','GradObj','on','DerivativeCheck','on','FinDiffType','central','Algorithm','interior-point','Display','Iter');
        % options = optimoptions('fmincon','GradObj','on','Algorithm','interior-point','Display','Iter');
    end
    myfun = @(params)fNegMLGrad(params,Resid,CrossResid,ResidCorr,period,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow);       
else
    if isempty(options)
        options = optimoptions('fmincon','Algorithm','interior-point','Display','Iter');
    end
    myfun = @(params) - nansum(fML(params,Resid,CrossResid,ResidCorr,period,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow));
end

% Use FMINCON for numerical optimization
try
    estParamsStep2 = fmincon(myfun,params0,[],[],[],[],lb,ub,[],options);
catch
    warning('FMINCON failed... Switch to FMINSEARCH.')
    estParamsStep2 = fminsearch(myfun,params0);
end    

% Compute MLE covariance matrix
[logL,Gradient] = fMLGrad(estParamsStep2,Resid,CrossResid,ResidCorr,period,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow); 
BHHH = Gradient' * Gradient;
if rcond(BHHH) < 1e-12
    warning('Covariance matrix of estimators cannot be computed precisely due to inversion difficulty.')
end
EstParamCovStep2 = inv(BHHH);

% Compute standard errors
se2 = diag(EstParamCovStep2);
if any(se2<0)     
    se2(se2<0) = NaN;
end
se = sqrt(se2);
zstat = estParamsStep2(:) ./ se;
pval = 0.5 * erfc(0.7071 * abs(zstat)) * 2;
pval(pval<1e-6) = 0;

% AIC and BIC
logLikeSum = nansum(logL);
numParam = length(params0);
aic = -2*logLikeSum + 2*numParam;
bic = -2*logLikeSum + numParam.*log(nobs);
adjustSampleSize = sum(logL~=0 & ~isnan(logL));

% Display the estimation results
if autodiffFlag
    fprintf('Method: Two-Step Maximum likelihood (Gradient)\n');
else
    fprintf('Method: Two-Step Maximum likelihood\n');
end
fprintf('Sample size: %d\n',nobs);
fprintf('Adjusted sample size: %d\n',adjustSampleSize);
fprintf('Logarithmic  likelihood: %12.6g\n',logLikeSum);
fprintf('Akaike   info criterion: %12.6g\n',aic);
fprintf('Bayesian info criterion: %12.6g\n',bic);
columnNames = {'Coeff','StdErr','tStat','Prob'};

if morePara
    if beta2Para
        rowNames = cell(2*dim+2,1);
        for m = 1:dim
            rowNames{m} = sprintf('a(%d)',m);
            rowNames{dim+m} = sprintf('b(%d)',m);            
        end
        rowNames{2*dim+1} = 'w1';
        rowNames{2*dim+2} = 'w2';
    else        
        rowNames = cell(2*dim+1,1);
        for m = 1:dim
            rowNames{m} = sprintf('a(%d)',m);
            rowNames{dim+m} = sprintf('b(%d)',m);            
        end
        rowNames{2*dim+1} = 'w';
    end
else
    if beta2Para
        rowNames = {'a';'b';'w2';'w1'};
    else
        rowNames = {'a';'b';'w'};
    end
end

try
    Table = table(estParamsStep2,se,zstat,pval,'RowNames',rowNames,'VariableNames',columnNames);    
catch
    Table = [estParamsStep2,se,zstat,pval];
end
disp(Table)

% Construct full-sample conditional variance 
nobsBig = size(DataFull,1);
Variance = zeros(nobsBig,dim);
LongRunVar = zeros(nobsBig,dim);
for m = 1:dim
    [~,~,Variance(:,m),LongRunVar(:,m)] = GarchMidas(DataFull(:,m),'Params',estParamsStep1(:,m),...
        'Period',period,'NumLags',nlagVar,'RollWindow',rollWindow,...
        'LogTau',logTau,'Beta2Para',beta2Para,...        
        'Gradient',autodiffFlag,'AdjustLag',adjustLag,'ThetaM',thetaM,'ZeroLogL',zeroLogLStep1);    
end

% Construct full-sample correlation matrix
Resid = zeros(nobsBig,dim);
for m = 1:dim
    Resid(:,m) = (DataFull(:,m) - estParamsStep1(1,m)) ./ sqrt(Variance(:,m));
end

CrossResid = zeros(dim,dim,nobsBig);
for t = 1:nobsBig
    CrossResid(:,:,t) = Resid(t,:)' * Resid(t,:);
end

ResidCorr = zeros(nobsBig,dim^2);
for t = period+1:nobsBig
    partialSum = sum(CrossResid(:,:,t-period:t-1),3);
    partialSd = sqrt(diag(partialSum));
    CorrTemp = partialSum ./ (partialSd * partialSd'); 
    ResidCorr(t,:) = CorrTemp(:);
end
sampleCorr = corr(Resid);
ResidCorr(1:period,:) = repmat(sampleCorr(:)',period,1);

[~,CorrMatrix,LongRunCorrMatrix] = fML(estParamsStep2,Resid,CrossResid,ResidCorr,period,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow);

end


%-------------------------------------------------------------------------
% Subfunction: MIDAS beta polynomial weights
function weights = midasBetaWeights(nlag,param1,param2)
seq = nlag:-1:1; 
if isempty(param2)    
    weights = (1-seq./nlag+10*eps).^(param1-1);    
else
    weights = (1-seq./nlag+10*eps).^(param1-1) .* (seq./nlag).^(param2-1);    
end
weights = weights ./ nansum(weights);
end


%-------------------------------------------------------------------------
% Subfunction: MIDAS beta polynomial weights
function weights = midasBetaWeightsCapital(nlag,param1,param2)
seq = nlag:-1:1;
if isempty(param2)
    weights = POWER(1-seq./nlag+10*eps, param1-1);
else
    weights = TIMES( POWER(1-seq./nlag+10*eps, param1-1), POWER(seq./nlag, param2-1)) ;
end
weights = RDIVIDE(weights,nansum(weights));
end


%-------------------------------------------------------------------------
% Subfunction: Compute log likelihood of GARCH-MIDAS
%              Also return conditional variance as a by-product
function [logL,CorrMatrix,LongRunCorrMatrix] = fML(params,Resid,CrossResid,ResidCorr,period,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow)

% Input Arguments:
%
% params:     Parameter vector [a;b;w] 
% Resid:      T-by-n standardized residuals obtained from univariate GARCH-MIDAS 
% CrossResid: n-by-b-by-T squared standardized residuals
% ResidCorr:  T-by-n^2 rolling-window sample correlation matrix from period t-freq to t-1
% freq:       How many days in a week/month/quarter/year
% nlagCorr:   How many weeks/months/quarters/years for the MIDAS average
% adjustLag:  Logical value on MIDAS lag adjustment for initial observations
% zeroLogL:   a vector of indices that force zeros of the likelihood values 
% beta2Para: logical value on two parameter MIDAS Beta polynomial
% rollWindow: Logical value on rolling window long-run component
%
% Output Arguments:
%
% logL       T-by-1 log likelihood
% CorrMatrix n-by-n-by-T conditional correlation matrix
% LongRunCorrMatrix n-by-n-by-T long-run component of correlation matrix
%
% Notes:
%
% Due to missing presample values, correlation matrix is initialized by the
% identity matrix

[nobs,dim] = size(Resid);
morePara = numel(params) > 4;

% Allocate parameters
if ~morePara    
    corrA0 = params(1);
    corrB0 = params(2);
    corrW10 = params(3);
    if beta2Para
        corrW20 = params(4);
    else
        corrW20 = [];
    end
else
    aSqrt = sqrt(params(1:dim));
    bSqrt = sqrt(params(dim+1:2*dim));    
    corrA0 = aSqrt * aSqrt';
    corrB0 = bSqrt * bSqrt';
    corrW10 = params(2*dim+1);
    if beta2Para
        corrW20 = params(2*dim+2);
    else
        corrW20 = [];
    end
end

% Pre-allocation and initialization
QuasiCorrMatrix = eye(dim);
CorrMatrix = zeros(dim,dim,nobs);
CorrMatrix(:) = repmat(eye(dim),1,nobs);
LongRunCorrMatrix = zeros(dim,dim,nobs);
LongRunCorrMatrix(:) = repmat(mean(ResidCorr,1),1,nobs);
logL = zeros(nobs,1);

% Positive constraint
intercept = 1 - corrA0 - corrB0;
if any(intercept < 0)
    logL = -Inf(nobs,1);    
    return
end

% Compute MIDAS weights
nlag = period*nlagCorr;
if rollWindow    
    weights = midasBetaWeights(nlag,corrW10,corrW20);
else
    weights = midasBetaWeights(nlagCorr,corrW10,corrW20);
end

% Conditional correlation matrix evolution
if rollWindow && adjustLag
    loopStart = 2;
elseif rollWindow && ~adjustLag
    loopStart = nlag + 1;
elseif ~rollWindow && adjustLag
    loopStart = period + 1;
else
    loopStart = nlag + 1;
end
  
for t = loopStart:nobs
    
    %---------------------
    % Long-run correlation
    %---------------------
    
    % Rolling window: long-run correlation will be updated every day
    if rollWindow
        if adjustLag && (t <= nlag+1)
            weights = midasBetaWeights(t-1,corrW10,corrW20);
            LongRunVec = weights * ResidCorr(1:t-1,:);
        else
            LongRunVec = weights * ResidCorr(t-nlag:t-1,:);
        end
        LongRun = reshape(LongRunVec,dim,dim);
        LongRunCorrMatrix(:,:,t) = LongRun;
    end
    
    % Fixed length: long-run correlation will be updated the first day of a week/month/quarter/year
    if ~rollWindow && mod(t,period) == 1
        if adjustLag && (t <= nlag+1)
            nlagReduce = (t-1)/period;
            weights = midasBetaWeights(nlagReduce,corrW10,corrW20);
            LongRunVec = weights * ResidCorr(1:period:t-1,:);
        else
            LongRunVec = weights * ResidCorr(t-nlag:period:t-1,:);
        end
        LongRun = reshape(LongRunVec,dim,dim);
        LongRunCorrMatrix(:,:,t:t+period-1) = reshape(repmat(LongRunVec,1,period),dim,dim,period);        
    end
    
    %-----------------------------
    % Correlation matrix recursion
    %-----------------------------
        
    % Compute conditional correlation matrix (unnormalised)
    % Refer to Eq (2.7) in Colacito et al. (2011)    
    QuasiCorrMatrix = LongRun .* intercept + corrA0 .* CrossResid(:,:,t-1) + corrB0 .* QuasiCorrMatrix;
    QuasiCorrMatrix = 0.5 .* (QuasiCorrMatrix + QuasiCorrMatrix');
    
    % Normalize conditional correlation matrix
    % Refer to Eq (2.9) and (2.10) in Colacito et al. (2011)
    stdQ = sqrt(diag(QuasiCorrMatrix));
    CorrMatrix(:,:,t) = QuasiCorrMatrix ./ (stdQ * stdQ'); 
        
    % Check positive definiteness
    [CholP,flag] = chol(CorrMatrix(:,:,t),'lower');
    if flag > 0
        logL = -Inf(nobs,1);
        return
    end
    
    %---------------
    % Log likelihood
    %---------------
    
    % Compute the log likelihood for the standardized residuals
    constant = -0.5 * dim * log(2*pi);
    detTerm = -sum(log(diag(CholP)));
    yTran = CholP \ Resid(t,:)';
    expTerm = -0.5 * (yTran' * yTran);
    logL(t) = constant + detTerm + expTerm;
    
end

% User specified zero logL
logL(zeroLogL) = 0;

% Cut the extra length if there are odd days in a week/month/quarter/year
LongRunCorrMatrix = LongRunCorrMatrix(:,:,1:nobs);

end

%-------------------------------------------------------------------------
function logL = fMLCapital(params,Resid,CrossResid,ResidCorr,period,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow) 
 
[nobs,dim] = size(Resid); 
morePara = numel(params) > 4; 
 
% Allocate parameters 
if ~morePara 
    corrA0 = params(1); 
    corrB0 = params(2); 
    corrW10 = params(3); 
    if beta2Para
        corrW20 = params(4);
    else
        corrW20 = [];
    end
else 
    aSqrt = SQRT(params(1:dim)); 
    bSqrt = SQRT(params(dim+1:2*dim));     
    corrA0 = MTIMES(aSqrt,aSqrt.'); 
    corrB0 = MTIMES(bSqrt,bSqrt.'); 
    corrW10 = params(2*dim+1);
    if beta2Para
        corrW20 = params(2*dim+2);
    else
        corrW20 = [];
    end    
end 

% Pre-allocation and initialization
QuasiCorrMatrix = eye(dim);
CorrMatrix = zeros(dim,dim,nobs);
CorrMatrix(:) = repmat(eye(dim),1,nobs);
logL = zeros(nobs,1);

% Positive constraint
intercept = 1 - corrA0 - corrB0;
if any(intercept < 0)
    logL = -Inf(nobs,1);    
    return
end

% Compute MIDAS weights 
nlag = period*nlagCorr;
if rollWindow    
    weights = midasBetaWeightsCapital(nlag,corrW10,corrW20); 
else
    weights = midasBetaWeightsCapital(nlagCorr,corrW10,corrW20);     
end

% Conditional correlation matrix evolution 
if rollWindow && adjustLag
    loopStart = 2;
elseif rollWindow && ~adjustLag
    loopStart = nlag + 1;
elseif ~rollWindow && adjustLag
    loopStart = period + 1;
else
    loopStart = nlag + 1;
end

for t = loopStart:nobs 
    
    %---------------------
    % Long-run correlation
    %---------------------
    
    % Rolling window: long-run correlation will be updated every day
    if rollWindow
        if adjustLag && (t <= nlag+1)            
            weights = midasBetaWeightsCapital(t-1,corrW10,corrW20);
            LongRunVec = MTIMES(weights,ResidCorr(1:t-1,:)); 
        else
            LongRunVec = MTIMES(weights,ResidCorr(t-nlag:t-1,:));            
        end
        LongRun = reshape(LongRunVec,dim,dim);
    end
    
    % Fixed length: long-run correlation will be updated the first day of a week/month/quarter/year
    if ~rollWindow && mod(t,period) == 1
        if adjustLag && (t <= nlag+1)
            nlagReduce = (t-1)/period;            
            weights = midasBetaWeightsCapital(nlagReduce,corrW10,corrW20);
            LongRunVec = MTIMES(weights,ResidCorr(1:period:t-1,:));             
        else
            LongRunVec = MTIMES(weights,ResidCorr(t-nlag:period:t-1,:));            
        end
        LongRun = reshape(LongRunVec,dim,dim);
    end
    
    %-----------------------------
    % Correlation matrix recursion
    %-----------------------------
    
    % Compute conditional correlation matrix (unnormalised) 
    % Refer to Eq (2.7) in Colacito et al. (2011) 
    QuasiCorrMatrix = TIMES(LongRun,intercept) + TIMES(corrA0,CrossResid(:,:,t-1)) + TIMES(corrB0,QuasiCorrMatrix); 
    QuasiCorrMatrix = 0.5 .* (QuasiCorrMatrix + QuasiCorrMatrix.'); 
 
    % Normalize conditional correlation matrix 
    % Refer to Eq (2.9) and (2.10) in Colacito et al. (2011) 
    stdQ = SQRT(diag(QuasiCorrMatrix)); 
    CorrMatrix(:,:,t) = RDIVIDE(QuasiCorrMatrix,MTIMES(stdQ,stdQ.')); 
 
    % Check positive definiteness 
    [CholP,flag] = CHOL(CorrMatrix(:,:,t),'lower'); 
    if flag > 0 
        logL = -Inf(nobs,1); 
        return 
    end 
 
    % Compute the log likelihood for the standardized residuals 
    constant = -0.5 * dim * log(2*pi); 
    detTerm = -sum(LOG(diag(CholP))); 
    yTran = MLDIVIDE(CholP,Resid(t,:).'); 
    expTerm = -0.5*MTIMES(yTran.',yTran); 
    logL(t) = constant + detTerm + expTerm;  
end

% User specified zero logL
logL(zeroLogL) = 0;

end


%-------------------------------------------------------------------------
% Subfunction: Compute log likelihood and gradient
function [logL,Gradient] = fMLGrad(params,Resid,CrossResid,ResidCorr,freq,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow)

logL = fML(params,Resid,CrossResid,ResidCorr,freq,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow);
if nargout == 1    
    return
end

fun = @(params) fMLCapital(params,Resid,CrossResid,ResidCorr,freq,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow);
x = params;

nfunVal = numel(logL);
nparams = numel(x);
fval = zeros(nfunVal,nparams);
Gradient = zeros(nfunVal,nparams);

for m = 1:nparams
    
    % One element of x carries an imaginary number
    xComplex = x;
    xComplex(m) = xComplex(m) + 1i;
    
    % The real component is the function value 
    % The imaginary component is the derivative
    complexValue = fun(xComplex);
    fval(:,m) = real(complexValue);
    Gradient(:,m) = imag(complexValue);
end

if any(norm(bsxfun(@minus,fval, logL))>1e-6)
    warning('Automatic gradients might be falsely computed. Capitalize arithmetic operator and utility function names.')
end

end


%-------------------------------------------------------------------------
% Subfunction: Compute negative log likelihood (scalar) and gradient
function [logL,Gradient] = fNegMLGrad(params,Resid,CrossResid,ResidCorr,freq,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow)
[logL,Gradient] = fMLGrad(params,Resid,CrossResid,ResidCorr,freq,nlagCorr,adjustLag,zeroLogL,beta2Para,rollWindow);
logL = - nansum(logL,1);
Gradient = - nansum(Gradient,1);
end


%-------------------------------------------------------------------------
function [A,psd] = CHOL(A, flag)
if nargin < 2
    flag = 'upper';
end
tol = 1e-8;
dim = size(A,1);

% First, zero the upper triangular portion of the matrix.
% If A is not symmetric, there is no error message. It is treated as if
% the upper triangular portion was a mirror to the lower triangular part.
A = tril(A);

% Recursively find the cholesky factor block by block
for r = 1:dim
    
    % When A(r,r) is in a neighborhood of zero, treat as if it was zero
    % When A(r,r) is too negative, it is not positive semi-definite.    
    if A(r,r) < tol      
        if A(r,r) < -tol            
            if nargout <= 1
                error('Not a positive semi-definite matrix.')
            else
                psd = 100;
                return
            end                
        end
        A(r:end,r) = 0;       
        continue
    end
    
    % Block (1,1), the scalar of the square root of the original value
    A(r,r) = SQRT(A(r,r));
    
    if r == dim
        break
    end
    
    % Block (1,2), the row vector of zeros   
    
    % Block (2,1), the column vector of deflated original values
    aInv = RDIVIDE(1.0, A(r,r));
    A(r+1:dim,r) = TIMES(A(r+1:dim,r), aInv);
    
    % Block (2,2), the matrix of rank 1 update
    Rank1 = tril(MTIMES(A(r+1:dim,r), A(r+1:dim,r).'));
    A(r+1:dim,r+1:dim) = A(r+1:dim,r+1:dim) - Rank1;
end

if strcmpi(flag,'upper')
    A = A.';
end

psd = 0;
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
function B = INV(A)
Areal = real(A);
Aimag = imag(A);
Breal = inv(Areal);
Bimag = -Breal * Aimag * Breal;
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
function C = MLDIVIDE(A,B)
Creal = real(A) \ real(B);
if size(A,1) == size(A,2)
    Cimag = imag(MTIMES(INV(A), B));
else
    Cimag = imag(MTIMES(INV(MTIMES(A.',A)), (MTIMES(A.',B))));  
end
C = complex(Creal, Cimag);
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
elseif all(A(:)>=0)
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
function B = SQRT(A)
Areal = real(A);
Aimag = imag(A);
Breal = sqrt(Areal);
Bimag = 0.5./Breal .* Aimag;
B = complex(Breal, Bimag);
end


%-------------------------------------------------------------------------
function C = TIMES(A,B)
C = A .* B + imag(A) .* imag(B);
end

