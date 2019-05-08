function varargout = enls1_adl_new(y,x,numofparams,output_type,polyconstr,a0,options,LB,UB)
% Exponential polinomial MiDaS estimation
% [exp_aopt, errors, weights ,R2,stderror,Tvalues,robustSE,robustT] =
% enls1_adl_new(y,x,numofparams,output_type,a0,options,LB,UB)
% Inputs:
% y - LHS for MiDaS regression
% x - RHS for MiDaS regression
% a0 - initial guess of parameters. If dropped, assumed to be [alpha beta -1 0], where alpha
% and beta - intercept and slope of OLS given these values.
% options - options for fmincon() minimization routine
% LB, UB - lower and upper bound for exp_aopt.
%Arthur Sinko, UNC, 2007
% last modified Arthur Sinko 12/7/2009

xmidas=x.xmidas;
xols=x.xols;

if isempty(xmidas)
    [aopt, error, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT]=ols(y,xols);
    weights=[];
%     aopt=exp_aopt;
    [aic bic]=aicbic(error,length(aopt));
    if exist('output_type','var') && isequal(output_type,'struct')
        weightsNN=[];
        aoptxols = aopt(end-size(x.xols,2)+1:end);
        aoptmidas = aopt(1:end-size(x.xols,2));
        varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error, weights, weightsNN, estvar, stderror,Tvalues,R2,robustSE,robustT)};
    else
        varargout={aopt, error, weights ,R2,stderror,Tvalues,robustSE,robustT};
    end
    return
end

numobs=length(y);

    exp0=[-1 -0 ]';
    a00=ols(y,[midas_X(x,'exp',exp0,polyconstr) xols]);
    LB0=[-10000 -10000 -100 -100];
    UB0=[10000 10000 10 0.5];
if numofparams==1
    exp0=exp0(1);
    LB0=LB0(1:3);
    UB0=UB0(1:3);
end
a00=[a00(1:2); exp0; a00(3:end)] ;
options0=optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',1000,'Jacobian','on');

if nargin<2
    error('insufficient number of variables')
end
if nargin<6
    a0=a00;
    options=options0;
    LB=LB0;
    UB=UB0;
end
if nargin==6
    options=options0;
    LB=LB0;
    UB=UB0;
end
if nargin==7
    LB=LB0;
    UB=UB0;
end
% a0=a00(1:3)
%[exp_aopt,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]=fmincon('mfrvobj1',a0,[],[],[],[],LB,UB,[],options,x,y)
[aopt,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA, JACOBIAN]=lsqnonlin('mfrvobj_adl',a0,LB,UB,options,x,y,polyconstr);
% if sum(abs(LAMBDA.lower))+sum(abs(LAMBDA.upper))
%     warning('Constrains are active. Standard errors are irrelevant')
% end

a=aopt;
[errors,jac,weights]=mfrvobj_adl(aopt,x,y,polyconstr);
HH=hessian('ssr_mfrvobj_adl',a,x,y,polyconstr);
if rcond(HH)<eps
    invHESSIAN=HH*NaN;
else
    invHESSIAN=inv(HH);
end
% H=HH/2;
R2=1-cov(errors)/cov((y));
% R2_total=1-cov(error)/cov((y));
n = size(y,1);
p = length(aopt)-1;
R3 = R2-(1-R2)*p/(n-p-1);
RSS = errors'*errors;

estvar=invHESSIAN*2*var(errors);%/(numobs-4);
stderror=sqrt(diag(estvar))';
Tvalues=aopt'./stderror;
%
%
%
jacob=2*jac.*repmat(errors,1,length(aopt));


M= floor(4*(numobs/100)^(2/9)); %NW window
SS=jacob'*jacob;

for kk=1:M
    %     for obs=kk+1:t
    SS=SS+(1-kk/(M+1))*((jacob(1+kk:end,:)'*jacob(1:end-kk,:))+...
        (jacob(1:end-kk,:)'*jacob(1+kk:end,:)));
    %     end
end

robustVCV = invHESSIAN*SS*invHESSIAN;

robustSE=sqrt(diag(invHESSIAN*SS*invHESSIAN))';
robustT=a'./robustSE;

error=errors;
[aic bic]=aicbic(errors,length(aopt));
if isequal(output_type,'struct')
    weightsNN=weights*aopt(2);
    aoptxols = aopt(end-size(x.xols,2)+1:end);
    aoptmidas = aopt(1:end-size(x.xols,2));
    varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error, RSS, weights, weightsNN, estvar, stderror,Tvalues,R2,R3, robustVCV, robustSE,robustT)};
else
    varargout={aopt, error, RSS, weights ,R2,R3, stderror,Tvalues,robustVCV, robustSE,robustT};
    
end
