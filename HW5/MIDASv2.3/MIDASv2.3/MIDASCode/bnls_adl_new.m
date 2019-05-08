function varargout = bnls_adl_new(y,x,numofparams,output_type,polyconstr,a0,options)
% [beta_aopt, error, weights, estvar, stderror,Tvalues,R2,robustSE,robustT] = bnls_adl_new(y,x,a0,options,LB,UB)
% Beta polinomial MiDaS estimation
% Inputs:
% a0 - initial guess of parameters. If dropped, assumed to be [alpha beta 1 5], where alpha
% and beta - intercept and slope of AR(1) process.
% options - options for fminunc() minimization routine
% Arthur Sinko, UNC, 2006
% last modified Arthur Sinko 12/17/2009

assert(isstruct(x),'x should be a structure with x.xmidas and x.xols fields');
xmidas=x.xmidas;
% xmidasd=x.xmidasd;
xols=x.xols;
if ~exist('output_type','var')
    output_type='array';
end

if isempty(xmidas)
    [beta_aopt, error, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT]=ols(y,xols);
    weights=[];
    [aic bic]=aicbic(error,length(beta_aopt));
    if isequal(output_type,'struct')
        aopt=beta_aopt;
        weightsNN=[];
        aoptxols = beta_aopt(end-size(x.xols,2)+1:end);
        aoptmidas = beta_aopt(1:end-size(x.xols,2));
        varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error, weights, weightsNN, estvar, stderror,Tvalues,R2,robustSE,robustT)};
    else
        varargout={beta_aopt, error, weights, estvar, stderror,Tvalues,R2,robustSE,robustT};
    end
    return
end

[numobs nlag]=size(xmidas);
beta0=[1; 5];
if numofparams==2
    a00 = ols(y,[midas_X(x,'beta',beta0,polyconstr) xols]);
    a00=[a00(1:2)' beta0' a00(3:end)']';
elseif numofparams==1
    a00 = ols(y,[midas_X(x,'beta',beta0,polyconstr) xols]);
    a00=[a00(1:2)' beta0(2)' a00(3:end)']';
else 
    error('numofparams should be either 1 or 2');
end

minimax=@(x)(max(1e-8,min(300,x)));
options0=optimset('LargeScale','off','Display','off','TolFun',1e-9,'TolX',1e-9,'MaxIter',5000,'MaxFunEvals',5000,'GradObj','off');
% LB0=[-5 -5 .5 .5];
% UB0=[300 300 300 300];
if nargin<2
    error('insufficient number of variables')
end
if nargin<6
    a0=a00;
    options=options0;
end
if nargin>5
    if sum((a0(3:4)<=eps)+(a0(3:4)>300))
        warning('wrong initial values')
        a0=a00;
    end
    if nargin==8
        options=options0;
    end
end
% LAMBDA.lower=0;
% LAMBDA.upper=0;
%[beta_aopt,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]=fmincon('ssr_r25',a0,[],[],[],[],LB,UB,[],options,x,y);
[a,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc('ssr_r25_adl_new',a0,options,x,y,polyconstr);
% [a,FVAL]=fminsearch('ssr_r25',a0,options,x,y);
HESSIAN=hessian('ssr_r25_adl_new',a,x,y,polyconstr);

if rcond(HESSIAN)<eps
    invHESSIAN=NaN*HESSIAN;
else
    invHESSIAN=inv(HESSIAN);
end

a(3)=minimax(a(3));
if numofparams==2
    a(4)=minimax(a(4));
end

beta_aopt=a;
% if sum(abs(LAMBDA.lower))+sum(abs(LAMBDA.upper))
%     warning('Constrains are active. Standard errors are irrelevant')
% end
[temp_5,temp_0,error,weights,jacob]=ssr_r25_adl_new(beta_aopt,x,y,polyconstr);
R2=1-cov(error)/cov((y));
n = size(y,1);
p = length(beta_aopt)-1;
R3 = R2-(1-R2)*p/(n-p-1);

RSS=temp_5;

estvar=invHESSIAN*2*var(error);%/(numobs-4);


stderror=sqrt(diag(estvar))';
Tvalues=beta_aopt(:)./stderror(:);



M= floor(4*(numobs/100)^(2/9)); %NW window
% disp('Outer gradient product.')
SS=jacob'*jacob;

for kk=1:M
    %     for obs=kk+1:t
    SS=SS+(1-kk/(M+1))*((jacob(1+kk:end,:)'*jacob(1:end-kk,:))+...
        (jacob(1:end-kk,:)'*jacob(1+kk:end,:)));
    %     end
end

HH=invHESSIAN;
robustVCV = HH*SS*HH;
robustSE=sqrt(diag(HH*SS*HH));
% disp('Robust t-value.')
robustT=a(:)./robustSE(:);
[aic bic]=aicbic(error,length(beta_aopt));
if isequal(output_type,'struct')
    aopt=beta_aopt;
    weightsNN=weights*aopt(2);
    aoptxols = beta_aopt(end-size(x.xols,2)+1:end);
    aoptmidas = beta_aopt(1:end-size(x.xols,2));
    varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error, RSS, weights, weightsNN, estvar, stderror,Tvalues,R2,R3, robustVCV ,robustSE,robustT)};
else
    varargout={beta_aopt, error, RSS, weights, estvar, stderror,Tvalues,R2,R3,robustVCV , robustSE,robustT};
end
