function varargout=almon_adl_new(y,x,poly,output_type,polyconstr)
%Almon lag function
%Arthur Sinko, UNC, 2007
% last modified Arthur Sinko 12/22/2009
assert(isstruct(x),'x should be a structure with x.xmidas and x.xols fields');
xmidas=x.xmidas;
xols=x.xols;

if isempty(xmidas)
    [aopt, error, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT]=ols(y,xols);
    weights=[];
    stderrors_almon = [];
    tstats_almon = [];
    robustSE_almon = [];
    robustT_almon = [];
    [aic bic]=aicbic(error,length(aopt));
    if isequal(output_type,'struct')
        weightsNN=[];
        aoptxols = aopt(end-size(x.xols,2)+1:end);
        aoptmidas = aopt(1:end-size(x.xols,2));
        varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error, weights, weightsNN, estvar, stderror,Tvalues,R2,robustV,robustSE,robustT)};
    else
        varargout={aopt, weights, stderrors_almon, stderror, tstats_almon, Tvalues, R2, error, robustV, robustSE, robustSE_almon, robustT, robustT_almon};
    end
    return
end



% X=xmidas*S;
[aopt, error, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT]=ols(y,[midas_X(x,'almon',poly,polyconstr) xols]);

RSS = error'*error;

weights = almon_weights(size(xmidas,2),poly);
%weights=[];%S*aopt(2:end-size(xols,2));
stderrors_almon = [];%sqrt(diag(S*estvar(2:end-size(xols,2),2:end-size(xols,2))*S'));
tstats_almon = [];%weights./stderrors_almon;

robustSE_almon = [];%sqrt(diag(S*robustV(2:end-size(xols,2),2:end-size(xols,2))*S'));
robustT_almon = [];%weights./robustSE_almon;

[aic bic]=aicbic(error,length(aopt));
if isequal(output_type,'struct')
    weightsNN=weights;
    aoptxols = aopt(end-size(x.xols,2)+1:end);
    aoptmidas = aopt(1:end-size(x.xols,2));
    varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error,RSS, weights, weightsNN, estvar, stderror,Tvalues,R2,R3, robustV,robustSE,robustT)};
else
    varargout={aopt,RSS, weights, stderrors_almon, stderror, tstats_almon, Tvalues, R2, R3, error, robustV, robustSE, robustSE_almon, robustT, robustT_almon};
end
