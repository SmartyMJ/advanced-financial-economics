function varargout = midas_sf_adl_new(y,x,str,output_type,polyconstr)
% [coeffs,error, estvar, stderror,Tvalues, weights_midas ,R2,IND,robustV,robustSE,robustT] = midas_sf_adl_new(y,x,str)
% Step polinomial MiDaS estimation
% Inputs:
% y - LHS for MiDaS regression
% x - RHS for MiDaS regression
% Arthur Sinko, UNC, 2007
% last modified Arthur Sinko 12/22/2009

% options0=optimset('LargeScale','off','Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'MaxFunEvals',500);
assert(isstruct(x),'x should be a structure with x.xmidas and x.xols fields');
xmidas=x.xmidas;
xols=x.xols;

if isempty(xmidas)
    [aopt, error, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT] =ols(y,xols);
    weights=[];
    [aic bic]=aicbic(error,length(aopt));
    if isequal(output_type,'struct')
        weightsNN=[];
        aoptxols = aopt(end-size(x.xols,2)+1:end);
        aoptmidas = aopt(1:end-size(x.xols,2));
        varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error, weights, weightsNN, estvar, stderror,Tvalues,R2,robustV,robustSE,robustT)};
    else
        varargout={aopt, error, weights, estvar, stderror,Tvalues,R2,robustV,robustSE,robustT};
    end
     
    return
end




[aopt, error, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT] =ols(y,[midas_X(x,'sfun',str,polyconstr) xols]);
% const=weights(1);
weights=[];%IND*aopt(2:end-size(xols,2));

RSS = error'*error;

    [aic bic]=aicbic(error,length(aopt));
    if isequal(output_type,'struct')
        weightsNN=weights;
        aoptxols = aopt(end-size(x.xols,2)+1:end);
        aoptmidas = aopt(1:end-size(x.xols,2));
        varargout={var2struct(aic,bic,aopt,aoptxols,aoptmidas, error,RSS, weights, weightsNN, estvar, stderror,Tvalues,R2,R3, robustV,robustSE,robustT)};
    else
        varargout={aopt, error,RSS, weights, estvar, stderror,Tvalues,R2,R3, robustV,robustSE,robustT};
    end
    


