%AIC & BIC criteria
function varargout=aicbic(err,k)
n=length(err);
nlogmse=n*log(mean(err.^2));
out.aic=nlogmse+2*k;
out.bic=nlogmse+k*log(n);
if nargout==1
    varargout={out};
else
    varargout={out.aic,out.bic};
end
