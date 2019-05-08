function [estimates, errors, estvar, stderror,Tvalues,R2,R3,robustV,robustSE,robustT]=ols(y,x,noconst)

% usage: [estimates, error, estvar, stderror,Tvalues,R2,R2adj]=ols(y,x,noconst)
%the ols(y,x, noconst) function takes a vector and
%a matrix (vector)as arguments.
%The rows of the arguments must be identical
%The y argument is the vector of of dependent variable
%The x argument is the matrix of independent variables
%The x matrix must NOT include a column of ones.
%The regression is run with a constant by default.
%IF the noconst parameter is 1, the regression is run WITHOUT a constant
%IF additional results such as errrors, estimated variance, or R2 are
%requested, they must be specified in the output matrix.

%copyright Rossen Valkanov, 02/22/96
%Last modified Arthur Sinko 9/12/2009
[r,t]=size(x); [v, w]=size(y);
if nargin<3 , noconst=0; end;
if noconst==0, xr=[ones(r,1) x]; 	  %run with a constant
else xr=x;end				              % NO constant
[r,t]=size(xr);
indic=0;
if r<=t
    indic=NaN;
end
if r~=v
    error('Incompatible matices');
end
invXpX=inv(xr'*xr);
P=invXpX*xr'+indic;	%projection matrix
estimates=P*y;       %ols estimates


if nargout>1 		%execute only if more than estimates are required
    degfr=r-t;
    yhat=xr*estimates;
    errors=y-yhat;
    RSS=errors'*errors;
    s2=RSS./degfr;
    estvar=s2*invXpX;
    stderror=sqrt(diag(estvar));
    TSS=(y-mean(y))'*(y-mean(y));			  %Total SS
    R2=1-RSS/TSS;                         %R2
    R3=1-((r-1)/(r-t))*(1-R2);            %adj R2--note: it might be negative or greater than zero, use only when many variables
    Tvalues=estimates./stderror; %tvalues for null coeff=0.
% end;
% 
% if nargout>7
    robustV=r/degfr*P*HAC_kernel(errors,'NW')*P';
    robustSE=sqrt(diag(robustV));
    robustT=estimates./robustSE; 
end