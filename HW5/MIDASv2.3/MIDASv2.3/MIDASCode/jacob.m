%Helper function for MIDAS toolbox. Computes Jacobian of function handler
%fname. num is a place of parameters in fname.
%Arthur Sinko July 2010

function out=jacob(fname,num,varargin)
params0=varargin;
out=[];
eps=1e-6;
for i=1:length(params0{num})
    params1=params0;
    params2=params0;
    params1{num}(i)=params0{num}(i)+eps/2;
    params2{num}(i)=params0{num}(i)-eps/2;
    temp=(feval(fname,params1{:})-feval(fname,params2{:}))/eps;
    out=[out temp];
end
end