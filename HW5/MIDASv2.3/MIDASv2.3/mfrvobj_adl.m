function [err,jac,w ] =mfrvobj_adl(atotal,x,y,polyconstr)
%Errors for Exponential MIDAS polynomial
% xmidas=x.xmidas;
% xols=x.ylagf
[dayLag]=size(x.xmidas,2);
a=atotal(1:end-size(x.xols,2));
alpha=a(1);
beta=a(2);
k1=a(3);
aols=atotal(end-size(x.xols,2)+1:end);
k2=0;
if length(a)==4
    k2=a(4);
end
params=[k1 k2];

% w=exp_weights(dayLag,k1,k2);
[xx w]=midas_X(x,'exp',params,polyconstr);
if ~isempty(x.xols)
err=(y-alpha-beta*xx-x.xols*aols(:));
else
err=(y-alpha-beta*xx);
end

jac=-[ones(size(xx)) xx a(2)*jacob(@midas_X,3,x,'exp',params,polyconstr) x.xols];

%RSS=err'*err;
%err'*err
%a
