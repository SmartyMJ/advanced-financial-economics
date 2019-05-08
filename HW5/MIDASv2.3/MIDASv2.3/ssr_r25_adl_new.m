function [SSR1,grad,error1,be_values,jac]=ssr_r25_adl_new(atotal,x,y,polyconstr)
% This function computes sum of squares of residuals for unrestricted/restricted beta
% polynomial.
% Arthur Sinko, UNC, 2006
% Modified Dec, 2009

xols=x.xols;
a=atotal(1:end-size(xols,2));
aols=atotal(end-size(xols,2)+1:end);
minimax=@(x)(max(1e-8,min(300,x)));
% disp(a')
if length(a)==3
    k1=1;
    k2=minimax(a(3));
    params=[k2];
else
    k1=minimax(a(3));
    k2=minimax(a(4));
    params=[k1 k2];
end


% smpl=size(xmidas,2);
% be_values=beta_weights(smpl,k1,k2);
% y_h=(xmidas*be_values);
[xx be_values]=midas_X(x,'beta',params,polyconstr);


if ~isempty(xols)
error1=((y-a(1)-a(2)*xx-xols*aols(:)));
else
error1=((y-a(1)-a(2)*xx));
end

SSR1=(error1'*error1);

jac_e=[ones(size(xx)) xx a(2)*jacob(@midas_X,3,x,'beta',params,polyconstr) xols];
jac=[];
for i=1:length(atotal)
    jac=[jac -2*jac_e(:,i).*error1];
end
grad=sum(jac,1)';
