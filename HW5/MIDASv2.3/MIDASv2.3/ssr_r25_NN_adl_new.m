function [SSR1,grad,error1,be_values,jac]=ssr_r25_NN_adl_new(atotal,x,y,polyconstr)
% Part of MIDAS example Computes sum of squares of residuals for Beta
% polynomial restricted/unrestricted MIDAS specification with a non-zero last
% coefficient.
% Arthur Sinko, UNC, 2007
%Last update 12/07/2010

xols=x.xols;
a=atotal(1:end-size(xols,2));
aols=atotal(end-size(xols,2)+1:end);
minimax=@(x)(max(1e-8,min(300,x))); 
% disp(a')
if length(a)==4
    k1=1;
    k2=minimax(a(3));
    theta=a(4);
params=[k2 theta];
else
    k1=minimax(a(3));
    k2=minimax(a(4));
    theta=a(5);
params=[k1 k2 theta];    
end

% a(3)=(a(3)*(a(3)>0)+1e-8*~(a(3)>0))*(a(3)<300)+300*~(a(3)<300);
% a(4)=(a(4)*(a(4)>0)+1e-8*~(a(4)>0))*(a(4)<300)+300*~(a(4)<300);

% be_values=betapdf(u,1,a(3));%u.^(a(3)-1).*(1-u).^(a(4)-1)./beta(a(3),a(4));
% be_values=beta_weights(smpl,k1,k2)+theta;

% be_values=be_values/sum(be_values)+theta;
[xx be_values]=midas_X(x,'betaNN',params,polyconstr);

if ~isempty(xols)
error1=y-a(1)-a(2)*xx-xols*aols(:);
else
error1=y-a(1)-a(2)*xx;
end

% error1=((y-a(1)-a(2)*y_h));

rat1=(error1'*error1);%/var(y);%./((de_long(y))'*(de_long(y)));

SSR1=rat1;

% jac_b=jac_bnlsNN_adl_new(a,xmidas,xols,be_values);
jac_e=[ones(size(xx)) xx a(2)*jacob(@midas_X,3,x,'betaNN',params,polyconstr) xols];
jac=[];
for i=1:length(atotal)
    jac=[jac -2*jac_e(:,i).*error1];
end
grad=sum(jac,1)';