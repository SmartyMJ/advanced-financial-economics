%This function allows for various midas polynomial specifications and
%various interpolation types.
function [X midas_weights]=midas_X(x,poly,params,comb)
xmidas=x.xmidas;
xmidasd=x.xmidasd;
ind0=xmidasd>0;

% X(size(xmidasd,1),1)=0;
switch comb
    case 'rt' %real time -- only non-zero elements count, weight proportional to the distance
        for i=size(xmidasd,1):-1:1
            ind=xmidasd(i,ind0(i,:));
            wts=weights(ind(:),params,poly);
            X(i,:)=xmidas(i,:)*[wts; zeros(size(xmidas,2)-length(ind),size(wts,2))]; 
            if isequal(poly,'sfun')||isequal(poly,'almon')
                midas_weights(:,:,i)=[wts; zeros(size(xmidas,2)-length(ind),size(wts,2))];
            else
                midas_weights(i,:)=[wts; zeros(size(xmidas,2)-length(ind),size(wts,2))];
            end
        end
        if ~exist('X','var')
            X=zeros(0,1);
        end
       
    case {'es'} %zeros at the end -- all zeros are at the end of the sample
        wts=weights(size(xmidasd,2),params,poly);
        X=xmidas*wts;
        midas_weights=repmat(wts,size(xmidasd,1),1);
        
    otherwise 
        error('Polynomial specification option is wrong.')
end



end
 
function w=weights(u,params,poly)
switch poly
    case 'beta'
        if length(params)==2
            w=b_weights(u,params(1),params(2));
        elseif length(params)==1
            w=b_weights(u,1,params(1));
        else
            error('Wrong parameterization of beta polynomial');
        end
    case 'betaNN'
        if length(params)==3
            w1=b_weights(u,params(1),params(2))+params(3);
            w = w1/sum(w1);
        elseif length(params)==2
            w1=b_weights(u,1,params(1))+params(2);
            w = w1/sum(w1);
        else
            error('Wrong parameterization of betaNN polynomial');
        end
    case 'exp'
        if length(params)==2
        w=exp_weights(u,params(1),params(2));
        elseif length(params)==1
            w=exp_weights(u,params(1),0);
        else 
            error('Wrong parameterization of exp polynomial');
        end
        
    case 'LOG'
        if length(params)==2
        w=log_weights(u,params(1),params(2));
        elseif length(params)==1
            w=log_weights(u,params(1),0);
        else 
            error('Wrong parameterization of exp polynomial');
        end
        w=flipud(w);
        
    case 'sfun'
        
        if isscalar(u)
            u=1:u;
        else
            u=abs(u-u(1)-1);
        end
        
        if u(end)>params(end)
            params=[params u(end)];
        end
        w=zeros(u(end),length(params));
        w(u<=params(1),1)=1;
        for j=2:length(params)
            w((u<=params(j))&(u>params(j-1)),j)=1;
        end
        
    case 'almon'
        w=[];
        if isscalar(u)
            u=1:u;
        end
                
        for i=0:(params)
            w=[w u(:).^i];
        end

    otherwise
        error('Polynomial specification is not correct')
end
end


function [be_values]=b_weights(ind,k1,k2)
%Beta weights function.
if length(ind)>1
    u=(ind(1)-ind)/(ind(1)-ind(end));
    u(1)=u(1)+eps;
    u(end)=u(end)-eps;
end

if isscalar(ind)
    u=linspace(eps,1-eps,ind)';
end

% k1=(k1*(k1>0)+1e-8*~(k1>0))*(k1<300)+300*~(k1<300);% 0<k1<300
% k2=(k2*(k2>0)+1e-8*~(k2>0))*(k2<300)+300*~(k2<300);% 0<k1<300

be_values=u.^(k1-1).*(1-u).^(k2-1);%u.^(k1-1).*(1-u).^(k2-1)./beta(k1,k2);
be_values=be_values/sum(be_values);
end


function poli=exp_weights(dayLag,k1,k2)

if isscalar(dayLag)
    iii=transpose(1:dayLag);
else
    iii=abs(dayLag-dayLag(1)-1);
end
poli=exp(k1*iii+k2*iii.^2)/sum(exp(k1*iii+k2*iii.^2));
if isnan(sum(poli))
    
    exppoly=k1*iii+k2*iii.^2;
    for i=dayLag:-1:1
        poli(i,1)=1/sum(exp(exppoly-k1*i-k2*i.^2));
    end
end
end

function poli2=log_weights(dayLag,k1,k2)

if isscalar(dayLag)
    iii=transpose(1:dayLag);
else
    iii=abs(dayLag-dayLag(1)-1);
end
poli2=log(1+abs(k1*iii+k2*iii.^2))/sum(log(1+abs(k1*iii+k2*iii.^2)));
if isnan(sum(poli2))
    
    logpoly=k1*iii+k2*iii.^2;
    for i=dayLag:-1:1
        poli2(i,1)=1/sum(log(logpoly-k1*i-k2*i.^2));
    end
end
end