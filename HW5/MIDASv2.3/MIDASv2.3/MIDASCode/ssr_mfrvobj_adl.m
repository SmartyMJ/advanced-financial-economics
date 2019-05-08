function SSR=ssr_mfrvobj_adl(a,x,y,polyconstr)
%Residual sum of squares for Exponential polynomial
SSR=mfrvobj_adl(a,x,y,polyconstr)'*mfrvobj_adl(a,x,y,polyconstr);
end
