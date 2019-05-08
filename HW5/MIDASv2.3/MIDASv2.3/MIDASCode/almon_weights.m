function poli=almon_weights(dayLag,coeffs)
pow=length(coeffs)-1;
temp=(1:dayLag)';
S(:,pow+1)=temp.^pow;
for i=pow:-1:1
    S(:,i)=temp.^(i-1);
end
poli=S*coeffs(:);
end