function var_est=HAC_kernel(errors,ktype)
% Computation of error component of HAC estimator 
% Should be either 'white' or 'NW'
var_est=diag(errors.^2);
if isequal(lower(ktype),'white')
    return
end

if ~isequal(lower(ktype),'nw')
    error('MIDAS_HAC_kernel:unknown_kernel','Unknown type of kernel')
end

T=length(errors);
q=floor(4*(T/100)^(2/9));
for i=1:q
    G=spdiags([zeros(i,1);errors(i+1:end).*errors(1:end-i)], i, T, T);
    var_est=var_est+(G+G')*(q+1-i)/(q+1);
end
    