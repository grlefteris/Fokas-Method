%=====Binomial Coefficient (Multiplicative Formula)=======
function b=bicoeff(n,k)
prod=1;
for i=1:k
prod=prod*((n+1-i)/i);
end
b=prod;