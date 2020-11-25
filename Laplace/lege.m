function p=lege(n,x)
p=0;
for k=0:n
    p=p+ (bicoeff(n,k).^2).*((x-1).^(n-k)).*((x+1).^k);
end
p=(1/(2^n))*p;
    