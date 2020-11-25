function S=legtran(z,l)
if z==0
    if l==0
        S=2;
    else
        S=0;
    end
%     S=0;
%     for k=0:l
%         fun=@(t) (exp(z*t)/(2^l)).*((bicoeff(l,k)).^2).*...
%             ((t-1).^(l-k)).*((t+1).^k);
%         S=S+quadgk(fun,-1,1);
%     end
else
    S=(sqrt(2*pi*z)/z)*besseli(l+0.5,z);
end