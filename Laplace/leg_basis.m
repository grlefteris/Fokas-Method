function phi=leg_basis(Ncol)
x=linspace(-1,1,Ncol)';
powers=repmat((0:Ncol-1),Ncol,1);
a=repmat(x-1,1,Ncol);
b=repmat(x+1,1,Ncol);
a=a.^(fliplr(powers));
b=b.^powers;

for n=0:Ncol-1
    for k=0:n
        c(n+1,k+1)=bicoeff(n,k)^2;
    end
    
    p(:,n+1)=(1/(2^n))*...
        sum((a(:,end-n:end).*b(:,1:n+1)).*repmat(c(n+1,1:n+1),Ncol,1),2);
    
end
phi=p;

