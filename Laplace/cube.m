function [X Y m h hk w k Nn]=cube(N,xa,xb,ya,yb)
Nn{1}=N;
k=(N+1)/2;
sx=abs(xb-xa);
sy=abs(yb-ya);
s=sqrt(sx^2+sy^2);
hs=(s/2)/(k-1);
hk=hs/(sqrt(2));
[zs zc ang m{1} h{1}]=irpolygon(4,N,0,[xa xb xb xa],[ya ya yb yb]);
X{1}=[linspace(xa,xb,N)' ; xb*ones(N,1); linspace(xb,xa,N)'; xa*ones(N,1)]; 
Y{1}=[ya*ones(N,1); linspace(ya,yb,N)' ; yb*ones(N,1); linspace(yb,ya,N)'];

for i=2:k
    N=N-2;
    Nn{i}=N;
    xa=xa+hk;xb=xb-hk;ya=ya+hk;yb=yb-hk;
    [zs zc ang m{i} h{i}]=irpolygon(4,N,0,[xa xb xb xa],[ya ya yb yb]);
    X{i}=[linspace(xa,xb,N)' ; xb*ones(N,1); linspace(xb,xa,N)'; xa*ones(N,1)]; 
    Y{i}=[ya*ones(N,1); linspace(ya,yb,N)' ; yb*ones(N,1); linspace(yb,ya,N)'];
end
w=ang-(pi/2);