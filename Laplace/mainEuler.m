clear;
clc;

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=4;
Nn=21;
M=n*Nn;
R=2*M;
% M=2;R=5;

[X Y k hs x y hx hy w N]=domain(n,Nn,1,0);

for i=1:k-1
U{i}=uan(X{i},Y{i});
W{i}=kron(w,ones(N{i},1));
q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
hx{i}=spdiags(hx{i},0,n*N{i},n*N{i});
hy{i}=spdiags(hy{i},0,n*N{i},n*N{i});
hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
end

[sides z ang m{1} h{1}]=irpolygon(n,N{1},0,x{1},y{1});

Q{1}=restrict(n,Nn,2);
P{1}=kron(eye(n),leg_basis(Nn));

A{1}=PLHS_Dirichlet(n,N{1},M,R,m{1},h{1},1);
Phat{1}=Ptransform(n,N{1},M,R,m{1},h{1});

K{1}=Q{1}*P{1};
L{1}=Q{1}*hk{1}*P{1};
a{1}=P{1}\U{1};
b{1}=P{1}\q{1};

for i=2:k-1
    M=n*N{i};
    R=2*M;
%     M=3;R=6;
    [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});
    
    A{i}=PLHS_Dirichlet(n,N{i},M,R,m{i},h{i},1);
    Phat{i}=Ptransform(n,N{i},M,R,m{i},h{i});
    
    P{i}=kron(eye(n),leg_basis(N{i}));
    Q{i}=restrict(n,N{i},2);
    K{i}=Q{i}*P{i};
    L{i}=Q{i}*hk{i}*P{i};
    
    a{i}=P{i}\U{i};
    b{i}=P{i}\q{i};
end
%-----------coefficient matrix size-----%
SI=0;
SJ=0;
for i=1:k-1
    SI=SI+size(P{i},1)+size(A{i},1);
    SJ=SJ+(2*n*N{i});
end
%---------------------------------------%
%---------------------------coefficient matrix assembly-------------------%
MAT=zeros(SI,SJ);
MAT(1:size(P{1},1),1:size(P{1},2))=P{1};
MAT(size(P{1},1)+1:size(P{1},1)+size(Phat{1},1),1:size(Phat{1},2))=-Phat{1};
MAT(size(P{1},1)+1:size(P{1},1)+size(Phat{1},1),...
    size(Phat{1},2)+1:size(Phat{1},2)+size(A{1},2))=A{1};
indi=size(P{1},1)+size(Phat{1},1);
indj=0;
for i=2:k-1
    MAT(indi+1:indi+size(K{i-1},1),indj+1:indj+size(K{i-1},2))=-K{i-1};
    
    MAT(indi+1:indi+size(K{i-1},1),...
        indj+size(K{i-1},2)+1:indj+size(K{i-1},2)+size(L{i-1},2))=-L{i-1};
    
    MAT(indi+1:indi+size(K{i-1},1),indj+size(K{i-1},2)+size(L{i-1},2)+1:...
        indj+size(K{i-1},2)+size(L{i-1},2)+size(P{i},2))=P{i};
    
    MAT(indi+size(K{i-1},1)+1:indi+size(K{i-1},1)+size(Phat{i},1),...
        indj+size(K{i-1},2)+size(L{i-1},2)+1:...
        indj+size(K{i-1},2)+size(L{i-1},2)+size(Phat{i},2))=-Phat{i};
    
    MAT(indi+size(K{i-1},1)+1:indi+size(K{i-1},1)+size(Phat{i},1),...
        indj+size(K{i-1},2)+size(L{i-1},2)+size(Phat{i},2)+1:...
        indj+size(K{i-1},2)+size(L{i-1},2)+size(Phat{i},2)+size(A{i},2))=A{i};
    
    indi=indi+size(K{i-1},1)+size(Phat{i},1);
    indj=indj+size(K{i-1},2)+size(L{i-1},2);
end
%-------------------------------------------------------------------------%
VEC=zeros(size(MAT,1),1);
VEC(1:n*N{1},1)=U{1};
% VEC(n*N{1}+1:n*N{1}+size(Phat{1}*a{1},1),1)=Phat{1}*a{1};

mat=sum(abs(MAT),2);
MAT=MAT./repmat(mat,1,size(MAT,2));
VEC=VEC./mat;
sol=MAT\VEC;
ind=0;
hold all
for i=1:k-1
    u{i}=real(sol(ind+1:ind+n*N{i}));
    ind=ind+2*n*N{i};
%     plot3(X{i},Y{i},P{i}*u{i});
    err(i)=norm(U{i}-P{i}*u{i},inf)/norm(U{i},inf);
%     err(i)=norm(a{i}-u{i},inf)/norm(a{i},inf);
    uap{i}=P{i}*u{i};
end
% figure
plot(err,'-*')
