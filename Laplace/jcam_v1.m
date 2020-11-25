clear;
clc;
warning off; % To avoid unnecessary warnings related to high condition number

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=4;
Nl=11;
Np=501;

[X Y k hs x y hx hy w N]=domain(n,Np,1,0);
% [sides z ang m h]=irpolygon(n,N,0,[0 1 1 0],[0 0 1 1]);
% X=real(sides);X=X(:);
% Y=imag(sides);Y=Y(:);
% w=ang-(pi/2);

M=n*Nl;
R=2*M;
U{1}=uan(X{1},Y{1});
W{1}=kron(w,ones(N{1},1));
q{1}=(cos(W{1}).*ux(X{1},Y{1}))+(sin(W{1}).*uy(X{1},Y{1}));
hk{1}=spdiags(-hs*ones(n*N{1},1),0,n*N{1},n*N{1});
V{1}=restrict(n,N{1});
P{1}=kron(eye(n),Lbasis(Nl,N{1}));
[sides z ang m{1} h{1}]=irpolygon(n,N{1},0,x{1},y{1});
G{1}=PLHS_Dirichlet(n,Nl,M,R,m{1},h{1},1);
H{1}=Ptransform(n,Nl,M,R,m{1},h{1});
b{1}=P{1}\U{1};
a{1}=P{1}\q{1};
s{1}=sum(abs(G{1}),2);
G{1}=G{1}./repmat(s{1},1,size(G{1},2));
H{1}=H{1}./repmat(s{1},1,size(H{1},2));

u{1}=U{1};
bb{1}=real(b{1});
% aa{1}=real(G{1}\(H{1}*bb{1}));
Nl=Nl*ones(k,1);%Nl(end-6:end)=[10 9 8 7 6 5 4]';


rhs=mRHS1(n,Np,M,R,m{1},h{1},u{1});
rhs=rhs./s{1};
aa{1}=real(G{1}\(rhs));


for i=2:k-1
    M=n*Nl(i);
    R=2*M;
    U{i}=uan(X{i},Y{i});
    W{i}=kron(w,ones(N{i},1));
    q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
    hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
    V{i}=restrict(n,N{i});
    P{i}=kron(eye(n),Lbasis(Nl(i),N{i}));
    [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});
    G{i}=PLHS_Dirichlet(n,Nl(i),M,R,m{i},h{i},1);
    H{i}=Ptransform(n,Nl(i),M,R,m{i},h{i});
    b{i}=P{i}\U{i};
    a{i}=P{i}\q{i};
    s{i}=sum(abs(G{i}),2);
    G{i}=G{i}./repmat(s{i},1,size(G{i},2));
    H{i}=H{i}./repmat(s{i},1,size(H{i},2));
    
  
    u{i}=V{i-1}*u{i-1}+V{i-1}*(hk{i-1}*P{i-1}*aa{i-1});
%     bb{i}=P{i}\u{i};
%     aa{i}=real(G{i}\(H{i}*bb{i}));
    rhs=mRHS1(n,N{i},M,R,m{i},h{i},u{i});
    rhs=rhs./s{i};
    aa{i}=real(G{i}\(rhs));
    err(i)=norm(U{i}-u{i},inf)/norm(U{i},inf);
end



