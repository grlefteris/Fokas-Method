clear;
clc;

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=30;
Nl=5;
Np=21;

M=n*Nl/30;
R=2*M;

% n=4;
% Nn=11;
% M=n*Nn/4;
% R=2*M;

[X Y k hs x y hx hy w N]=domain(n,Np,1,0);
k=2;
for i=1:k-1
U{i}=uan(X{i},Y{i});
W{i}=kron(w,ones(N{i},1));
q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
hx{i}=spdiags(hx{i},0,n*N{i},n*N{i});
hy{i}=spdiags(hy{i},0,n*N{i},n*N{i});
hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
end
[sides z ang m{1} h{1}]=irpolygon(n,N{1},0,x{1},y{1});
Q{1}=restrict(n,N{1});
P{1}=kron(eye(n),Lbasis(Nl,N{1}));
A{1}=PLHS_Dirichlet(n,Nl,M,R,m{1},h{1},1);
Phat{1}=Ptransform(n,Nl,M,R,m{1},h{1});
K{1}=Q{1}*P{1};
L{1}=Q{1}*hk{1}*P{1};
a{1}=P{1}\U{1};
b{1}=P{1}\q{1};
sj(1)=n*Nl;
for i=2:k-1
    M=n*N{i};
    R=2*M;
    [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});    
    A{i}=PLHS_Dirichlet(n,Nl,M,R,m{i},h{i},1);
    Phat{i}=Ptransform(n,Nl,M,R,m{i},h{i});
    P{i}=kron(eye(n),Lbasis(Nl,N{i}));  
    Q{i}=restrict(n,N{i});
    K{i}=Q{i}*P{i};
    L{i}=hk{i}*P{i};
    a{i}=P{i}\U{i};
    b{i}=P{i}\q{i};
    sj(i)=n*Nl;
end
%-----------coefficient matrix size-----%
SI=0;
SJ=0;
for i=1:k-1
    SI=SI+size(P{i},1)+size(A{i},1);
    SJ=SJ+(2*n*Nl);
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
        
    MAT(indi+1:indi+size(K{i-1},1),indj+size(K{i-1},2)+size(L{i-1},2)+1:...
        indj+size(K{i-1},2)+size(L{i-1},2)+size(P{i},2))=P{i};
            
    MAT(indi+1:indi+size(K{i-1},1),...
        indj+2*size(K{i-1},2)+size(P{i},2)+1:indj+2*size(K{i-1},2)+size(P{i},2)+size(L{i},2))=-L{i};           
    
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
mat=sum(abs(MAT),2);
MAT=MAT./repmat(mat,1,size(MAT,2));
VEC=VEC./mat;

%------adaptive mesh refinement interior levels------------%
k=k-1;
hs1=0.005;
% hs1=0.01;
V=restrict(n,N{k});
Pl=kron(eye(n),Lbasis(Nl,N{k+1})); %kron(eye(n),leg_basis(N{k+1}));
hxx=kron(cos(w),ones(N{k+1},1));
hyy=kron(sin(w),ones(N{k+1},1));
Sx=@(x,y)-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));
Sy=@(x,y)-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));

jt=size(MAT,2);
for i=1:(hs/hs1)
    ii=size(MAT,1);
    jj=size(MAT,2);
    Xt{i}=V*X{k}+((V*Sx(X{k},Y{k})).*abs(i*hs1*hxx));
    Yt{i}=V*Y{k}+((V*Sy(X{k},Y{k})).*abs(i*hs1*hyy));
    h1{i}=spdiags(-i*hs1*ones(n*N{k},1),0,n*N{k},n*N{k});
    mati=zeros(n*N{k+1},jj+n*N{k+1});
    matj=zeros(ii,n*N{k+1});
    MAT=[MAT matj;mati];

    MAT(ii+1:ii+n*N{k+1},jj+1:jj+n*Nl)=Pl;
    MAT(ii+1:ii+n*N{k+1},jt-2*n*Nl+1:jt-n*Nl)=-V*P{k};
    MAT(ii+1:ii+n*N{k+1},jt-n*Nl+1:jt)=-V*(h1{i}*P{k});
    
    VEC=[VEC;zeros(n*N{k+1},1)];
end
k=k+1;
%------adaptive mesh refinement interior levels------------%

sol=MAT\VEC;
% K=sparse(MAT);
% tic;Pc=colamd(K);K2=K(:,Pc);R2=qr(K2,0);xx2=R2\(R2'\(K2'*VEC));xx2(Pc)=xx2;toc
% sol=xx2;
% K=sparse(MAT);
% sol=lsqr(K,VEC,1e-8);

ind=0;
ind2=2*sum(sj);
hold all
for i=1:k-1
    u{i}=real(sol(ind+1:ind+n*Nl));
    un{i}=real(sol(ind+n*N{i}+1:ind+2*n*Nl));
    ind=ind+2*n*Nl;
%     plot3(X{i},Y{i},P{i}*u{i},'-k*');
%     mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),reshape(P{i}*u{i},n,N{i}))
    mesh(reshape(X{i},N{i},n),reshape(Y{i},N{i},n),abs(reshape(U{i},N{i},n)-reshape(P{i}*u{i},N{i},n))/norm(U{i},inf));
    err(i)=norm(U{i}-P{i}*u{i},inf)/norm(U{i},inf);
    err(i)=norm(a{i}-u{i},inf)/norm(a{i},inf);
%     err(i)=norm(b{i}-un{i},inf)/norm(b{i},inf);
end
for i=1:(hs/hs1)
            
    ut{i}=real(sol(ind2+1:ind2+n*Nl));
    
    ind2=ind2+n*Nl;
    
    XT=reshape(Xt{i},N{k},n);YT=reshape(Yt{i},N{k},n);UT=reshape(Pl*ut{i},N{k},n);UA=reshape(uan(Xt{i},Yt{i}),N{k},n);
%     tx=XT(1,2:end);tx=[tx XT(1,1)];XT=[XT;tx];
%     ty=YT(1,2:end);ty=[ty YT(1,1)];YT=[YT;ty];
%     tu=UT(1,2:end);tu=[tu UT(1,1)];UT=[UT;tu];
%     ta=UA(1,2:end);ta=[ta UA(1,1)];UA=[UA;ta];
%      mesh(XT,YT,UT);
    
    
%     mesh(reshape(Xt{i},n,N{k}),reshape(Yt{i},n,N{k}),reshape(Pl*ut{i},n,N{k}))
    mesh(XT,YT,abs(UT-UA)./norm(ut{i},inf));
%     plot3(XT,YT,UT,'k*')
    
    
%     mesh(reshape(Xt{i},N{k},n),reshape(Yt{i},N{k},n),...
%     abs(reshape(Pl*ut{i},N{k},n)-reshape(uan(Xt{i},Yt{i}),N{k},n))/norm(uan(Xt{i},Yt{i}),inf));
%     plot3(Xt{i},Yt{i},Pl*ut{i},'-r.');
    ee(i)=norm(uan(Xt{i},Yt{i})-Pl*ut{i},inf)/norm(uan(Xt{i},Yt{i}),inf);
end

% semilogy(2:k-1,err(2:end),'-k*','MarkerSize',12)
% xlabel('Levels','fontsize',16,'FontName','Times New Roman');
% ylabel('Maximum Error','fontsize',16,'FontName','Times New Roman');
