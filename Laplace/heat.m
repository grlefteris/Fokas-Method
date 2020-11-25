clear;
clc;

% xn=@(x,y) (x*cos(pi/4)-y*sin(pi/4));
% yn=@(x,y) (x*sin(pi/4)+y*cos(pi/4));
% 
% uan= @(x,y)      exp(-2*(pi^2)*0.2).*sin(pi*xn(x,y)).*sin(pi*yn(x,y));
% ux = @(x,y)   pi*exp(-2*(pi^2)*0.2).*cos(pi*xn(x,y)).*sin(pi*yn(x,y));
% uy = @(x,y)   pi*exp(-2*(pi^2)*0.2).*sin(pi*xn(x,y)).*cos(pi*yn(x,y));

uan= @(x,y)      exp(-2*(pi^2)*0.2).*sin(pi*x).*sin(pi*y);
ux = @(x,y)   pi*exp(-2*(pi^2)*0.2).*cos(pi*x).*sin(pi*y);
uy = @(x,y)   pi*exp(-2*(pi^2)*0.2).*sin(pi*x).*cos(pi*y);

n=4;
Nn=21;

[X Y k hs x y hx hy w N]=domain(n,Nn,1,0);
for i=1:k-1
    M=n*N{i};
    R=2*M;
    U{i}=uan(X{i},Y{i});
    W{i}=kron(w,ones(N{i},1));
    q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
    hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
    Q{i}=restrict(n,N{i},2);
    P{i}=kron(eye(n),leg_basis(N{i}));
    [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});
    A{i}=LHS_Dirichlet(n,N{i},M,R,m{i},h{i});
    Phat{i}=Ptransform(n,N{i},M,R,m{i},h{i});
    a{i}=P{i}\U{i};
    b{i}=P{i}\q{i};
    sia(i)=n*N{i};
    sib(i)=2*n*M;
     sj(i)=n*N{i};
    MAT{i}=blkdiag(P{i},A{i});
end
for i=1:k-2
    ma{i}=blkdiag(-Phat{i},zeros(sia(i+1),sj(i)));
end
ma=blkdiag(ma{:},-Phat{i+1});
MAT=blkdiag(MAT{:});
MAT(sia(1)+1:end,1:end-sj(end))=MAT(sia(1)+1:end,1:end-sj(end))+ma;
% MAT=zeros(sum(sia)+sum(sib),2*sum(sj));


%-------------Pre-Processing Interpolation, Initialize BDF4---------------%
Sx=@(x,y)-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));
Sy=@(x,y)-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));

hs=0.000000001;
V2=restrict(n,N{1},1);
V3=restrict(n,N{2},1);
V4=restrict(n,N{3},1);
hxx2=kron(cos(w),ones(N{2},1));
hyy2=kron(sin(w),ones(N{2},1));
hxx3=kron(cos(w),ones(N{3},1));
hyy3=kron(sin(w),ones(N{3},1));
hxx4=kron(cos(w),ones(N{4},1));
hyy4=kron(sin(w),ones(N{4},1));


X2=V2*X{1}+((V2*Sx(X{1},Y{1})).*abs(hs*hxx2));
Y2=V2*Y{1}+((V2*Sy(X{1},Y{1})).*abs(hs*hyy2));
X3=V3*X2+((V3*Sx(X2,Y2)).*abs(hs*hxx3));
Y3=V3*Y2+((V3*Sy(X2,Y2)).*abs(hs*hyy3));
X4=V4*X3+((V4*Sx(X3,Y3)).*abs(hs*hxx4));
Y4=V4*Y3+((V4*Sy(X3,Y3)).*abs(hs*hyy4));

AA=A{1};PP=Phat{1};
s=sum(abs(AA),2);
AA=AA./repmat(s,1,size(AA,2));
PP=PP./repmat(s,1,size(PP,2));
nvec=real(AA\(PP*a{1}));

h2=spdiags(-1*hs*ones(n*N{1},1),0,n*N{1},n*N{1});
h3=spdiags(-2*hs*ones(n*N{2},1),0,n*N{2},n*N{2});
h4=spdiags(-3*hs*ones(n*N{3},1),0,n*N{3},n*N{3});

U2=V2*U{1}+V2*(h2*P{1}*nvec);
U3=V3*V2*U{1}+V3*(h3*V2*(P{1}*nvec));
U4=V4*V3*V2*U{1}+V4*(h4*V3*(V2*(P{1}*nvec)));
norm(uan(X2,Y2)-U2,inf)/norm(uan(X2,Y2),inf)
norm(uan(X3,Y3)-U3,inf)/norm(uan(X3,Y3),inf)
norm(uan(X4,Y4)-U4,inf)/norm(uan(X3,Y3),inf)


a2=P{2}\U2;
a3=P{3}\U3;
a4=P{4}\U4;
a{2}=P{2}\uan(X2,Y2);
a{3}=P{3}\uan(X3,Y3);
a{4}=P{4}\uan(X4,Y4);
%-------------Pre-Processing Interpolation, Initialize BDF4---------------%


%---------------------------coefficient matrix assembly-------------------%
% MAT(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),1:sj(1))=-Q{1}*P{1};
% % MAT(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),sj(1)+1:2*sj(1))=-Q{1}*(hk{1}*P{1});%explicit euler
% MAT(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),2*sj(1)+sj(2)+1:2*sj(1)+2*sj(2))=-hk{2}*P{2};%implicit euler

% MAT(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),1:sj(1))=(1/3)*Q{2}*Q{1}*P{1};
% MAT(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),2*sj(1)+1:2*sj(1)+sj(2))=-(4/3)*Q{2}*P{2};
% MAT(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),2*sj(1)+2*sj(2)+sj(3)+1:2*sj(1)+2*sj(2)+2*sj(3))=-(2/3)*hk{3}*P{3};

% MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),1:sj(1))=-(2/11)*Q{3}*Q{2}*Q{1}*P{1};
% MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),2*sj(1)+1:2*sj(1)+sj(2))=(9/11)*Q{3}*Q{2}*P{2};
% MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),2*sj(1)+2*sj(2)+1:2*sj(1)+2*sj(2)+sj(3))=-(18/11)*Q{3}*P{3};
% MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),2*sj(1)+2*sj(2)+2*sj(3)+sj(4)+1:2*sj(1)+2*sj(2)+2*sj(3)+2*sj(4))=-(6/11)*hk{4}*P{4};

MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+sia(5),1:sj(1))=(3/25)*Q{4}*Q{3}*Q{2}*Q{1}*P{1};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+sia(5),2*sj(1)+1:2*sj(1)+sj(2))=-(16/25)*Q{4}*Q{3}*Q{2}*P{2};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+sia(5),2*sj(1)+2*sj(2)+1:2*sj(1)+2*sj(2)+sj(3))=(36/25)*Q{4}*Q{3}*P{3};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+sia(5),2*sj(1)+2*sj(2)+2*sj(3)+1:2*sj(1)+2*sj(2)+2*sj(3)+sj(4))=-(48/25)*Q{4}*P{4};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4)+sib(4)+sia(5),2*sj(1)+2*sj(2)+2*sj(3)+2*sj(4)+sj(5)+1:2*sj(1)+2*sj(2)+2*sj(3)+2*sj(4)+2*sj(5))=-(12/25)*hk{5}*P{5};

indi=sum(sia(1:5))+sum(sib(1:5));
indj=2*sj(1);
for i=6:k-1
    
    MAT(indi+1:indi+sia(i),indj+1:indj+sj(i-4))=(3/25)*Q{i-1}*Q{i-2}*Q{i-3}*Q{i-4}*P{i-4};
    MAT(indi+1:indi+sia(i),indj+2*sj(i-4)+1:indj+2*sj(i-4)+sj(i-3))=-(16/25)*Q{i-1}*Q{i-2}*Q{i-3}*P{i-3};
    MAT(indi+1:indi+sia(i),indj+2*sj(i-4)+2*sj(i-3)+1:indj+2*sj(i-4)+2*sj(i-3)+sj(i-2))=(36/25)*Q{i-1}*Q{i-2}*P{i-2};
    MAT(indi+1:indi+sia(i),indj+2*sj(i-4)+2*sj(i-3)+2*sj(i-2)+1:indj+2*sj(i-4)+2*sj(i-3)+2*sj(i-2)+sj(i-1))=-(48/25)*Q{i-1}*P{i-1};
    MAT(indi+1:indi+sia(i),indj+2*sj(i-4)+2*sj(i-3)+2*sj(i-2)+2*sj(i-1)+sj(i)+1:indj+2*sj(i-4)+2*sj(i-3)+2*sj(i-2)+2*sj(i-1)+2*sj(i))=-(12/25)*hk{i}*P{i};
    
    indi=indi+sia(i)+sib(i);
    indj=indj+2*sj(i-4);
end
%-------------------------------------------------------------------------%
VEC=zeros(size(MAT,1),1);
VEC(1:n*N{1},1)=U{1};
VEC(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),1)=U2;
VEC(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),1)=U3;
VEC(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),1)=U4;

mat=sum(abs(MAT),2);
MAT=MAT./repmat(mat,1,size(MAT,2));
VEC=VEC./mat;
sol=MAT\VEC;
ind=0;
hold all
X{2}=X2;
X{3}=X3;
X{4}=X4;
Y{2}=Y2;
Y{3}=Y3;
Y{4}=Y4;
for i=1:k-1
    u{i}=real(sol(ind+1:ind+n*N{i}));
%     plot3(X{i},Y{i},P{i}*u{i});
%     mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),reshape(P{i}*u{i},n,N{i}))
%      mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),abs(reshape(U{i},n,N{i})-reshape(P{i}*u{i},n,N{i}))/norm(U{i},inf));
    ind=ind+2*n*N{i};
%     err(i)=norm(U{i}-P{i}*u{i},inf)/norm(U{i},inf);
    err(i)=norm(a{i}-u{i},inf);
    uap{i}=P{i}*u{i};
end
plot(1:k-1,err,'-k.','MarkerSize',10)
xlabel('Levels','fontsize',12);ylabel('Maximum Error','fontsize',12);