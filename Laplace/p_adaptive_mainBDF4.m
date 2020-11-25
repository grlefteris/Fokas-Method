clear;
clc;

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=4;
Nn=11;

N1=11;
N2=11;
N3=11;
N4=11;


[X Y k hs x y hx hy w N]=domain(n,Nn,1,0);clear X Y
[X1 Y1 k1 hs xx yy hx hy w NN]=domain(n,N1,1,0);
[X2 Y2 k2 hs xx yy hx hy w NN]=domain(n,N2,1,0);
[X3 Y3 k3 hs xx yy hx hy w NN]=domain(n,N3,1,0);
[X4 Y4 k4 hs xx yy hx hy w NN]=domain(n,N4,1,0);
% X{1}(1:N1)=X1{1}(1:N1);
% X{1}(N1+1:N1+N2)=X2{1}(N2+1:2*N2);
% X{1}(N1+N2+1:N1+N2+N3)=X3{1}(2*N3+1:3*N3);
% X{1}(N1+N2+N3+1:N1+N2+N3+N4)=X4{1}(3*N4+1:4*N4);
% Y{1}(1:N1)=Y1{1}(1:N1);
% Y{1}(N1+1:N1+N2)=Y2{1}(N2+1:2*N2);
% Y{1}(N1+N2+1:N1+N2+N3)=Y3{1}(2*N3+1:3*N3);
% Y{1}(N1+N2+N3+1:N1+N2+N3+N4)=Y4{1}(3*N4+1:4*N4);
% X{1}=X{1}';Y{1}=Y{1}';
% k=k4;
for i=1:k-1
    N1=Nn-2*(i-1);
    N2=Nn-2*(i-1);
    N3=Nn-2*(i-1);
    N4=Nn-2*(i-1);
%     
    X{i}(1:N1)=X1{1}(1:N1);
    X{i}(N1+1:N1+N2)=X2{1}(N2+1:2*N2);
    X{i}(N1+N2+1:N1+N2+N3)=X3{1}(2*N3+1:3*N3);
    X{i}(N1+N2+N3+1:N1+N2+N3+N4)=X4{1}(3*N4+1:4*N4);
    Y{i}(1:N1)=Y1{1}(1:N1);
    Y{i}(N1+1:N1+N2)=Y2{1}(N2+1:2*N2);
    Y{i}(N1+N2+1:N1+N2+N3)=Y3{1}(2*N3+1:3*N3);
    Y{i}(N1+N2+N3+1:N1+N2+N3+N4)=Y4{1}(3*N4+1:4*N4);
    X{i}=X{i}';Y{i}=Y{i}';   
    t1=linspace(-1,1,N1)';t2=linspace(-1,1,N2)';t3=linspace(-1,1,N3)';t4=linspace(-1,1,N4)';
%     r1=[-1 -0.9 0 0.9 1]';r2=[-1 -0.9 0 0.9 1]';r3=[-1 -0.9 0 0.9 1]';r4=[-1 -0.9 0 0.9 1]';
    M=n*Nn;
    R=2*M;
    U{i}=uan(X{i},Y{i});
    W{i}=[w(1)*ones(N1,1);w(2)*ones(N2,1);w(3)*ones(N3,1);w(4)*ones(N4,1)];
%     q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
    hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
    Q{i}=restrict(n,N{i},2);
    P{i}=blkdiag(Lbasis(N1,N1,t1),Lbasis(N2,N2,t2),Lbasis(N3,N3,t3),Lbasis(N4,N4,t4));
%     p{i}=blkdiag(Lbasis(N1,5,r1),Lbasis(N2,5,r2),Lbasis(N3,5,r3),Lbasis(N4,5,r3));
    [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});
%     A{i}=LHS_Dirichlet_TEST(n,N1,N2,N3,N4,M,R,m{i},h{i});
%     Phat{i}=Ptransform_TEST(n,N1,N2,N3,N4,M,R,m{i},h{i});
    A{i}=LHS_Dirichlet_TEST(n,Nn,Nn,Nn,Nn,M,R,m{i},h{i});
    Phat{i}=Ptransform_TEST(n,Nn,Nn,Nn,Nn,M,R,m{i},h{i});
    a{1}=P{1}\U{1};
    %     a{i}=P{i}\U{i};
%     b{i}=P{i}\q{i};

    sia(i)=n*N{i};
    sib(i)=2*n*M;
     sj(i)=n*Nn;

    MAT{i}=blkdiag(P{i},A{i});
end
for i=1:k-2
    ma{i}=blkdiag(-Phat{i},zeros(sia(i+1),sj(i)));
end
ma=blkdiag(ma{:},-Phat{i+1});
MAT=blkdiag(MAT{:});
MAT(sia(1)+1:end,1:end-sj(end))=MAT(sia(1)+1:end,1:end-sj(end))+ma;
% MAT=zeros(sum(sia)+sum(sib),2*sum(sj));


%---------------------------coefficient matrix assembly-------------------%
MAT(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),1:sj(1))=-Q{1}*P{1};
% MAT(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),sj(1)+1:2*sj(1))=-Q{1}*(hk{1}*P{1});%explicit euler
MAT(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),2*sj(1)+sj(2)+1:2*sj(1)+2*sj(2))=-hk{2}*P{2};%implicit euler

MAT(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),1:sj(1))=(1/3)*Q{2}*Q{1}*P{1};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),2*sj(1)+1:2*sj(1)+sj(2))=-(4/3)*Q{2}*P{2};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),2*sj(1)+2*sj(2)+sj(3)+1:2*sj(1)+2*sj(2)+2*sj(3))=-(2/3)*hk{3}*P{3};

MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),1:sj(1))=-(2/11)*Q{3}*Q{2}*Q{1}*P{1};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),2*sj(1)+1:2*sj(1)+sj(2))=(9/11)*Q{3}*Q{2}*P{2};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),2*sj(1)+2*sj(2)+1:2*sj(1)+2*sj(2)+sj(3))=-(18/11)*Q{3}*P{3};
MAT(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),2*sj(1)+2*sj(2)+2*sj(3)+sj(4)+1:2*sj(1)+2*sj(2)+2*sj(3)+2*sj(4))=-(6/11)*hk{4}*P{4};

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

mat=sum(abs(MAT),2);
MAT=MAT./repmat(mat,1,size(MAT,2));
VEC=VEC./mat;
sol=MAT\VEC;
ind=0;
hold all
for i=1:k-1
    u{i}=real(sol(ind+1:ind+n*Nn));
    un{i}=real(sol(ind+n*N{i}+1:ind+2*n*N{i}));
    plot3(X{i},Y{i},P{i}*u{i},'-k.','MarkerSize',16);
%     mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),reshape(P{i}*u{i},n,N{i}))
%      mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),abs(reshape(U{i},n,N{i})-reshape(P{i}*u{i},n,N{i}))/norm(U{i},inf));
    ind=ind+2*n*N{i};
    err(i)=norm(U{i}-P{i}*u{i},inf)/norm(U{i},inf);
%     err(i)=norm(a{i}-u{i},inf)/norm(a{i},inf);
%     err(i)=norm(b{i}-un{i},inf)/norm(b{i},inf);
%     uap{i}=P{i}*u{i};
end
% semilogy(2:k-1,err(2:end),'-ksquare','MarkerSize',12)
xlabel('Levels','fontsize',12);ylabel('Maximum Error','fontsize',12);