clear;
clc;

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

% uan = @(x,y)  sinh(3*x).*sin(3*y);
% ux  = @(x,y)  3*cosh(3*x).*sin(3*y);
% uy  = @(x,y)  3*sinh(3*x).*cos(3*y);

n=4;
Nl=11;
Np=81;


[X Y k hs x y hx hy w N]=domain(n,Np,1,0);
k=k-1;
Nl=Nl*ones(k,1);Nl(end-6:end)=[9 8 7 6 5 4 3]';
ang=w+(pi/2);
for i=1:k-1
    M=n*Nl(i);
    R=2*M;
    U{i}=uan(X{i},Y{i});
    W{i}=kron(w,ones(N{i},1));
    q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
    hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
    Q{i}=restrict(n,N{i});
    P{i}=kron(eye(n),Lbasis(Nl(i),N{i}));
    [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});
    A{i}=LHS_Dirichlet(n,Nl(i),M,R,m{i},h{i});
    Phat{i}=Ptransform(n,Nl(i),M,R,m{i},h{i});
    a{i}=P{i}\U{i};
    b{i}=P{i}\q{i};
    sia(i)=n*N{i};
    sib(i)=2*n*M;
     sj(i)=n*Nl(i);
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
hs1=0.000000001;
hs2=0.000000002;
V=restrict(n,N{1});
hxx=kron(cos(w),ones(N{2},1));
hyy=kron(sin(w),ones(N{2},1));
Sx=@(x,y)-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));
Sy=@(x,y)-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));
Xn1=V*X{1}+((V*Sx(X{1},Y{1})).*abs(hs1*hxx));
Yn1=V*Y{1}+((V*Sy(X{1},Y{1})).*abs(hs1*hyy));
Xn2=V*X{1}+((V*Sx(X{1},Y{1})).*abs(hs2*hxx));
Yn2=V*Y{1}+((V*Sy(X{1},Y{1})).*abs(hs2*hyy));

AA=A{1};PP=Phat{1};
s=sum(abs(AA),2);
AA=AA./repmat(s,1,size(AA,2));
PP=PP./repmat(s,1,size(PP,2));
nvec=real(AA\(PP*a{1}));

h1=spdiags(-hs1*ones(n*N{1},1),0,n*N{1},n*N{1});
h2=spdiags(-hs2*ones(n*N{1},1),0,n*N{1},n*N{1});
Ut1=V*U{1}+V*(h1*P{1}*nvec);
Ut2=V*U{1}+V*(h2*P{1}*nvec);
norm(uan(Xn1,Yn1)-Ut1,inf)/norm(uan(Xn1,Yn1),inf)
norm(uan(Xn2,Yn2)-Ut2,inf)/norm(uan(Xn2,Yn2),inf)

or=7;
Pol=polyFit2D([U{1};Ut1;Ut2],[X{1};Xn1;Xn2],[Y{1};Yn1;Yn2],or,or);
Utt=polyVal2D(Pol,X{2},Y{2},or,or);
Uttt=polyVal2D(Pol,X{3},Y{3},or,or);
Utttt=polyVal2D(Pol,X{4},Y{4},or,or);
norm(U{2}-Utt,inf)/norm(U{2},inf)
norm(U{3}-Uttt,inf)/norm(U{3},inf)
norm(U{4}-Utttt,inf)/norm(U{4},inf)
a2=P{2}\Utt;
a3=P{3}\Uttt;
a4=P{4}\Utttt;
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
VEC(sia(1)+sib(1)+1:sia(1)+sib(1)+sia(2),1)=Utt;
VEC(sia(1)+sib(1)+sia(2)+sib(2)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3),1)=Uttt;
VEC(sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+1:sia(1)+sib(1)+sia(2)+sib(2)+sia(3)+sib(3)+sia(4),1)=Utttt;

mat=sum(abs(MAT),2);
MAT=MAT./repmat(mat,1,size(MAT,2));
VEC=VEC./mat;
sol=(MAT)\(VEC);

% K=sparse(MAT);
% tic;P=colamd(K);K2=K(:,P);R2=qr(K2,0);xx2=R2\(R2'\(K2'*VEC));xx2(P)=xx2;toc
% sol=xx2;

ind=0;
hold all
for i=1:k-1
    u{i}=real(sol(ind+1:ind+n*Nl(i)));
    un{i}=real(sol(ind+n*Nl(i)+1:ind+2*n*Nl(i)));
%     plot3(X{i},Y{i},P{i}*u{i});
%     mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),reshape(P{i}*u{i},n,N{i}))
%      mesh(reshape(X{i},n,N{i}),reshape(Y{i},n,N{i}),abs(reshape(U{i},n,N{i})-reshape(P{i}*u{i},n,N{i}))/norm(U{i},inf));
    ind=ind+2*n*Nl(i);
%     err(i)=norm(U{i}-P{i}*u{i},inf)/norm(U{i},inf);
    err(i)=norm(a{i}-u{i},inf)/norm(a{i},inf);
    errn(i)=norm(b{i}-un{i},inf)/norm(b{i},inf);
%     uap{i}=P{i}*u{i};
end
semilogy(2:k-1,err(2:end),'-k*','MarkerSize',10)
semilogy(2:k-1,errn(2:end),'-ko','MarkerSize',10)
xlabel('Levels','fontsize',16,'FontName','Times New Roman');
ylabel('Maximum Error','fontsize',16,'FontName','Times New Roman');