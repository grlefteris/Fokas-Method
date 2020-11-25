clear;
clc;

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=4;
N=10;


[sides z1 ang m h]=irpolygon(n,N,0,[0 1 1 0 0],[0 0 1 1 0]);


X=reshape(real(sides),n*N,1);
Y=reshape(imag(sides),n*N,1);


w=ang-(pi/2);


M=2;
R=4;
M=n*N;
R=0.25*M;

k=2;
for i=1:k-1
U=uan(X,Y);
W=kron(w,ones(N,1));
q=(cos(W).*ux(X,Y))+(sin(W).*uy(X,Y));
end
% U(N+1:2*N)=ones(N,1);%%%dummy solution at side 2

P=kron(eye(n),leg_basis(N));

G=LHS_Dirichlet(n,N,M,R,m,h);
H=Ptransform(n,N,M,R,m,h);
a1=P\U;
b1=P\q;

[rhs]=mRHS4(n,N,M,R,m,h,reshape(U,N,n));

% %impose mixed boundary conditions%
% bcD1=H(:,1:N);
% bcD3=H(:,2*N+1:3*N);
% bcN1=G(:,1:N);
% bcN3=G(:,2*N+1:3*N);
% G(:,1:N)=-bcD1;
% G(:,2*N+1:3*N)=-bcD3;
% H(:,1:N)=-bcN1;
% H(:,2*N+1:3*N)=-bcN3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
s1=sum(abs(G),2);
G=G./repmat(s1,1,size(G,2));
H=H./repmat(s1,1,size(H,2));
rhs=rhs./s1;

% c=[a1(1:N);b1(N+1:2*N);a1(2*N+1:3*N);b1(3*N+1:4*N)];
% d=[b1(1:N);a1(N+1:2*N);b1(2*N+1:3*N);a1(3*N+1:4*N)];
norm(G*b1-(H*a1),inf)
sol=real(G\(H*a1));
% norm(sol-b1,inf)/norm(b1,inf)
norm(sol-b1,inf)/norm(b1,inf)

% p=leg_basis(N);
% for i=1:n
%     ec(i,1)=norm(sol((i-1)*N+1:i*N)-c((i-1)*N+1:i*N),inf);
%     es(i,1)=norm(p*sol((i-1)*N+1:i*N)-p*c((i-1)*N+1:i*N),inf);
% end




% %------adaptive mesh refinement interior levels------------%
% k=k-1;
% hs1=0.001;
% V=restrict(n,N,1);
% Pl=kron(eye(n),leg_basis(N-2));
% hxx=kron(cos(w),ones(N-2,1));
% hyy=kron(sin(w),ones(N-2,1));
% % Sx=@(x,y)-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));%-((x==0)&(y>0))+((x==0)&(y<0));
% % Sy=@(x,y)-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));%-((x>0)&(y==0))+((x<0)&(y==0));
% Sx=@(x,y)  -(x==max(x)) +(x==min(x));
% Sy=@(x,y) ((x>min(x))&(x<max(x))&(y<0.5))-((x>min(x))&(x<max(x))&(y>0.5));
% 
% Xe=X;Ye=Y;Ue=U;UT1=U;
% nvec=real(G\(H*d+v1));
% 
% c=[a1(1:N);b1(N+1:2*N);a1(2*N+1:3*N);b1(3*N+1:4*N)];
% 
% nvec(1:N)=b1(1:N);
% nvec(2*N+1:3*N)=b1(2*N+1:3*N);
% 
% % hold all
% for i=1:3
% 
%     Xt{i}=V*X+((V*Sx(X,Y)).*abs(i*hs1*hxx));
%     Yt{i}=V*Y+((V*Sy(X,Y)).*abs(i*hs1*hyy));
% 
%     Xe=[Xe;Xt{i}];
%     Ye=[Ye;Yt{i}];
%     
%     h1{i}=spdiags(-i*hs1*ones(n*N,1),0,n*N,n*N);
%     
%     Ut1{i}=V*U+V*(h1{i}*P*nvec);
% %     plot3(Xt{i},Yt{i},Ut1{i},'-b')
%     
%     UT1=[UT1;Ut1{i}];
%     
%     ee(i)=norm(uan(Xt{i},Yt{i})-Ut1{i},inf);
%     
% end
% 
% y1=@(x) x+1;
% 
% NN=10000;
% xx=rand(NN,1);
% yy=2*rand(NN,1);
% xn=xx( (yy<y1(xx)) );
% yn=yy( (yy<y1(xx)) );
% 
% o=3;
% Pol=polyFit2D(UT1,Xe,Ye,o,o);
% UP=polyVal2D(Pol,xn,yn,o,o);
% norm(uan(xn,yn)-UP,inf)%/norm(uan(xn,yn),inf)
% % plot3(xn,yn,abs(uan(xn,yn)-UP)/norm(uan(xn,yn),inf),'.','MarkerSize',5)
% % figure
% % UP=abs(uan(xn,yn)-UP)/norm(uan(xn,yn),inf);
% 
% rbf=rbfcreate([Xe';Ye'],UT1','RBFFunction', 'multiquadric');
% zm =rbfinterp([xn';yn'],rbf);
% 
% 
% % F = TriScatteredInterp(Xe,Ye,UT1,'linear');
% % % zm=interpn(Xe,Ye,UT1,xn,yn,'cubic');
% % zm=F(xn,yn);
% norm(uan(xn,yn)-zm',inf)/norm(uan(xn,yn),inf)
% plot3(xn,yn,abs(uan(xn,yn)-zm')/norm(uan(xn,yn),inf),'.','MarkerSize',5)
% 
% 
% % mesh(xm,ym,zm);
% % 
% % xp=xm( ((ym>=y2(xm))&(ym<=y3(xm))) );%|  ((ym>y6(xm))&(ym>y1(xm))&(ym<yg(xm))) );
% % yp=ym( ((ym>=y2(xm))&(ym<=y3(xm))) );%|  ((ym>y6(xm))&(ym>y1(xm))&(ym<yg(xm))) );
% % zp=zm( ((ym>=y2(xm))&(ym<=y3(xm))) );
% % hold all
% % plot3(xp,yp,zp,'w.','MarkerSize',10)
% % grid off
