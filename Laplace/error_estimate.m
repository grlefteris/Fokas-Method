clear;
clc;
for ii=3:5
clearvars -except ii c3 c4 c5
kk=1;
for ee=3:20
   
uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=ii;
No=ee;
Np=No;

[X Y k hs x y hx hy w Nn]=domain(n,Np,1,0);
X=X{1};x=x{1};
Y=Y{1};y=y{1};

% [sides z ang m h]=irpolygon(n,Np,0,[0 0.8 0.5 1 -1 -1 0],[-1 -0.3 0.5 1 1 0 -1]);
% %[0 0.8 0.9 -0.85 -0.85 0],[-0.85 -0.3 0.9 1 0 -1]
% X=real(sides);X=X(:);
% Y=imag(sides);Y=Y(:);
% w=ang-(pi/2);

M=n*No;
R=2*M;
U=uan(X,Y);
W=kron(w,ones(Np,1));
q=(cos(W).*ux(X,Y))+(sin(W).*uy(X,Y));
P=kron(eye(n),Lbasis(No,Np));
[sides z ang m h]=irpolygon(n,Np,0,x,y);

%-------exact coefficients----%
% x=@(t,m,h) real(m+t*h);
% y=@(t,m,h) imag(m+t*h);
% fu=@(t,m,h) uan(x(t,m,h),y(t,m,h));
% fq=@(t,m,h,w)(cos(w).*ux(x(t,m,h),y(t,m,h)))+(sin(w).*uy(x(t,m,h),y(t,m,h)));
% bn=@(n,m,h)(n+0.5)*quad(@(t)fu(t,m,h).*lege(n,t),-1,1,1e-13);
% an=@(n,m,h,w)(n+0.5)*quad(@(t)fq(t,m,h,w).*lege(n,t),-1,1,1e-13);
% for j=1:n
%     for nn=0:No-1
%         bb((j-1)*No+(nn+1),1)=bn(nn,m(j),h(j));
%         aa((j-1)*No+(nn+1),1)=an(nn,m(j),h(j),w(j));
%     end
% end
G=LHS_Dirichlet(n,No,M,R,m,h);
H=Ptransform(n,No,M,R,m,h);
bb=P\U;
aa=P\q;
%-------exact coefficients----%

%--------mixed boundary conditions-------%
% temp1=G(:,[No+1:2*No 3*No+1:4*No]);
% G(:,[No+1:2*No 3*No+1:4*No])=-H(:,[No+1:2*No 3*No+1:4*No]);
% H(:,[No+1:2*No 3*No+1:4*No])=-temp1;
% t1=aa([No+1:2*No 3*No+1:s4*No]);
% aa([No+1:2*No 3*No+1:4*No])=bb([No+1:2*No 3*No+1:4*No]);
% % bb([No+1:2*No 3*No+1:4*No])=t1;
% temp1=G(:,[No+1:2*No 3*No+1:4*No 5*No+1:6*No]);
% G(:,[No+1:2*No 3*No+1:4*No 5*No+1:6*No])=-H(:,[No+1:2*No 3*No+1:4*No 5*No+1:6*No]);
% H(:,[No+1:2*No 3*No+1:4*No 5*No+1:6*No])=-temp1;
% t1=aa([No+1:2*No 3*No+1:4*No 5*No+1:6*No]);
% aa([No+1:2*No 3*No+1:4*No 5*No+1:6*No])=bb([No+1:2*No 3*No+1:4*No 5*No+1:6*No]);
% bb([No+1:2*No 3*No+1:4*No 5*No+1:6*No])=t1;
% s=sum(abs(G),2);
% G=G./repmat(s,1,size(G,2));
% H=H./repmat(s,1,size(H,2));
% sol=G\(H*bb);
% est(ee-2,1)=norm(G*sol-H*bb,inf);%error estimate
% ex(ee-2,1)= norm(sol-aa,inf)/norm(aa,inf);%exact error
%--------mixed boundary conditions-------%


%--------dirichlet boundary conditions-------%
s=sum(abs(G),2);
G=G./repmat(s,1,size(G,2));
H=H./repmat(s,1,size(H,2));
a=G\(H*bb);
% norm(G*a-H*bb,inf)
% norm(a-aa,inf)/norm(aa,inf)
% est(ee-2,1)=norm(G*a-H*bb,inf);
% ex(ee-2,1)= norm(a-aa,inf)/norm(aa,inf);
cn(kk,1)=norm(a-aa,inf)/norm(aa,inf);
kk=kk+1;
clearvars -except cn ii jj kk c3 c4 c5
%--------dirichlet boundary conditions-------%

end
% semilogy((3:13),est,'--ko','MarkerSize',10)
% hold all
% semilogy((3:13),ex,'-r*','MarkerSize',10)
% legend('Estimated','Exact')
% xlabel('N');ylabel('Error');
if ii==3
    c3=cn;
else if ii==4
    c4=cn;
else 
    c5=cn;
    end
end

end