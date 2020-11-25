clear;
clc;

% for ii=1:20

uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);

n=4;
Nn=201;
% M=n*Nn;
% R=2*M;
M=3;R=5;

[X Y k hs x y hx hy w N]=domain(n,Nn,1,0);

% for i=1:k-1
% U{i}=uan(X{i},Y{i});
% W{i}=kron(w,ones(N{i},1));
% q{i}=(cos(W{i}).*ux(X{i},Y{i}))+(sin(W{i}).*uy(X{i},Y{i}));
% hx{i}=spdiags(hx{i},0,n*N{i},n*N{i});
% hy{i}=spdiags(hy{i},0,n*N{i},n*N{i});
% hk{i}=spdiags(-hs*ones(n*N{i},1),0,n*N{i},n*N{i});
% Q{i}=restrict(n,N{i},2);
% P{i}=kron(eye(n),leg_basis(N{i}));
% end
% 
% [sides z ang m{1} h{1}]=irpolygon(n,N{1},0,x{1},y{1});
% 
% 
% A{1}=PLHS_Dirichlet(n,N{1},M,R,m{1},h{1},1);
% Phat{1}=Ptransform(n,N{1},M,R,m{1},h{1});
% a{1}=P{1}\U{1};
% b{1}=P{1}\q{1};
% for i=2:k-1
%     M=n*N{i};
%     R=6*M;
% %     M=3;R=5;
%     [sides z ang m{i} h{i}]=irpolygon(n,N{i},0,x{i},y{i});
%     
%     A{i}=PLHS_Dirichlet(n,N{i},M,R,m{i},h{i},1);
%     Phat{i}=Ptransform(n,N{i},M,R,m{i},h{i});
%     a{i}=P{i}\U{i};
%     b{i}=P{i}\q{i};
% end
% 
% K{1}= (4/3)*Q{2}*P{2};         
% L{1}= (1/3)*Q{2}*Q{1}*P{1};   
% T{1}=-(2/3)*hk{3}*P{3}; 
% for i=2:k-3   
%     K{i}= (4/3)*Q{i+1}*P{i+1};    
%     L{i}= (1/3)*Q{i+1}*Q{i}*P{i};   
%     T{i}=-(2/3)*hk{i+2}*P{i+2}; 
% end
% %-----------coefficient matrix size-----%
% SI=size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1)+size(Phat{2},1);
% SJ=size(P{1},2)+size(A{1},2)+2*size(P{2},2);
% for i=3:k-1
%     SI=SI+size(L{i-2},1)+size(Phat{i},1);
%     SJ=SJ+size(P{i},2)+size(T{i-2},2);
% end
% %---------------------------------------%
% %---------------------------coefficient matrix assembly-------------------%
% MAT=zeros(SI,SJ);
% MAT(1:size(P{1},1),1:size(P{1},2))=P{1};
% MAT(size(P{1},1)+1:size(P{1},1)+size(Phat{1},1),1:size(Phat{1},2))=-Phat{1};
% MAT(size(P{1},1)+1:size(P{1},1)+size(Phat{1},1),...
%     size(Phat{1},2)+1:size(Phat{1},2)+size(A{1},2))=A{1};
% 
% MAT(size(P{1},1)+size(Phat{1},1)+1:size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1),...
%     1:size(Q{1}*P{1},2))=-Q{1}*P{1};
% 
% MAT(size(P{1},1)+size(Phat{1},1)+1:size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1),...
%     2*size(Q{1}*P{1},2)+1:2*size(Q{1}*P{1},2)+size(P{2},2))=P{2};
% 
% MAT(size(P{1},1)+size(Phat{1},1)+1:size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1),...
%     2*size(Q{1}*P{1},2)+size(P{2},2)+1:2*size(Q{1}*P{1},2)+size(P{2},2)+size(P{2},2))=-hk{2}*P{2};
% 
% MAT(size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1)+1:size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1)+size(Phat{2},1),...
%     2*size(Q{1}*P{1},2)+1:2*size(Q{1}*P{1},2)+size(Phat{2},2))=-Phat{2};
% 
% MAT(size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1)+1:size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1)+size(Phat{2},1),...
%     2*size(Q{1}*P{1},2)+size(Phat{2},2)+1:2*size(Q{1}*P{1},2)+size(Phat{2},2)+size(A{2},2))=A{2};
% 
% indi=size(P{1},1)+size(Phat{1},1)+size(Q{1}*P{1},1)+size(Phat{2},1);
% indj=0;
% for i=3:k-1
%     MAT(indi+1:indi+size(L{i-2},1),indj+1:indj+size(L{i-2},2))=L{i-2};
%         
%     MAT(indi+1:indi+size(L{i-2},1),indj+2*size(L{i-2},2)+1:...
%         indj+2*size(L{i-2},2)+size(K{i-2},2))=-K{i-2};
%     
%     MAT(indi+1:indi+size(L{i-2},1),indj+2*size(L{i-2},2)+2*size(K{i-2},2)+1:...
%         indj+2*size(L{i-2},2)+2*size(K{i-2},2)+size(P{i},2))=P{i};        
%     
%     MAT(indi+1:indi+size(L{i-2},1),indj+2*size(L{i-2},2)+2*size(K{i-2},2)+size(P{i},2)+1:...
%         indj+2*size(L{i-2},2)+2*size(K{i-2},2)+size(P{i},2)+size(T{i-2},2))=T{i-2};  
%     
% 
%     
%     MAT(indi+size(L{i-2},1)+1:indi+size(L{i-2},1)+size(Phat{i},1),...
%         indj+2*size(L{i-2},2)+2*size(K{i-2},2)+1:...
%         indj+2*size(L{i-2},2)+2*size(K{i-2},2)+size(Phat{i},2))=-Phat{i};
%     
%     MAT(indi+size(L{i-2},1)+1:indi+size(L{i-2},1)+size(Phat{i},1),...
%         indj+2*size(L{i-2},2)+2*size(K{i-2},2)+size(Phat{i},2)+1:...
%         indj+2*size(L{i-2},2)+2*size(K{i-2},2)+size(Phat{i},2)+size(A{i},2))=A{i};
%     
%     indi=indi+size(L{i-2},1)+size(Phat{i},1);
%     indj=indj+2*size(P{i-2},2);
% end
% %-------------------------------------------------------------------------%
% VEC=zeros(size(MAT,1),1);
% VEC(1:n*N{1},1)=U{1};
% 
% mat=sum(abs(MAT),2);
% MAT=MAT./repmat(mat,1,size(MAT,2));
% VEC=VEC./mat;
% sol=MAT\VEC;
% ind=0;
% hold all
% for i=1:k-1
%     u{i}=real(sol(ind+1:ind+n*N{i}));
%     ind=ind+2*n*N{i};
% %     plot3(X{i},Y{i},P{i}*u{i});
%     err(i)=norm(U{i}-P{i}*u{i},inf)/norm(U{i},inf);
%     uap{i}=P{i}*u{i};
% end
% % plot(err,'-*')
% e(ii)=norm(err,inf);
% % end