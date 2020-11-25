clear;clc;
inner_iterations=0;
Npol=4;
Ncol=11;
k=(Ncol+1)/2;
M1=2;%2*41;
R1=5;%8*M1;


rad_z=1;
rot_z=0;
uan=@(x,y) exp(1+x).*cos(2+y);
%===========================LEVEL 1========================================
[A q dom z ang]=G_D2N_Laplace(Npol,Ncol,M1,R1,rad_z,rot_z,'dirichlet');
axis_x=real(z(1:end-1));
axis_y=imag(z(1:end-1));
qq=reshape(q,Npol*Ncol,1);
%==============Decompose Normal derivative===============%
for j=1:Npol
    Qx(:,j)=-abs(cos(ang(j)-(pi/2)))*q(:,j);
    Qy(:,j)=-abs(sin(ang(j)-(pi/2)))*q(:,j);
end
QQx=kron(-abs(cos(ang-(pi/2))),ones(Ncol,1)).*qq;
QQy=kron(-abs(sin(ang-(pi/2))),ones(Ncol,1)).*qq;
QQx=reshape(QQx,Ncol,Npol);
QQy=reshape(QQy,Ncol,Npol);
%========================================================%
%==========================================================================
for j=1:Npol-1
    x(:,j)=linspace(axis_x(j),axis_x(j+1),Ncol);
    y(:,j)=linspace(axis_y(j),axis_y(j+1),Ncol);
end
x(:,Npol)=linspace(axis_x(Npol),axis_x(1),Ncol);
y(:,Npol)=linspace(axis_y(Npol),axis_y(1),Ncol);
X=reshape(x,Npol*Ncol,1);
Y=reshape(y,Npol*Ncol,1);
u=reshape(uan(x,y),Ncol,Npol);
U=reshape(u,Npol*Ncol,1);
hk=sqrt((x(1,1)-x(2,1))^2+(y(1,1)-y(2,1))^2)/tan(ang(1)-(pi/2));
w=ang-pi/2;
hx=abs(cos(w)*hk);
hy=abs(sin(w)*hk);
hx=repmat(hx',Ncol,1);
hy=repmat(hy',Ncol,1);
err(1)=norm(reshape(u,Npol*Ncol,1)-reshape(uan(x,y),Npol*Ncol,1),2);
% plot3(x,y,u,'b*');
% hold all
%=============================NEXT LEVELS==================================
j=1;
for i=2:1%k-1
    v=u;
    Sx=-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));
    Sy=-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));
    Sx=Sx(2:end-1,:);
    Sy=Sy(2:end-1,:);    
    hx=hx(2:end-1,:);
    hy=hy(2:end-1,:);    
    x=x(2:end-1,:);
    y=y(2:end-1,:);
    hx=hx.*Sx;
    hy=hy.*Sy;    
    x=x+hx;
    y=y+hy;
    X=[X;reshape(x,size(x,1)*Npol,1)];
    Y=[Y;reshape(y,size(y,1)*Npol,1)];
    xx=x(1,:);
    yy=y(1,:);
    hx=abs(hx);
    hy=abs(hy);
    u=u(2:end-1,:);
    Qx=Qx(2:end-1,:);
    Qy=Qy(2:end-1,:);
    u(2:end-1,:)=u(2:end-1,:)+...
        (hx(2:end-1,:).*Qx(2:end-1,:))+(hy(2:end-1,:).*Qy(2:end-1,:));
    for l=2:Npol
        u(1,l)=((v(end-1,l-1)+...
            hx(end-1,l-1).*Qx(end-1,l-1)+hy(end-1,l-1).*Qy(end-1,l-1))+...
            (v(2,l)+...
            hx(2,l).*Qx(2,l)+hy(2,l).*Qy(2,l)))/2;
        u(end,l-1)=u(1,l);   
    end
    u(1,1)=((v(end-1,Npol)+...
        hx(end-1,Npol).*Qx(end-1,Npol)+hy(end-1,Npol).*Qy(end-1,Npol))+...
        (v(2,1)+...
        hx(2,1).*Qx(2,1)+hy(2,1).*Qy(2,1)))/2;
    u(end,Npol)=u(1,1);      
               
    QXold=Qx;
    QYold=Qy;
    
    xmax=xx(1);
    scale=1;
    MMM=Npol*(Ncol-2*(i-1));
    RRR=MMM/512;
    [Q dom z ang E(i)]=PDtN_Laplace(Npol,Ncol-2*(i-1),...
                  2,5,rad_z,rot_z,xx*scale,yy*scale,'dirichlet',scale,u);
    for t=1:Npol
    Qx(:,t)=-abs(cos(w(t)))*Q(:,t);
    Qy(:,t)=-abs(sin(w(t)))*Q(:,t);
    end
    
%     %===Implicit Step===%
%     V=v;
%     v=v(2:end-1,:);   
%     for iter=1:inner_iterations
%     u=v+((hx/2).*(Qx+QXold))+((hy/2).*(Qy+QYold));
%     for l=2:Npol
%         u(1,l)=((V(end-1,l-1)+...
%             (hx(end-1,l-1)/2).*(Qx(end-1,l-1)+QXold(end-1,l-1))+...
%             (hy(end-1,l-1)/2).*(Qy(end-1,l-1)+QYold(end-1,l-1)))+...
%             (V(2,l)+...
%             (hx(2,l)/2).*(Qx(2,l)+QXold(2,l))+...
%             (hy(2,l)/2).*(Qy(2,l)+QYold(2,l))))/2;
%         u(end,l-1)=u(1,l);   
%     end
%     u(1,1)=((V(end-1,Npol)+...
%         (hx(end-1,Npol)/2).*(Qx(end-1,Npol)+QXold(end-1,Npol))+...
%         (hy(end-1,Npol)/2).*(Qy(end-1,Npol)+QYold(end-1,Npol)))+...
%         (V(2,1)+...
%         (hx(2,1)/2).*(Qx(2,1)+QXold(2,1))+...
%         (hy(2,1)/2).*(Qy(2,1)+QYold(2,1))))/2;
%     u(end,Npol)=u(1,1);  
%     [Q dom z ang]=DtN_Laplace(Npol,Ncol-2*(i-1),...
%                   M,R,rad_z,rot_z,xx,yy,'dirichlet',u);
%     for t=1:Npol
%     Qx(:,t)=-abs(cos(w(t)))*Q(:,t);
%     Qy(:,t)=-abs(sin(w(t)))*Q(:,t);
%     end    
%     end
%     %==================%
    
    U=[U;reshape(u,size(u,1)*Npol,1)];
    err(i)=(norm(reshape(u,Npol*(Ncol-2*(i-1)),1)-...
           reshape(uan(x,y),Npol*(Ncol-2*(i-1)),1),inf))/...
           norm(reshape(uan(x,y),Npol*(Ncol-2*(i-1)),1),inf);
%    MM(j)=getframe;
%    plot3(x,y,u,'*');
    j=j+1;    
end
% Sx=-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));
% Sy=-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));
% Sx=Sx(2:end-1,:);Sy=Sy(2:end-1,:);hx=hx(2:end-1,:);hy=hy(2:end-1,:);    
% x=x(2:end-1,:);y=y(2:end-1,:);hx=hx.*Sx;hy=hy.*Sy;x=x+hx;y=y+hy;
% xx=x(1,:);yy=y(1,:);hx=abs(hx);hy=abs(hy);u=u(2:end-1,:);
% Qx=Qx(2:end-1,:);Qy=Qy(2:end-1,:);
% u(2:end-1,:)=u(2:end-1,:)+...
%         (hx(2:end-1,:).*Qx(2:end-1,:))+(hy(2:end-1,:).*Qy(2:end-1,:));
% u=sum(u)/Npol;   
% err(k)=norm(u-uan(0,0),inf)/norm(uan(0,0),inf); 
% % plot3(0,0,u,'*');  
% norm(err,inf)

% Pol=polyFit2D(U(size(X)-Npol*(15+13+11)+1:end),...
%               X(size(X)-Npol*(15+13+11)+1:end),...
%               Y(size(Y)-Npol*(15+13+11)+1:end),4,4);

%  Pol=polyFit2D(U,X,Y,2,2);
% for ii=i+1:k
%     Sx=-((x>0)&(y>0))+((x<0)&(y>0))+((x<0)&(y<0))-((x>0)&(y<0));
%     Sy=-((x>0)&(y>0))-((x<0)&(y>0))+((x<0)&(y<0))+((x>0)&(y<0));
%     Sx=Sx(2:end-1,:);
%     Sy=Sy(2:end-1,:);    
%     hx=hx(2:end-1,:);
%     hy=hy(2:end-1,:);    
%     x=x(2:end-1,:);
%     y=y(2:end-1,:);
%     hx=hx.*Sx;
%     hy=hy.*Sy;    
%     x=x+hx;
%     y=y+hy;
%     co(ii-i)=length(X);
%     X=[X;reshape(x,size(x,1)*Npol,1)];
%     Y=[Y;reshape(y,size(y,1)*Npol,1)];
%     cn(ii-i)=length(X);
%     hx=abs(hx);
%     hy=abs(hy);
% end
% un=polyVal2D(Pol,X,Y,2,2);
% for ii=i+1:k
%     err(ii)=norm(un(co(ii-i)+1:cn(ii-i))-...
%         uan(X(co(ii-i)+1:cn(ii-i)),Y(co(ii-i)+1:cn(ii-i))),inf)/...
%         norm(uan(X(co(ii-i)+1:cn(ii-i)),Y(co(ii-i)+1:cn(ii-i))),inf);
% end
