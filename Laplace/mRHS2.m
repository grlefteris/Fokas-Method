function y=mRHS2(n,Np,M,R,m,h,u)
u=reshape(u,length(u)/n,n);
t=linspace(-1,1,Np);
%======================Interpolate known boundary data====================%
tt=linspace(-1,1,100); %evaluation accuracy
% p_order=20;%Npoints-1; % maximum = Ncol-1
tic
for j=1:n
    P(j,:)=polyfit(t',u(:,j),Np);
    w(j,:)=polyval(P(j,:),tt);
end
t1=toc
%=========================================================================%
tic
for p=1:n
    for j=1:n
        for r=1:M
            
            F(p,j)=-1i*h(j)*(-conj(h(p))*(R/M));       %frequency matrix
            Fc(p,j)=1i*conj(h(j))*(-conj(h(p))*(R/M)); %conjugate frequency matrix                  
            
            I(((p-1)*M)+r,j)=h(j)*(-conj(h(p))*(R/M)*r)*...
                exp(-1i*(-conj(h(p))*(R/M)*r)*m(j))*...
                simpson(tt,exp(F(p,j)*r*tt).*w(j,:),2,'1/3');
           
            Ic(((p-1)*M)+r,j)=conj(h(j))*(-conj(h(p))*(R/M)*r)*...
                exp(1i*(-conj(h(p))*(R/M)*r)*conj(m(j)))*...               
                simpson(tt,exp(Fc(p,j)*r*tt).*w(j,:),2,'1/3');
        end
    end
end
t2=toc
I=-sum(I,2);
Ic=-sum(Ic,2);
y=[I;Ic];

        