function G=mRHS3(Npol,Ncol,M,R,m,h,u)
t=linspace(-1,1,Ncol);
%======================Interpolate known boundary data====================%
tt=linspace(-1,1,10000); %evaluation accuracy
p_order=Ncol-1; % maximum = Ncol-1
p1=polyfit(t',u(:,1),p_order);
p2=polyfit(t',u(:,2),p_order);
p3=polyfit(t',u(:,3),p_order);
p4=polyfit(t',u(:,4),p_order);

w1=polyval(p1,tt);
w2=polyval(p2,tt);
w3=polyval(p3,tt);
w4=polyval(p4,tt);

% w1=spline(t,u(:,1),tt);
% w2=spline(t,u(:,2),tt);
% w3=spline(t,u(:,3),tt);
% w4=spline(t,u(:,4),tt);
w=[w1;w2;w3;w4];
%=========================================================================%

for p=1:Npol
    for j=1:Npol
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
I=-sum(I,2);
Ic=-sum(Ic,2);
G=[I;Ic];
            
            
            
            
       


























