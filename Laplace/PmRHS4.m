function G=PmRHS4(Npol,Ncol,M,R,m,h,scale,u)
t=linspace(-1,1,Ncol);
%======================Interpolate known boundary data====================%
tt=linspace(-1,1,1000); %evaluation accuracy
p_order=3;%Ncol-1; % maximum = Ncol-1

for j=1:Npol
    P(j,:)=polyfit(t',u(:,j),p_order);
    w(j,:)=polyval(P(j,:),tt);
end

%=========================================================================%

for p=1:Npol
    for j=1:Npol
        for r=1:M
            
            F(p,j)=-1i*(h(j)/scale)*(-conj(h(p)/scale)*(R/M));       %frequency matrix
            Fc(p,j)=1i*conj((h(j)/scale))*(-conj(h(p)/scale)*(R/M)); %conjugate frequency matrix                  
            
            I(((p-1)*M)+r,j)=(h(j)/scale)*(-conj(h(p)/scale)*(R/M)*r)*...
                exp(-1i*(-conj(h(p)/scale)*(R/M)*r)*(m(j)/scale))*...
                simpson(tt,exp(F(p,j)*r*tt).*w(j,:),2,'1/3');
           
            Ic(((p-1)*M)+r,j)=conj((h(j)/scale))*(-conj(h(p)/scale)*(R/M)*r)*...
                exp(1i*(-conj(h(p)/scale)*(R/M)*r)*conj((m(j)/scale)))*...               
                simpson(tt,exp(Fc(p,j)*r*tt).*w(j,:),2,'1/3');
        end
    end
end
I=-sum(I,2);
Ic=-sum(Ic,2);
G=[I;Ic];