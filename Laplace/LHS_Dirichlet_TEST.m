function A=LHS_Dirichlet_TEST(Npol,N1,N2,N3,N4,M,R,m,h)
Ncol=N1+N2+N3+N4;
m=conj(m');
h=conj(h');
for p=1:Npol
     L((p-1)*M+1:p*M,1)=-conj(h(p))*(R/M)*((1:M)');  
end
hj=[h(1)*ones(1,N1) h(2)*ones(1,N2) h(3)*ones(1,N3) h(4)*ones(1,N4)];
mj=[m(1)*ones(1,N1) m(2)*ones(1,N2) m(3)*ones(1,N3) m(4)*ones(1,N4)];
chj=[conj(h(1))*ones(1,N1) conj(h(2))*ones(1,N2) conj(h(3))*ones(1,N3) conj(h(4))*ones(1,N4)];
cmj=[conj(m(1))*ones(1,N1) conj(m(2))*ones(1,N2) conj(m(3))*ones(1,N3) conj(m(4))*ones(1,N4)];

L=repmat([L;L],1,Ncol);
hj=[repmat(hj,Npol*M,1);repmat(chj,Npol*M,1)];
mj=[repmat(mj,Npol*M,1);repmat(cmj,Npol*M,1)];
    
ll=[0:N1-1 0:N2-1 0:N3-1 0:N4-1];

he=abs(hj).*[exp(-1i*L(1:Npol*M,:).*mj(1:Npol*M,:));...
             exp(1i*L(Npol*M+1:end,:).*mj(Npol*M+1:end,:))];
S=legendre_fft([-1i*L(1:Npol*M,:);1i*L(Npol*M+1:end,:)].*...
    hj,repmat(ll,2*Npol*M,1));
A=he.*S;

    