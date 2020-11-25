function A=PLHS_Dirichlet(Npol,Ncol,M,R,m,h,scale)
for p=1:Npol
     L((p-1)*M+1:p*M,1)=-conj(h(p)/scale)*(R/M)*((1:M)');
     hj(1,(p-1)*Ncol+1:p*Ncol)=h(p);
    chj(1,(p-1)*Ncol+1:p*Ncol)=conj(h(p));  
     mj(1,(p-1)*Ncol+1:p*Ncol)=m(p);
    cmj(1,(p-1)*Ncol+1:p*Ncol)=conj(m(p));    
end
L=repmat([L;L],1,Npol*Ncol);
hj=[repmat(hj,Npol*M,1);repmat(chj,Npol*M,1)];
mj=[repmat(mj,Npol*M,1);repmat(cmj,Npol*M,1)];
    
he=abs(hj/scale).*[exp(-1i*L(1:Npol*M,:).*(mj(1:Npol*M,:)/scale));...
             exp(1i*L(Npol*M+1:end,:).*(mj(Npol*M+1:end,:)/scale))];
S=legendre_fft([-1i*L(1:Npol*M,:);1i*L(Npol*M+1:end,:)].*...
    (hj/scale),repmat(0:Ncol-1,2*Npol*M,Npol));
A=he.*S;

    