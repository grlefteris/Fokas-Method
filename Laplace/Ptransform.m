function B=Ptransform(n,N,M,R,m,h)


lamda=kron((-conj(h)),(R/M)*(1:M)');
l=repmat(0:N-1,n*M,n);
h=conj(h');%conjugate transpose attention!
m=conj(m');%conjugate transpose attention!
w= repmat(-1i*lamda,1,n*N).*repmat(kron(h,      ones(1,N)),n*M,1);
ws=repmat( 1i*lamda,1,n*N).*repmat(kron(conj(h),ones(1,N)),n*M,1);
Ph= legendre_fft(w, l);
Phs=legendre_fft(ws,l);
E= repmat(kron(h,      ones(1,N)),n*M,1).*repmat(lamda,1,n*N).*...
  exp(-1i.*repmat(lamda,1,n*N).*repmat(kron(m,      ones(1,N)),n*M,1));
Es=repmat(kron(conj(h),ones(1,N)),n*M,1).*repmat(lamda,1,n*N).*...
  exp( 1i.*repmat(lamda,1,n*N).*repmat(kron(conj(m),ones(1,N)),n*M,1));
B=-[Ph.*E;Phs.*Es];

    