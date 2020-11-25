function B=Ptransform_TEST(n,N1,N2,N3,N4,M,R,m,h)

N=N1+N2+N3+N4;
lamda=kron((-conj(h)),(R/M)*(1:M)');
l=repmat([0:N1-1 0:N2-1 0:N3-1 0:N4-1],n*M,1);
h=conj(h');%conjugate transpose attention!
m=conj(m');%conjugate transpose attention!
w= repmat(-1i*lamda,1,N).*repmat([h(1)*ones(1,N1) h(2)*ones(1,N2) h(3)*ones(1,N3) h(4)*ones(1,N4)],n*M,1);
ws=repmat( 1i*lamda,1,N).*repmat([conj(h(1))*ones(1,N1) conj(h(2))*ones(1,N2) conj(h(3))*ones(1,N3) conj(h(4))*ones(1,N4)],n*M,1);
Ph= legendre_fft(w, l);
Phs=legendre_fft(ws,l);

E= repmat([h(1)*ones(1,N1) h(2)*ones(1,N2) h(3)*ones(1,N3) h(4)*ones(1,N4)],n*M,1).*repmat(lamda,1,N).*...
  exp(-1i.*repmat(lamda,1,N).*repmat([m(1)*ones(1,N1) m(2)*ones(1,N2) m(3)*ones(1,N3) m(4)*ones(1,N4)],n*M,1));
Es=repmat([conj(h(1))*ones(1,N1) conj(h(2))*ones(1,N2) conj(h(3))*ones(1,N3) conj(h(4))*ones(1,N4)],n*M,1).*repmat(lamda,1,N).*...
  exp( 1i.*repmat(lamda,1,N).*repmat([conj(m(1))*ones(1,N1) conj(m(2))*ones(1,N2) conj(m(3))*ones(1,N3) conj(m(4))*ones(1,N4)],n*M,1));
B=-[Ph.*E;Phs.*Es];

    