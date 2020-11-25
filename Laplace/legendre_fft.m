function S=legendre_fft(z,l)
% S=zeros(1,length(z));
% S(find((z==0)&(l==0)))=2;
% S(find(z~=0))=(sqrt(2*pi*z(find(z~=0)))./z(find(z~=0))).*...
%     besseli(l(find(z~=0))+0.5,z(find(z~=0)));
[n m]=size(z);
S=zeros(n,m);
S(find((z==0)&(l==0)))=2;
S(find(z~=0))=(sqrt(2*pi*z(find(z~=0)))./z(find(z~=0))).*...
    besseli(l(find(z~=0))+0.5,z(find(z~=0)));