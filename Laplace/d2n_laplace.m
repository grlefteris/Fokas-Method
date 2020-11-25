function [unapp z_side z_corner ang]=d2n_laplace(Npol,Ncol,M,R,rad_z,rot_z,axis_x,axis_y,problem,u)
%[z_side z_corner m h ang]=polygon(Npol,Ncol,rad_z,rot_z);
[z_side z_corner ang  m h]=irpolygon2(Npol,Ncol,rot_z,axis_x,axis_y);
A=cmatrix3(Npol,Ncol,M,R,m,h,problem);
Am = sum(abs(A),2);
A = A./repmat(Am,1,Npol*Ncol);
if nargin<10
[G lamda]=rhs3(Npol,M,R,m,h,ang,problem);
else 
    G=mRHS4(Npol,Ncol,M,R,m,h,u);
end
G=G./Am;
c=A\G; 
c=real(c);
c=reshape(c,Ncol,Npol);
phi=leg_basis(Ncol);
unapp=zeros(Ncol,Npol);
for j=1:Npol
    for i=1:Ncol
        
        summ=0;
        for l=1:Ncol
            summ=summ + (c(l,j)*phi(i,l));
        end
        
        unapp(i,j)=summ;
        
    end
end     

%===========Display D2N Errors on every polygon's side=====================
if strcmp(problem,'dirichlet') 
for j=1:Npol
    x=linspace(real(z_corner(j)),real(z_corner(j+1)),Ncol);
    y=linspace(imag(z_corner(j)),imag(z_corner(j+1)),Ncol);
    angx=ang(j);
    angy=angx+(pi/2);
    q=@(x,y) cos(angx)*(exp(1+x).*cos(2+y)) +...
        cos(angy)*(-exp(1+x).*sin(2+y));
    err=max(abs(unapp(:,j)-q(x,y)'));
    fprintf('Error of side %i = %e\n',j,err);
    err=norm(unapp(:,j)-(q(x,y))',inf)/norm(q(x,y),inf);
    fprintf('NormError of side %i = %e\n',j,err);
end
else if strcmp(problem,'neumann')
for j=1:Npol
    x=linspace(real(z_corner(j)),real(z_corner(j+1)),Ncol);
    y=linspace(imag(z_corner(j)),imag(z_corner(j+1)),Ncol);
    
    q=@(x,y) exp(1+x).*cos(2+y);
    err=max(abs(unapp(:,j)-q(x,y)'));
    fprintf('Error of side %i = %e\n',j,err);
    err=norm(unapp(:,j)-q(x,y)',inf)/norm(q(x,y),inf);
    fprintf('NormError of side %i = %e\n',j,err);
end
else
        disp('Error: Choose Dirichlet or Neumann Boundary Conditions');
    end
end
disp('\n');
%==========================================================================