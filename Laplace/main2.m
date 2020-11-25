% function [u,err,hs]=main(n,Nl,Np)

% This function returns the approximated solution u for each spatial-level
% k, in the interior of the computational domain, and the corresponding
% numerical error err. In all cases, Dirichlet boundary conditions are used. 
% The computational domain is a regular polygon inscribed in the unit 
% circle with n sides. Nl is defined as the order of the Legendre basis 
% functions, whereas Np is defined as the number of grid points per side of
% the polygon on the boundaries. The corresponding numerical errors are
% computed according to a given analytical solution U.
clear;clc
warning off; 

%------ANALYTICAL SOLUTION------%  % An analytical solution for benchmarking
uan = @(x,y)  exp(1+x).*cos(2+y);
ux  = @(x,y)  exp(1+x).*cos(2+y);
uy  = @(x,y) -exp(1+x).*sin(2+y);
%------ANALYTICAL SOLUTION------%
n=4;
Nl=11;
Np=21;
M=n*Nl; % Predefined collocation parameters according to Fokas heuristic 
R=2*M;



[X,Y,k,hs,xx,yy,w,N]=domain(n,Np,1,0); % A function returning the spatial coordinates X, Y, for each spatial-level k
                                       % the spatial-stepping hs, the coordinates of the corners of the polygon xx, yy, 
                                       % the angles w between the normal vector and the x-axis for each side of the
                                       % polygon, and the number of grid points for each spatial-level k.
                                       
                                       
%-------- Compute unknown Neumann values on the boundary------------------%

U{1}=uan(X{1},Y{1}); % Given Dirichlet boundary conditions
u{1}=U{1};           % The solution on the boundary (level 1) is known, 
                     % thus only the Neumann values are unknown
                     
W{1}=kron(w,ones(N{1},1));
q{1}=(cos(W{1}).*ux(X{1},Y{1}))+(sin(W{1}).*uy(X{1},Y{1}));
P{1}=kron(eye(n),Lbasis(Nl,N{1}));                   % Legendre basis matrix
[m{1},h{1}]=sides(n,0,xx{1},yy{1});
A{1} = LHS(n,Nl,M,R,m{1},h{1});                      % Coefficient matrix A
s{1} = RHS(n,N{1},M,R,m{1},h{1},u{1});               % Right hand side s
scale{1}=sum(abs(A{1}),2);                           % scaling factor
A{1}=A{1}./repmat(scale{1},1,size(A{1},2));
s{1}=s{1}./scale{1};
x{1}=real(A{1}\s{1});    % Solving for the Legendre expansion coefficients 
                         % of the unknown Neumann values
norm(q{1}-P{1}*x{1},inf)/norm(q{1},inf)                         
%-------- Compute unknown Neumann values on the boundary------------------%

