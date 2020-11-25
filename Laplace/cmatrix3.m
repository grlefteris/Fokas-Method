function A=cmatrix3(Npol,Ncol,M,R,m,h,problem)
A=zeros(2*Npol*M,Npol*Ncol);
Temp=zeros(2*Npol,Npol);

switch problem
    
    case 'dirichlet'

    for p=1:Npol
        for j=1:Npol
        
            for l=0:Ncol-1
                for r=1:M


                lamda=(-conj(h(p)))*(R/M)*r;
                Block(r,l+1)=  abs(h(j)) * exp(-1i*( lamda  )*m(j))*...
                    legtran((-1i*( lamda )*h(j)),l) ; 
                SBlock(r,l+1)= abs(conj(h(j))) * exp(1i*( lamda  )*conj(m(j)))*...
                    legtran((1i*( lamda )*conj(h(j))),l)  ; %Schwartz Conjugate
                     
                            
                end
            end
        
        
        Temp(p,j)=1;
        B=kron(Temp,Block);
        A=A+B;
        Temp(p,j)=0;
        Temp(Npol+p,j)=1;
        B=kron(Temp,SBlock);
        A=A+B;
        Temp(Npol+p,j)=0;
        
        
        end
    end

    case 'neumann'
        for p=1:Npol
            for j=1:Npol
        
                for l=0:Ncol-1
                    for r=1:M
   
                                                
                    lamda=(-conj(h(p)))*(R/M)*r;
                    Block(r,l+1)=  lamda * h(j) * exp(-1i*( lamda  )*m(j))*...
                      legtran((-1i*( lamda )*h(j)),l) ; 
                    SBlock(r,l+1)= lamda * conj(h(j)) * exp(1i*( lamda  )*conj(m(j)))*...
                      legtran((1i*( lamda )*conj(h(j))),l)  ; %Schwartz Conjugate
                      
               
                    end
                end
        
        
                    Temp(p,j)=1;
                    B=kron(Temp,Block);
                    A=A+B;
                    Temp(p,j)=0;
                    Temp(Npol+p,j)=1;
                    B=kron(Temp,SBlock);
                    A=A+B;
                    Temp(Npol+p,j)=0;
        
        
            end
        end
        
    case 'mixed'% else if strcmp(problem,'mixed')
        
        %=======Dirichlet Side===========%
        for p=1:Npol
        for j=[1 2 4]
        
            for l=0:Ncol-1
                for r=1:M
                
                lamda=(-conj(h(p)))*(R/M)*r;
                Block(r,l+1)=  abs(h(j)) * exp(-1i*( lamda  )*m(j))*...
                    legtran((-1i*( lamda )*h(j)),l) ; 
                SBlock(r,l+1)= abs(conj(h(j))) * exp(1i*( lamda  )*conj(m(j)))*...
                    legtran((1i*( lamda )*conj(h(j))),l)  ; %Schwartz Conjugate
                      
               
                end
            end
        
        
        Temp(p,j)=1;
        B=kron(Temp,Block);
        A=A+B;
        Temp(p,j)=0;
        Temp(Npol+p,j)=1;
        B=kron(Temp,SBlock);
        A=A+B;
        Temp(Npol+p,j)=0;
        
        
        end
        end
        %================================%
        
        %=========Neumann Side===========%
        for p=1:Npol
        for j=3
        
                for l=0:Ncol-1
                    for r=1:M
   
                                                
                    lamda=(-conj(h(p)))*(R/M)*r;
                    Block(r,l+1)=  lamda * h(j) * exp(-1i*( lamda  )*m(j))*...
                      legtran((-1i*( lamda )*h(j)),l) ; 
                    SBlock(r,l+1)= lamda * conj(h(j)) * exp(1i*( lamda  )*conj(m(j)))*...
                      legtran((1i*( lamda )*conj(h(j))),l)  ; %Schwartz Conjugate
                      
               
                    end
                end
        
        
        Temp(p,j)=1;
        B=kron(Temp,Block);
        A=A+B;
        Temp(p,j)=0;
        Temp(Npol+p,j)=1;
        B=kron(Temp,SBlock);
        A=A+B;
        Temp(Npol+p,j)=0;
        
        
        end
        end
        %================================%
        
    
    otherwise 
        disp('Error: Choose Dirichlet or Neumann Boundary Conditions');
end

        

         