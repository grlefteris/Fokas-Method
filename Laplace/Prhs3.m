function [G colpoints]=Prhs3(Npol,M,R,m,h,ang,problem,scale)
G=zeros(2*Npol*M,1);

switch problem

    case 'dirichlet'
        col=1; 
    for p=1:Npol
    for r=1:M
        
        sum1=0;
        sum2=0;
        lamda=(-conj(h(p)/scale))*(R/M)*r;
        colpoints(col)=lamda;
        col=col+1;
        
        for i=1:Npol
            
            integrand1= @(t) exp(-1i*lamda*(h(i)/scale)*t).*...
                exp(1+real((m(i)+(t*h(i)))/scale)).*cos(2+imag((m(i)+(t*h(i)))/scale));
                
            sum1=sum1 + ( lamda*(h(i)/scale)*exp(-1i*lamda*(m(i)/scale))*...
                quadgk(integrand1,-1,1) );
            
            integrand2= @(t) exp(1i*lamda*conj((h(i)/scale))*t).*...
                exp(1+real((m(i)+(t*h(i)))/scale)).*cos(2+imag((m(i)+(t*h(i)))/scale));
                
            sum2=sum2 + ( lamda*conj((h(i)/scale))*exp(1i*lamda*conj((m(i)/scale)))*...
                quadgk(integrand2,-1,1) );
            
        end
            
            
        G(((p-1)*M)+r)=-sum1;    
        G(((p-1)*M)+r+(Npol*M))=-sum2;        
    end
    end

    case 'neumann'
        col=1; 
    for p=1:Npol
    for r=1:M
        
        sum1=0;
        sum2=0;
        lamda=(-conj(h(p)))*(R/M)*r;
        colpoints(col)=lamda;
        col=col+1;
        
        for i=1:Npol
            
            angx=ang(i);
            angy=angx+(pi/2);
            
            integrand1= @(t) exp(-1i*lamda*h(i)*t).*...
                (cos(angx)*(exp(1+real(m(i)+(t*h(i)))).*cos(2+imag(m(i)+(t*h(i)))))+...
                cos(angy)*(-exp(1+real(m(i)+(t*h(i)))).*sin(2+imag(m(i)+(t*h(i))))));
            
                
            sum1=sum1 + ( abs(h(i))*exp(-1i*lamda*m(i))*...
                quadgk(integrand1,-1,1) );

            
            integrand2= @(t) exp(1i*lamda*conj(h(i))*t).*...
                (cos(angx)*(exp(1+real(m(i)+(t*h(i)))).*cos(2+imag(m(i)+(t*h(i)))))+...
                cos(angy)*(-exp(1+real(m(i)+(t*h(i)))).*sin(2+imag(m(i)+(t*h(i))))));
                
            sum2=sum2 + ( abs(conj(h(i)))*exp(1i*lamda*conj(m(i)))*...
                quadgk(integrand2,-1,1) );
            
            
        end
            
            
        G(((p-1)*M)+r)=-sum1;    
        G(((p-1)*M)+r+(Npol*M))=-sum2;        
    end
    end

    case 'mixed'
       
       %=======Dirichlet Side==========% 
       col=1; 
       for p=[1 2 4]
       for r=1:M
        
        sum1=0;
        sum2=0;
        lamda=(-conj(h(p)))*(R/M)*r;
        colpoints(col)=lamda;
        col=col+1;
        
        for i=[1 2 4]%1:Npol
            
            integrand1= @(t) exp(-1i*lamda*h(i)*t).*...
                exp(1+real(m(i)+(t*h(i)))).*cos(2+imag(m(i)+(t*h(i))));
                
            sum1=sum1 + ( lamda*h(i)*exp(-1i*lamda*m(i))*...
                quadgk(integrand1,-1,1) );
            
            integrand2= @(t) exp(1i*lamda*conj(h(i))*t).*...
                exp(1+real(m(i)+(t*h(i)))).*cos(2+imag(m(i)+(t*h(i))));
                
            sum2=sum2 + ( lamda*conj(h(i))*exp(1i*lamda*conj(m(i)))*...
                quadgk(integrand2,-1,1) );
            
        end
            
            
        G(((p-1)*M)+r)=-sum1;    
        G(((p-1)*M)+r+(Npol*M))=-sum2;        
       end
       end
       %===============================% 
        
       %=========Neumann Side==========%
       %col=1; 
       for p=3
       for r=1:M
        
        sum1=0;
        sum2=0;
        lamda=(-conj(h(p)))*(R/M)*r;
        colpoints(col)=lamda;
        col=col+1;
        
        for i=3%1:Npol
            
            angx=ang(i);
            angy=angx+(pi/2);
            
            integrand1= @(t) exp(-1i*lamda*h(i)*t).*...
                (cos(angx)*(exp(1+real(m(i)+(t*h(i)))).*cos(2+imag(m(i)+(t*h(i)))))+...
                cos(angy)*(-exp(1+real(m(i)+(t*h(i)))).*sin(2+imag(m(i)+(t*h(i))))));
            
                
            sum1=sum1 + ( abs(h(i))*exp(-1i*lamda*m(i))*...
                quadgk(integrand1,-1,1) );

            
            integrand2= @(t) exp(1i*lamda*conj(h(i))*t).*...
                (cos(angx)*(exp(1+real(m(i)+(t*h(i)))).*cos(2+imag(m(i)+(t*h(i)))))+...
                cos(angy)*(-exp(1+real(m(i)+(t*h(i)))).*sin(2+imag(m(i)+(t*h(i))))));
                
            sum2=sum2 + ( abs(conj(h(i)))*exp(1i*lamda*conj(m(i)))*...
                quadgk(integrand2,-1,1) );

      
            
        end
            
            
        G(((p-1)*M)+r)=-sum1;    
        G(((p-1)*M)+r+(Npol*M))=-sum2;        
       end
       end
       %===============================% 


    otherwise
        disp('Error: Choose Dirichlet or Neumann Boundary Conditions');
end
