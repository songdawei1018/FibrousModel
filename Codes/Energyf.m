function W = Energyf(kappa0,d0,ds,eps_s,F,CF,JF,eps_norm,wf,mu1,mu2,flag)

%kappa0: fiber density x fiber stiffness
%d0:strain scale for softening
%ds:strain scale for stiffening
%eps_s:critical strain for stiffening
%F:deformation gradient 3x3
%Mf:discretizations for fibers
%wf:associated weights for fibers
%mu1 and mu2: Lame constants of the isotropic components

kk = length(wf); 
W = 0;
lam_s = 1 + eps_s;

for p = 1:kk
   
    eps = eps_norm(p);
    
    % flag = 1
    if abs(flag - 1) < 1e-2
    
      
    if eps < 0                
        x = d0*(d0*exp(eps/d0) - eps - d0);   
    elseif eps >= eps_s        
        x = 1/2*(eps_s)^2 + (eps_s - ds)*(eps - eps_s) + ds^2*(exp((eps-eps_s)/ds)-1);        
    else        
        x = 1/2*eps^2;             
    end
    
    
    end
    
    
    % flag = 2
    if abs(flag - 2) < 1e-2
        
        lamn2 = ((1 + eps)/(lam_s))^2;        
        if lamn2 >= 1           
            x = 1/(2*ds)*(exp(ds*(lamn2 - 1)^2) - 1);           
        else           
            x = 0;            
        end
         
    end
    

    
    if abs(flag - 3) < 1e-2
        
        lam = 1 + eps;        
        if lam >= lam_s           
            x = 1/2*(lam - lam_s)^2;            
        else           
            x = 0;            
        end
         
    end    
    
    W = W + x*wf(p);
       
end

W = W*kappa0/(2*pi);

I1 = trace(CF);
W = W + mu1/2*(I1 - 3 - 2*log(JF)) + mu2/2*(JF - 1)^2;


end