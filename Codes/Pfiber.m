function Pf = Pfiber(kappa0,d0,ds,eps_s,F,Finv,CF,JF,Mf,mf,eps_norm,wf,lame,alpha)

%kappa0: fiber density x fiber stiffness
%d0:strain scale for softening
%ds:strain scale for stiffening
%eps_s:critical strain for stiffening
%F:deformation gradient 3x3
%Mf:discretizations for fibers
%wf:associated weights for fibers

kk = length(wf);
Pf = zeros(3,3);
I1 = trace(CF);
lam_s = eps_s + 1;


for i = 1:3
for j = 1:3
                
    
for p = 1:kk  %sum over all fiber orientations
   
    M = Mf(:,p); % 3 x 1 vector
    m = mf(:,p);   % The p-th fiber orientation in the current config.
    
    eps = eps_norm(p);
    lambda = eps + 1;
           
    if eps < 0 

        wpp = exp(eps/d0);
        wp = d0*(wpp - 1);            
          
    elseif eps >= eps_s
        
        wpp = exp((eps-eps_s)/ds);
        wp = eps_s + ds*(wpp - 1);
        
    else
        
        wp = eps;
               
    end
    
        Pf(i,j) = Pf(i,j) + wp/lambda*m(i)*M(j)*wf(p);
              
end
    
              
end 
end


%Pf = Pf*kappa0/(2*pi) + lame*log(JF)*exp(alpha*(log(JF))^2)*transpose(Finv);

Pf = Pf*kappa0/(2*pi) + (lame/alpha)*sinh(alpha*(JF - 1))*JF*transpose(Finv);


end