function sse = sseval(x,K,Mf,wf,xdata1,xdata2,ydata1,ydata2)  %Compute the data difference

d0 = x(1); 
ds = x(2); 
lame =x(3);
alpha = x(4); 
eps = x(5);
%K = x(5);



ndat1 = length(xdata1);
ndat2 = length(xdata2);
norm1 = norm(ydata1);
norm2 = norm(ydata2);
factor = norm1^2/norm2^2;


kk = length(Mf);
eps_norm = zeros(kk,1);
mf = zeros(3,kk);
%stretchx1 = zeros(ndat1,1);
%stretchy1 = zeros(ndat1,1);


options = optimset('MaxFunEvals',10000);
sse = 0;

%y0 = [1,1];

% Compute the data match for the 1st data set.

for i=1:ndat1 % Loop over all the data points. 
    
    sz = xdata1(i);  %axial strain
        
    %fun = @(y)def2(y,sz,K,d0,ds,eps,Mf,wf,mu1,mu2,1);
    
    %if i > 1      
    %    y0 = yopt;     
    %end
    
    %yopt = fsolve(fun,y0,options); 
    %sx = yopt(1);
    %sy = yopt(2);

    sx = 1;
    sy = 1;

    %stretchx(i,1) = sx;
    %stretchy(i,1) = sy;
    
    F1 = [sx 0 0; 0 sy 0; 0 0 sz];
    F1inv = inv(F1);
    F1T = transpose(F1);
    C1 = F1T*F1;
    J1 = det(F1);
    
    %Compute the deformed fiber orientation
    for j = 1:length(Mf)      
        mf(:,j) = F1*Mf(:,j);
        eps_norm(j,1) = norm(mf(:,j)) - 1;
    end
    

    %Compute S1 for the matrix phase
    S1 = Pfiber(K,d0,ds,eps,F1,F1inv,C1,J1,Mf,mf,eps_norm,wf,lame,alpha);
      
    sse = sse + (ydata1(i)-S1(3,3))^2; % Compute the difference in data. 
    
end


% Compute the data match for the 2nd data set.

%sse = 0;

gamm = 0.01;
%y0 = [1,1];

for i=1:ndat2 % Loop over all the data points. 
    
    F33 = xdata2(i);  %axial strain
    %F11 = stretchx(i);
    %F22 = stretchy(i);
    
    %fun = @(y)def2(y,F33,K,d0,ds,eps,Mf,wf,mu1,mu2,1);
    
    %if i > 1      
    %    y0 = yopt;     
    %end
    
    %yopt = fsolve(fun,y0,options); 
    %F11 = yopt(1);
    %F22 = yopt(2);
    
    F11 = 1;
    F22 = 1;
          
    F1 = [F11 0 gamm*F33; 0 F22 0; 0 0 F33];
    F1inv = inv(F1);
    F1T = transpose(F1);
    C1 = F1T*F1;
    J1 = det(F1);

    %Compute the deformed fiber orientation
    for j = 1:length(Mf)      
        mf(:,j) = F1*Mf(:,j);
        eps_norm(j,1) = norm(mf(:,j)) - 1;
    end
    
    %Compute S1 for the matrix phase
    S1 = Pfiber(K,d0,ds,eps,F1,F1inv,C1,J1,Mf,mf,eps_norm,wf,lame,alpha); 
    
   
    modulus = S1(1,3)/gamm;
    sse = sse + factor*(ydata2(i)-modulus)^2; % Compute the difference in data. 
    
end


end


