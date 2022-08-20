function EQ = def2(x,lambda,kappa0,d0,ds,eps_s,Mf,wf,mu1,mu2,flag)

latx = x(1);
laty = x(2);


F1 = [latx 0 0; 0 laty 0; 0 0 lambda];
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
S1 = Pfiber(kappa0,d0,ds,eps_s,F1,F1inv,C1,J1,Mf,mf,eps_norm,wf,mu1,mu2,flag);


EQ(1) = S1(1,1);
EQ(2) = S1(2,2);


end