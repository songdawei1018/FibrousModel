function W = Energyinc(F,CF,JF,mu1,Im)

%kappa0: fiber density x fiber stiffness
%d0:strain scale for softening
%ds:strain scale for stiffening
%eps_s:critical strain for stiffening
%F:deformation gradient 3x3
%Mf:discretizations for fibers
%wf:associated weights for fibers
%mu1 and mu2: Lame constants of the isotropic components


I1 = trace(CF);
W = -Im*mu1/2*log(1-(I1-3)/Im);

end