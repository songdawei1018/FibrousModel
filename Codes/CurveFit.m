clear all;
clc;

% Read fiber orientations and weights
load('./M_404.mat');
load('./weight_404.mat');

Mf = transpose(n);
wf = w2; 
%Mf = transpose(direction);
%wf = area_sph;
kk = length(Mf);


data1=load('./Data/sample1.txt');
data2=load('./Data/sample2.txt');
data3=load('./Data/sample3.txt');
data123 = [data1 data2 data3];
data_tension = mean(data123,2);


data4=load('./Data/sample4.txt');
data5=load('./Data/sample5.txt');
data6=load('./Data/sample6.txt');
data456 = [data4 data5 data6];
data_compression = mean(data456,2);


% Compute the axial stress and axial strain.
stress = [data_tension;data_compression];
strain = 1/100*[0;1;2;3;4;5;6;7;8;9;10;-1;-2;-3;-4;-5;-6;-7;-8;-9;-10];
stretch = strain + 1;

%--------------------------------------------------

shear1=load('./Data/shear1.txt');
shear2=load('./Data/shear2.txt');
shear3=load('./Data/shear3.txt');
shear123 = [shear1 shear2 shear3];
shear_tension = mean(shear123,2);
shear_tension(1)=48.2966666667;


shear4=load('./Data/shear4.txt');
shear5=load('./Data/shear5.txt');
shear6=load('./Data/shear6.txt');
shear456 = [shear4 shear5 shear6];
shear_compression = mean(shear456,2);


storage = [shear_tension;shear_compression];


data_srp1 = dlmread('./Data/collagen_axial.csv');
data_srp2 = dlmread('./Data/collagen_storage.csv');



%Extract axial data
stretch1 = data_srp1(:,1) + 1;
stress = data_srp1(:,2);

%Extract shear data
stretch2 = data_srp2(:,1) + 1;
storage = data_srp2(:,2);


norm1 = norm(stress);
norm2 = norm(storage);

factor = (norm1/norm2)^2;

sum = norm1^2 + factor*norm2^2;



%return

%--------------------------------------------------
eps_s = 0; % The critical strain for fiber stiffening
K = 15*1146.282974;

%Prescribe the initial guesses
d0 = 10; 
ds = 10; 
lame = 0;
alpha = 1;


funopt=@(x)sseval(x,K,Mf,wf,stretch1,stretch2,stress,storage); 
ninit = 50;
costarray = zeros(ninit,1);
Parray = zeros(ninit,5);

for i = 1:ninit

d0 = 10*rand;
ds = 10*rand;
lame = 100*rand;
alpha = 100*rand;
eps_s = rand;
x0 = [d0;ds;lame;alpha;eps_s]; % Initial guess

AA = [];
bb = [];
Aeq = [];
beq = [];
lb = [0,0,0,0,0];
ub = [+Inf,+Inf,+Inf,+Inf,1];
[bestx, cost] = fmincon(funopt,x0,AA,bb,Aeq,beq,lb,ub)

costarray(i) = cost;
Parray(i,1) = bestx(1);
Parray(i,2) = bestx(2);
Parray(i,3) = bestx(3);
Parray(i,4) = bestx(4);
Parray(i,5) = bestx(5);


end

[cost_min, I] = min(costarray)
d0 = Parray(I,1)
ds = Parray(I,2)
lame = Parray(I,3)
alpha = Parray(I,4)
eps_s = Parray(I,5)

ds = 1e4;
d0 = 1e4;

s1 = zeros(100,1);
s2 = zeros(100,1);
s3 = zeros(100,1);
eps_norm = zeros(kk,1);
mf = zeros(3,kk);
latx = zeros(100,1);
laty = zeros(100,1);

y0 = [1,1];    

for i=1:100 % plot the stress-strain curve
    
    F33 = 0.5+(i-1)/100; % prescribe the axial stretch
    s1(i) = F33-1; 
    
    F11 = 1;
    F22 = 1;
    latx(i) = F11;
    laty(i) = F22;
    
    
    F1 = [F11 0 0; 0 F22 0; 0 0 F33];
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
    S1 = Pfiber(K,d0,ds,eps_s,F1,F1inv,C1,J1,Mf,mf,eps_norm,wf,lame,alpha);
    s2(i) = S1(3,3);
    
end


plot(s1,s2,'r','LineWidth',6);
xlim([-0.1 0.1]);
hold on;
scatter(stretch1 - 1,stress,200,'r','s');
xlabel('Axial Strain','FontSize',24);
ylabel('Axial Stress (kPa)','FontSize',24);
set(gca,'FontSize',20);
legend({'Model','Experiment'},'FontSize',20,'Location','Northwest');


%-------------------------------------------------------------
%Compute storage modulus as a function of axial strain
t1 = zeros(100,1);
t2 = zeros(100,1);
t3 = zeros(100,1);
gamm = 0.01;


for i=1:100 % plot the stress-strain curve
    
    F33 = 0.8 + (i-1)/100; % prescribe the axial stretch in x1 direction
    t1(i) = F33-1; 
    
    F11 = latx(i);
    F22 = laty(i);

    
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
    S1 = Pfiber(K,d0,ds,eps_s,F1,F1inv,C1,J1,Mf,mf,eps_norm,wf,lame,alpha); 
    
    %modulus = S1(1,3)/(gamm*F33);
    modulus = S1(1,3)/gamm;
    t2(i) = modulus;
    
end

figure;
plot(t1,t2,'b','LineWidth',6);
xlim([-0.1 0.1]);
hold on;
scatter(stretch2 - 1,storage,200,'b','s');
xlabel('Axial Strain','FontSize',24);
ylabel('Storage Modulus (kPa)','FontSize',24);
set(gca,'FontSize',20);
legend({'Model','Experiment'},'FontSize',24,'Location','Northwest');





