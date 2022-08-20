data1 = dlmread('./output0_w1w2=1.txt','',1,0);
stretch = data1(:,1);
stress_num = data1(:,3);
stress_LCC = data1(:,4);

plot(stretch,stress_num,'r','LineWidth',3);
hold on;
plot(stretch,stress_LCC,'b:','LineWidth',3);

legend({'Numerical derivative','LCC'},'FontSize',24);

xlabel('Stretch','LineWidth',24);
ylabel('PK-Stress','LineWidth',24);
xlim([0.4 1.6])

set(gca,'FontSize',24);
