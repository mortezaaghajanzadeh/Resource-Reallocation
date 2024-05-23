dlmwrite('figures\t_hat_co2.csv',mean(t_hat,1));

t_hat_co   = importdata('figures\t_hat_co.csv').*100;
t_hat_nox  = importdata('figures\t_hat_nox.csv').*100;
t_hat_pm10 = importdata('figures\t_hat_pm10.csv').*100;
t_hat_pm25 = importdata('figures\t_hat_pm25.csv').*100;
t_hat_so2  = importdata('figures\t_hat_so2.csv').*100;
t_hat_voc  = importdata('figures\t_hat_voc.csv').*100;
t_hat_co2  = importdata('figures\t_hat_co2.csv').*100;

nei_yrs = [1990 1996 1999 2002 2005 2008];
nei_yrsIndex = [1 7 10 13 16 19];

plot(nei_yrs,t_hat_co(  nei_yrsIndex),'rd','LineWidth',1.5); hold all;   
plot(nei_yrs,t_hat_nox( nei_yrsIndex),'rv','LineWidth',1.5); hold all;   
plot(nei_yrs,t_hat_pm10(nei_yrsIndex),'rs','LineWidth',1.5); hold all;   
plot(nei_yrs,t_hat_pm25(nei_yrsIndex),'b+','LineWidth',1.5); hold all;   
plot(nei_yrs,t_hat_so2( nei_yrsIndex),'bd','MarkerFaceColor','b','LineWidth',1.5); hold all;   
plot(nei_yrs,t_hat_voc( nei_yrsIndex),'bv','MarkerFaceColor','b','LineWidth',1.5); hold all;   
plot(nei_yrs,t_hat_co2( nei_yrsIndex),'bs','MarkerFaceColor','b','LineWidth',1.5); hold all;   
plot(1990:2008,t_hat_co  ,'--r','LineWidth',1.5); hold all;    
plot(1990:2008,t_hat_nox ,'--r','LineWidth',1.5); hold all;    
plot(1990:2008,t_hat_pm10,'--r','LineWidth',1.5); hold all;    
plot(1990:2008,t_hat_pm25,'--b','LineWidth',1.5); hold all;    
plot(1990:2008,t_hat_so2 ,'--b','LineWidth',1.5); hold all;    
plot(1990:2008,t_hat_voc ,'--b','LineWidth',1.5); hold all;    
plot(1990:2008,t_hat_co2 ,'--b','LineWidth',1.5); hold all;    
set(gca,'FontSize',11); 
set(gca,'YTick',[50 100 150 200 250 300 350 400]); 
axis([1990 2008 50 400]); set(gca,'XTick',[1990 1995 2000 2005 2010]); 
xlabel('Year'); ylabel('1990=100'); box off;
legend('CO','NO_{x}','PM_{10}','PM_{2.5}','SO_{2}','VOCs','CO_{2}','Location','northwest');
saveas(gcf,'figures/f7a.eps','epsc');
close

