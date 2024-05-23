nei_yrs = [1990 1996 1999 2002 2005 2008];
nei_yrsIndex = [1 7 10 13 16 19];

plot(1990:2008,Z_hat_o(:,mainpoll)   ,'-','LineWidth',2.5); hold all;    
plot(nei_yrs,Z_hat_o_cf(nei_yrsIndex,1),'rp','LineWidth',1.5); hold all; 
plot(nei_yrs,Z_hat_o_cf(nei_yrsIndex,2),'ro','LineWidth',1.5); hold all; 
plot(nei_yrs,Z_hat_o_cf(nei_yrsIndex,3),'rs','LineWidth',1.5); hold all; 
plot(nei_yrs,Z_hat_o_cf(nei_yrsIndex,4),'rv','LineWidth',1.5); hold all; 
plot(1990:2008,Z_hat_o_cf(:,1),'--r','LineWidth',1.5); hold all;    
plot(1990:2008,Z_hat_o_cf(:,2),'--r','LineWidth',1.5); hold all;    
plot(1990:2008,Z_hat_o_cf(:,3),'--r','LineWidth',1.5); hold all;    
plot(1990:2008,Z_hat_o_cf(:,4),'--r','LineWidth',1.5); hold all;    
set(gca,'FontSize',13); 
set(gca,'YTick',[0 30 60 90 120 150]); 
axis([1990 2008 0 150]); set(gca,'XTick',[1990 1995 2000 2005 2010]); 
xlabel('Year'); ylabel('1990=100'); box off;
if mainpoll == 1 | mainpoll == 7
    legend('Actual Data (All Shocks)','Foreign Competitiveness Shocks Only','U.S. Competitiveness Shocks Only','U.S. Regulation Shocks Only','U.S. Expenditure Share Shocks Only','Location','Southwest')
end

cd 'figures'
filename =  ['cf_' pollname '.eps'];
saveas(gcf,filename,'epsc')
close

filename =  ['t_hat_' pollname '.csv'];
dlmwrite(filename,wmean(t_hat,repmat(bl.Rds0(1,:)',[1 Y]),1));

cd ..

