weight1=repmat(bl.Rds0(1,pm.dirty)',[1 Y]);
weight2=repmat(bl.Rds0(1,pm.clean)',[1 Y]);
plot(1990:2008,100.*squeeze(wmean(t_hat(pm.dirty,:),weight1,1)),'-' ,'LineWidth',1.5); hold on;  plot(1990:2008,100.*squeeze(wmean(t_hat(pm.clean,:),weight2,1)),'--','LineWidth',1.5);  set(gca,'FontSize',17);  axis([1990 2011 0 300]);  ylabel('1990=100');  set(gca,'XTick',[1990 1995 2000 2005 2010]);  set(gca,'YTick',[0 100 200 300]);  xlabel('Year');  box off;  legend('Dirty Industries','Clean Industries','Location','SouthEast'); saveas(gcf,'figures\f4.eps','epsc'); close;

