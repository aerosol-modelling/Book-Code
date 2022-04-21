% Set the line style for the figures - handy when plotting different simulations
% in the same figure by running the script again when the figures are open
lstyle = '-';               %  '-', '--', '-.', or ':'

% Concentrations as a function of time
figure(1)
semilogy(T/60,C*1e-6,'LineStyle',lstyle)
hold on; set(gcf,'Color','white'); set(gca,'ColorOrderIndex',1)
xlabel('Time (min.)')
ylabel('{\itC} (cm^{-3})')
lgd = legend(string(1:nmax),'Location','best');
title(lgd,'Cluster no.','FontWeight','normal')

% Final concentrations as a function of cluster size
figure(2)
semilogy(1:nmax,C(end,:)*1e-6,'Marker','o','LineStyle',lstyle)
hold on; set(gcf,'Color','white'); set(gca,'ColorOrderIndex',1)
xlabel('Number of molecules in cluster')
xlim([0.5 nmax+0.5]); xticks(1:nmax)
ylabel('Final {\itC} (cm^{-3})')

% Formation rate as a function of time
figure(3)
semilogy(T/60,J*1e-6,'LineStyle',lstyle)
hold on; set(gcf,'Color','white'); set(gca,'ColorOrderIndex',1)
xlabel('Time (min.)')
ylabel('{\itJ}_{out} (cm^{-3} s^{-1})')