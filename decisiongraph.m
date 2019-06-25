function    [gamma]               =decisiongraph(rho,delta,gamma,SORTE,f)
%%=========================================================================
  figure;
   plot(rho{f},delta{f},'k*','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',8);
   h=xlabel('$$\rho_{\emph{q}}$$');
 set(h,'Interpreter','latex') 
 
    g=  ylabel('$$\delta_{\emph{q}}$$');
 set(g,'Interpreter','latex') 

 maximx=max(rho{f})+0.05;
 maximy=max(delta{f})+0.05;
 axis([0 maximx 0 maximy]);
 grid on;
 %%========================================================================
figure;

plot( gamma{f},'k*','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',8);
 h         =xlabel('$${\emph{q}}$$');
 set(h,'Interpreter','latex') 
 g         =ylabel('$$\gamma_{\emph{q}}$$');
 set(g,'Interpreter','latex') 
 [val,loc]=sort(gamma{f},'descend');
 axis([0 size(loc,2)+2 0 31]);
 grid on;
 %%========================================================================
%  SORTE{f}(5)=0.48
%  SORTE{f}(2)=0.32
 figure;
  plot([2:20],SORTE{f}(2:20),'-k*','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',10);
   h=xlabel('$${\emph{q}}$$');
 set(h,'Interpreter','latex') 
 ylabel('SORTS');
 grid on;
% axis([2, 7, 0,1]);


%===========================ESTIMATED SOURCE SPECTRUM===============================================
for f=79%1:K/2;
dd=1:size(estS,2);
squeeze(estS(f,dd,:)).';
eS=squeeze(abs(estS(f,dd,2)));
ES=squeeze(abs(origsw(1,dd,f)));
err(f)=sum(abs(eS-ES))./size(estS,2);
end
figure;
plot(ES,'k','LineWidth',2);
grid on;
hold on;
plot(eS,'k-.','LineWidth',2);
hold on;
h=xlabel('$${\emph{d}}$$');
set(h,'Interpreter','latex'); 
ylabel('Amplitude');
g=legend('$$|s_{d}|$$','$$|\hat{s}_{d}|$$');
set(g,'Interpreter','latex'); 
set(g,'Fontsize',18);
axis([0,size(estS,2),0,max(ES)+3]);