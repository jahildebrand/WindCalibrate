function [TFFig] = tffigmakel(TFFig,TFnew,TFold,tf_file,dBaseName,col,coh)
%
% make TF Fig low
%JAh Sept 2020
figure(TFFig);
semilogx(TFnew(1:end-1,1),TFnew(1:end-1,2),'r','LineWidth',3); %
hold on
semilogx(TFold(1:end,1),TFold(1:end,2),'b:','LineWidth',3); %
%
tflegend = replace(tf_file,'_','.');
legend('WindTF',tflegend,'AutoUpdate','off','Location','Northwest');
title([dBaseName,' Hydrophone ',tf_file(1:3)])
xlabel('Frequency [Hz]')
ylabel('Inverse Sensitivity [dB re uPa//counts]')
% add shading
% figure(TFFig)
yl = ylim;
v = [10 TFnew(end,1) yl(1) yl(2)];
axis(v);
xl = xlim;
gray = [0.8 0.8 0.8];
patch([xl(1) col*10 col*10 xl(1)],[yl(1) yl(1) yl(2) yl(2)],gray)
patch([coh*10 xl(2) xl(2) coh*10],[yl(1) yl(1) yl(2) yl(2)],gray)
set(gca,'children',flipud(get(gca,'children')))
semilogx(TFold(1:end,1),TFold(1:end,2),'b:','LineWidth',3); %
grid on