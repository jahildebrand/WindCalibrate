function [TFFig] = tffigmake(TFFig,TFnew,TFold,tf_file,dBaseName,col,coh)
%
% make TF Fig
%JAh Sept 2020
figure(TFFig);
semilogx(TFnew(10:end,1),TFnew(10:end,2),'r','LineWidth',3); %
hold on
semilogx(TFold(10:end,1),TFold(10:end,2),'b:','LineWidth',3); %
tflegend = replace(tf_file,'_','.');
legend('WindTF',tflegend,'AutoUpdate','off','Location','Northwest');
title([dBaseName,' Hydrophone ',tf_file(1:3)])
xlabel('Frequency [Hz]')
ylabel('Inverse Sensitivity [dB re uPa//counts]')
% add shading
% figure(TFFig)
xl = xlim;
yl = ylim;
gray = [0.8 0.8 0.8];
patch([xl(1) col*100 col*100 xl(1)],[yl(1) yl(1) yl(2) yl(2)],gray)
patch([coh*100 xl(2) xl(2) coh*100],[yl(1) yl(1) yl(2) yl(2)],gray)
set(gca,'children',flipud(get(gca,'children')))
semilogx(TFold(10:end,1),TFold(10:end,2),'b:','LineWidth',3); %
grid on