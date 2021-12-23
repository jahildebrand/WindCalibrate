function [TFold, TFnew] = tfmake49(dfreq,nf,col, coh, freq, Ptf, MTFCorr,tfn)
% create new TF
%JAH Sept 2020
global p
iPtf = find(p.tf.freq > 0 & p.tf.freq < 100); % Get TF below 100 Hz
miP = max(iPtf);
TFnew = [p.tf.freq(iPtf)',p.tf.uppc(iPtf)'];% adds < 100 Hz data
TFnew(miP+1:miP+nf-1,:) = [freq(2:nf)',Ptf(2:nf)']; %  100 Hz - 100 kHz
TFold = TFnew;
TFnew(1:miP+col-1,2) = TFnew(1:miP+col-1,2) + MTFCorr(col); % for 100 Hz - 400 Hz
TFnew(miP+col:miP+coh,2) = TFnew(miP+col:miP+coh,2) + MTFCorr(col:coh)'; % for 500 Hz - 20 kHz
% if (tfn >= 697 && tfn < 780 ) % 20 kHz xover
%     for i = 1:80
%     TFnew(miP+coh+i,2) = TFnew(miP+coh+i,2) + (1 - i/100) * MTFCorr(coh); % taper over 8 kHz
%     end
% else
TFnew(miP+coh+1:end,2) = TFnew(miP+coh+1:end,2) + MTFCorr(coh); % for > 20 kHz
end
%