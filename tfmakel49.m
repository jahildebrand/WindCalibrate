function [TFold, TFnew] = tfmakel49(dfreql,nfl,col, coh, freq, Ptf, MTFCorr)
% create new TF for low and mid
%JAH Sept 2020
iPtf = find(freq >= 10 & freq <= (nfl-1)*dfreql); % part below between 10 and end Hz
if10 = find(freq < 10 + dfreql/2 & freq > 10 - dfreql/2); % 10 Hz
TFnew = [freq(iPtf)',Ptf(iPtf)'];% adds 10 Hz to nfl-1 Hz data
TFold = TFnew;
TFnew(1:col-if10,2) = TFnew(1:col-if10,2) + MTFCorr(col)'; % for < col
TFnew(col-if10+1:coh-if10+1,2) = TFnew(col-if10+1:coh-if10+1,2) + MTFCorr(col-if10+1:coh-if10+1)'; % for col < i < coh
liPtf = length(iPtf);
if liPtf > coh-if10+1
    TFnew(coh-if10+2:liPtf,2) = TFnew(coh-if10+2:liPtf,2) + MTFCorr(coh-if10+1); % for < col
end
end
%