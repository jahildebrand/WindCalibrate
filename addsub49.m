function [zTD] = addsub49(addsub,zTD)
% Get data to add or subtract from Wind v Noise plot
%JAH 12-21
if strcmp(addsub,'a')
    pl = selectdataA('selectionmode','brush');
    zTD = [zTD, pl'];
    zTD = unique(zTD); % remove duplicate
end
if strcmp(addsub,'s')
    pl = selectdataS('selectionmode','brush');
    zTD = setdiff(zTD, pl);
end