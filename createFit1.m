function [fitresult, gof] = createFit1(logws, lmpwrtf)
%CREATEFIT1(LOGWS,LMPWRTF)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : logws
%      Y Output: lmpwrtf
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 20-Nov-2019 17:17:47


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( logws, lmpwrtf );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedPoints = excludedata( xData, yData, 'Indices', [655 656 657 1398] );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
plot( fitresult, xData, yData, excludedPoints );
% Label axes
xlabel logws
ylabel lmpwrtf
grid off


