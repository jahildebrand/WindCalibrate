function p = getWindParams(varargin)
% JAH 7/2021
% FILE PARAMETERS
% p.harp.Proj = 'SOCAL';
% p.harp.Proj = 'Antarc';
p.harp.Proj = 'GofMX';
% p.harp.Proj = 'OTSG';
% p.harp.Proj = 'OCNMS';
% p.harp.Proj = 'WAT';
p.harp.Proj = 'NRS';
% p.harp.Proj = 'Hawaii';
% p.harp.Site = '06';
% p.harp.Site = '_CINMS_C';
% p.harp.Site = 'HATB';
p.harp.Site = '06';
p.harp.Depl = '14';
% parameters for Wind LTSA v49
p.calavg = 'n'; % calculate LTSA Average ? y or n othewise read in previous
p.calavgm = 'n'; % calculate LTSAm Average ? y or n othewise read in previous
p.calavgl = 'y'; % calculate LTSAl Average ? y or n othewise read in previous
p.use = 'n'; % use full band data ~100kHz and 100 Hz bin
p.usem = 'n'; % use mid band data ~10kHz and 10 Hz bin
p.usel = 'y'; % use low band data ~1kHz and 1 Hz bin
p.RegPlt = 'on'; % show regression plots off or on
p.FrePlt = 'on'; % show freq plots for each sea state off or on
p.NMPlt = 'on'; % show Noise Model plot for each sea state off or on
p.WndPlt = 'on'; % show wind speed plots for each frequency off or on
p.SaveTF = 'yes'; % save  new TF
% p.revise = 'r'; % allow revision of TF cutoff frequencies
p.tf.freq = [];
p.tf.uppc = [];
%
% p.harp.band = 'high'; %  high mid or low
% p.harp.harpDataSummaryCSV = 'L:\Shared drives\Wind_deltaTF\HARPdataSummaryWIND.csv';
% p.harp.harpDataSummaryCSV = 'E:\Wind_TF\HARPdataSummaryWIND.csv';
p.harp.harpDataSummaryCSV = 'H:\Shared drives\MBARC_All\Code\WindCode\HARPdataSummaryWIND.csv';
p.harp.harpDataSummary = readtable(p.harp.harpDataSummaryCSV);
p.harp.WindFolder = ['Z:\Wind_Data\',p.harp.Proj,'\'];
% p.ltsa.LTSAFolder = ['Z:\LTSA\',p.harp.Proj,'\',p.harp.Proj,...
%     '_',p.harp.Site,'_',p.harp.Depl];
p.ltsa.LTSAFolder = ['Z:\LTSA\GofMX\'];
p.tf.TFsFolder = 'H:\Shared drives\MBARC_TF';
% p.tf.TFsFolder = 'Z:\Wind_TF\Output\Kona_TF\TF_Wind'; %uses wind TF
p.tf.TFsFolderOld = 'G:\Harp_TF';
% p.tf.TFsFolderOld = 'G:\Harp_TF\OLD\';
p.harp.OutFolder = 'Z:\Wind_TF\Output\All3_TF';
p.harp.Short = p.harp.Depl;
%
p.harp.NA = 5;     % number of time slices (spectral averages) to read per raw file
p.harp.tres = 2;   % time bin resolution 0 = month, 1 = days, 2 = hours
%
p.harp.OutName1 = [p.harp.Proj,p.harp.Site,p.harp.Depl,...
    '_WindNoise.mat'];
p.harp.OutName2 = [p.harp.Proj,p.harp.Site,p.harp.Depl,...
    '_WindNoise_mid.mat'];
p.harp.OutName3 = [p.harp.Proj,p.harp.Site,p.harp.Depl,...
    '_WindNoise_low.mat'];
%
p.harp.OutWinFig = ['Z:\Wind_TF\Output\All3_TF\VsFreqVsWind\'];
p.harp.OutVsFreq = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_VsFreq.fig'];
p.harp.OutVsFreqm = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_VsFreqm.fig'];
p.harp.OutVsFreql = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_VsFreql.fig'];
p.harp.OutVsWind1 = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_VsWind.fig'];
p.harp.OutVsWind2 = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_VsWind_mid.fig'];
p.harp.OutVsWind3 = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_VsWind_low.fig'];
p.harp.OutTFCorr = [p.harp.Proj,p.harp.Depl,p.harp.Site,...
    '_TFCorr.mat'];
p.harp.OutBad =  ['Z:\Wind_TF\Output\All3_TF\BAD'];