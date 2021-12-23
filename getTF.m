function p = getTF()
% find TF number from XLS and then read in TF 
global p PARAMS HANDLES
% from LTSA name, determine deployment and PreAmp number and Time Period
deplMatch = [];
iFile = 1;

% if iscell(fn_files)
%     fnFileForTF = fn_files{1};
% else
%     fnFileForTF = fn_files;
% end
while isempty(deplMatch) && iFile<=size(p.harp.harpDataSummary,1)
    p.harp.dBaseName = strrep(p.harp.harpDataSummary.Data_ID{iFile},'-','');
    p.harp.dBaseName = strrep(p.harp.dBaseName,'_','');
    if (strfind(lower(p.harp.dBaseName),lower([p.harp.Proj,p.harp.Site,p.harp.Depl])))
        deplMatch = 1;
    elseif (strfind(lower(p.harp.dBaseName),lower([p.harp.Proj,p.harp.Depl,p.harp.Site])))
        deplMatch = 1;
%     elseif (strfind(lower(p.harp.dBaseName),lower([p.harp.Proj,p.harp.Site,p.harp.Short,p.harp.Depl])))
%         deplMatch = 1;
    end
    iFile = iFile+1;
end
if isempty(deplMatch)
    p.harp.dBaseName = [p.harp.Proj,p.harp.Site,p.harp.Depl];
    p.harp.dBaseName = strrep(p.harp.dBaseName,'_','-');
    disp('No matching deployment in HARP database)')
    prompt = ' Enter PreAmp Number: ';
    tfNum = input(prompt,'s');
    p.tf.tfn = str2double(tfNum);
    prompt = ' Enter Depth in m: ';
    p.tf.depth = input(prompt);
else
    deplMatchIdx = iFile-1;    %Find the TF
    tfNum = p.harp.harpDataSummary.PreAmp(deplMatchIdx);
    p.tf.tfn = tfNum;
    if isnan(p.tf.tfn)
        prompt = ' Enter PreAmp Number: ';
        tfNum = input(prompt,'s');
        p.tf.tfn = str2double(tfNum);
    end
    p.tf.depth = (p.harp.harpDataSummary.Depth_m(deplMatchIdx));
    if (isnan(p.tf.depth))
        prompt = ' Enter Depth in m: ';
        p.tf.depth = input(prompt);
    end
end
%
tfd = floor(p.tf.tfn/100)*100;
if tfd <= 300
        % TFsFold = [TFsFolder,'\100-399','*'];
    suggestedTFPath = 'P:\Shared drives\MBARC_TF\100-399\320';
else
%     TFsFold = [TFsFolder,num2str(tfd),'_series'];
    TFsFold = [p.tf.TFsFolder,'\',num2str(tfd),'*'];
    %Search TFs folder for the appropriate preamp
    tfdir = dir(TFsFold);
    tfList = dir(fullfile(tfdir.folder,tfdir.name));
    tfMatch = [];
    iTF = 1;
    while isempty(tfMatch) && iTF<=size(tfList,1)
        tfMatch = strfind(tfList(iTF).name,tfNum);
        iTF = iTF+1;
    end
    if ~isempty(tfMatch)
        tfMatchIdx = iTF-1;
        suggestedTFPath = fullfile(tfList(tfMatchIdx).folder,tfList(tfMatchIdx).name);
    else
%         try Old TF versions
        warning('No matching New TF in new TFs folder)')
        TFsFold = [p.tf.TFsFolderOld,'\',num2str(tfd),'_series'];
%         Search TFs folder for the appropriate preamp
        tfList = dir(TFsFold);
        tfMatch = [];
        iTF = 1;
        while isempty(tfMatch) && iTF<=size(tfList,1)
            tfMatch = strfind(tfList(iTF).name,tfNum);
            iTF = iTF+1;
        end
        if ~isempty(tfMatch)
            tfMatchIdx = iTF-1;
            suggestedTFPath = fullfile(p.tf.TFsFold,tfList(tfMatchIdx).name);
        end
    end
end
if ~exist('suggestedTFPath','var')
     suggestedTFPath = 'F:\Shared drives\MBARC_TF\';
end
[tf_file, tf_pathname ] = ...
    uigetfile(fullfile(suggestedTFPath,'*.tf'),'Pick Transfer Function');

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([ 'Transfer function: ' tf_file  ]);
tf = fullfile(tf_pathname, tf_file);
loadTFwind(tf); % open and read Transfer function file:
% correct TF for more than one measurement at a specific frequency
[~,ia,ic] = unique(p.tf.freq);
if length(ia) ~= length(ic)
    disp(['Error: TF file ',tf,' is not monotonically increasing'])
end
p.tf.freq = p.tf.freq(ia);
p.tf.uppc = p.tf.uppc(ia);
p.tf.tffile = tf_file;
