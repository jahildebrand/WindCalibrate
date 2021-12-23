function [] = calLTSA()
global p
Proj = PARAMS.harp.Proj ;
Site = PARAMS.harp.Site ;
Short = PARAMS.harp.Short ;
Depl = PARAMS.harp.Depl;
band = PARAMS.harp.band; % high mid or low
harpDataSummaryCSV = PARAMS.harp.harpDataSummaryCSV ;
harpDataSummary = PARAMS.harp.harpDataSummary ;
WindFolder = PARAMS.harp.WindFolder ;
LTSAFolder = PARAMS.ltsa.LTSAFolder ;
TFsFolder = PARAMS.tf.TFsFolder ;
TFsFolderOld = PARAMS.tf.TFsFolderOld ;
OutFolder = PARAMS.harp.OutFolder ;
OutName = PARAMS.harp.OutName1;
NA = PARAMS.harp.NA ;     % number of time slices (spectral averages) to read per raw file
tres = PARAMS.harp.tres ;
% Establish location of LTSA files
if exist(fullfile(LTSAFolder,[Proj,Depl,Site,'*.ltsa']),'file')
    ltsafile = fullfile(LTSAFolder,[Proj,Depl,Site,'*.ltsa']);
elseif exist(LTSAFolder,'dir')
    ltsafile = LTSAFolder;
else
    ltsafile = LTSAFolder;
end
[fn_files, fn_pathname] = uigetfile(ltsafile,'Pick LTSA(s)',...
    'MultiSelect','on');
PARAMS.ltsa.LTSAFolder = fn_pathname; % corrects ltsa folder
%
% if strcmp(PARAMS.harp.band,'mid')
%     fn_fileh = fn_files; fn_pathh = fn_pathname;  % save LTSA high names
%     [fn_filem, fn_pathm] = uigetfile('X:\*.ltsa','Pick LTSA(s)',...
%         'MultiSelect','on');
%     ndisk = str2num(fn_pathm(end-7:end-6));
%     fn_pathnew = [fn_pathm(1:end-8),'0',...
%         num2str(ndisk),fn_pathm(end-5:end)];
%     fn_filenew = [fn_filem(1:end-20),'0',...
%         num2str(ndisk),fn_filem(end-17:end)];
%     i = 1;
%     while (exist(fn_pathnew) > 0)
%         fn_files{1,i} = fn_filenew; i = i + 1;
%         fnew = fullfile(fn_pathnew,fn_filenew);
%         copyfile(fnew,fn_pathh) ; % source destination
%         % increment disk number
%         ndisk = ndisk + 1;
%         fn_pathnew = [fn_pathnew(1:end-8),'0',...
%             num2str(ndisk),fn_pathnew(end-5:end)];
%         fn_filenew = [fn_filenew(1:end-20),'0',...
%             num2str(ndisk),fn_filenew(end-17:end)];
%     end
% end
% sort if multiple files selected
sfn = size(fn_files);
if sfn(2) > 1
    fn_files = sort(fn_files);
end
% from LTSA name, determine deployment and PreAmp number and Time Period
deplMatch = [];
iFile = 1;
% if iscell(fn_files)
%     fnFileForTF = fn_files{1};
% else
%     fnFileForTF = fn_files;
% end
while isempty(deplMatch) && iFile<=size(harpDataSummary,1)
    dBaseName = strrep(harpDataSummary.Data_ID{iFile},'-','');
    dBaseName = strrep(dBaseName,'_','');
    if (strfind(lower(dBaseName),lower([Proj,Site,Depl])))
        deplMatch = 1;
    elseif (strfind(lower(dBaseName),lower([Proj,Depl,Site])))
        deplMatch = 1;
    elseif (strfind(lower(dBaseName),lower([Proj,Site,Short,Depl])))
        deplMatch = 1;
    end
    iFile = iFile+1;
end
if isempty(deplMatch)
    dBaseName = [Proj,Site,Depl];
    dBaseName = strrep(dBaseName,'_','-');
    disp('No matching deployment in HARP database)')
    prompt = ' Enter PreAmp Number: ';
    tfNum = input(prompt,'s');
    tfn = str2double(tfNum);
    prompt = ' Enter Depth in m: ';
    depth = input(prompt);
else
    deplMatchIdx = iFile-1;    %Find the TF
    tfNum = harpDataSummary.PreAmp{deplMatchIdx};
    tfn = str2double(tfNum);
    if isnan(tfn)
        prompt = ' Enter PreAmp Number: ';
        tfNum = input(prompt,'s');
        tfn = str2double(tfNum);
    end
    depth = str2double(harpDataSummary.Depth_m{deplMatchIdx});
    if (isnan(depth))
        prompt = ' Enter Depth in m: ';
        depth = input(prompt);
    end
end
%
tfd = floor(tfn/100)*100;
if tfd <= 300
        % TFsFold = [TFsFolder,'\100-399','*'];
    suggestedTFPath = 'L:\Shared drives\MBARC_TF\100-399\320';
else
%     TFsFold = [TFsFolder,num2str(tfd),'_series'];
    TFsFold = [TFsFolder,'\',num2str(tfd),'*'];
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
        TFsFold = [TFsFolderOld,'\',num2str(tfd),'_series'];
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
            suggestedTFPath = fullfile(TFsFold,tfList(tfMatchIdx).name);
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
loadTF(tf); % open and read Transfer function file:
% correct TF for more than one measurement at a specific frequency
[~,ia,ic] = unique(PARAMS.tf.freq);
if length(ia) ~= length(ic)
    disp(['Error: TF file ',tf,' is not monotonically increasing'])
end
tf_freq = PARAMS.tf.freq(ia);
tf_uppc = PARAMS.tf.uppc(ia);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating Averages for:')
if iscell(fn_files) % then it's a cell full of filenames
    fn = cell(length(fn_files),1);
    for k = 1:length(fn_files)
        fn{k} = fullfile(fn_pathname, fn_files{k});
        disp(fn{k});
    end
else % it's just one file, but put into cell
    fn = cell(1,1);
    fn{1} = fullfile(fn_pathname,fn_files);
    disp(fn);
end
nltsas = length(fn);    % number of LTSA files
mnum2secs = 24*60*60;

% read ltsa headers and sum the total number of raw files for
% pre-allocating vectors/matrices
nrftot = 0;
for k = 1:nltsas    % loop over files
    PARAMS.ltsa = [];   % clear
    PARAMS.ltsa.ftype = 1;
    [PARAMS.ltsa.inpath,infile,ext] = fileparts(fn{k});
    PARAMS.ltsa.infile = [infile,ext];
    disp(infile)
    read_ltsahead  % better would be to just read nrftot instead of whole header
    nrf = PARAMS.ltsa.nrftot;
    nrftot = nrftot + nrf;
    if k == 1 % read and set some useful parameters that should be the same across all ltsas
        nf = PARAMS.ltsa.nf;
        nave = PARAMS.ltsa.nave(1);
        freq = PARAMS.ltsa.freq;
        dfreq = PARAMS.ltsa.dfreq;
        fs0 = PARAMS.ltsa.fs;
    end
end

% Transfer function correction vector
Ptf = interp1(tf_freq,tf_uppc,freq,'linear','extrap');

% fill up header matrix H
PARAMS.ltsa = [];   % clear
H = zeros(nrftot,3);    % 3 columns: filenumber, datenumber, byteloc in filenumber
cnt1 = 1;
cnt2 = 0;
doff = datenum([2000 0 0 0 0 0]);   % convert ltsa time to millenium time
eltsa = [];
for k = 1:nltsas    % loop over files
    PARAMS.ltsa.ftype = 1;
    [PARAMS.ltsa.inpath,infile,ext] = fileparts(fn{k});
    PARAMS.ltsa.infile = [infile,ext];
    read_ltsahead  % better would be to just read nrftot instead of whole header
    nrf = PARAMS.ltsa.nrftot;
    fs0 = PARAMS.ltsa.fs;
    knave = PARAMS.ltsa.nave;
    cnt2 = cnt2 + nrf;
    H(cnt1:cnt2,1) = k.*ones(nrf,1);
    H(cnt1:cnt2,2) = PARAMS.ltsa.dnumStart + doff;  % add doff JAH to make real datenum
    H(cnt1:cnt2,3) = PARAMS.ltsa.byteloc;
    cnt1 = cnt2 + 1;
    eltsa(k) = cnt2;
end

dvec = datevec(H(:,2));

if tres == 0
    mnum = dvec(:,1).*12 + dvec(:,2);   % month number where 1 = Jan 2000
elseif tres == 1
    mnum = floor(datenum(dvec));   % day number where 1 = Jan 2000
elseif tres == 2
    mnum = floor(datenum(dvec) * 24);   % hour
else
    disp(['Error: unknown time resolution = ',num2str(tres)])
end
mnumMin = min(mnum);
mnumMax = max(mnum);
%
dur = unique(mnum);   % unique averaging time bins
nm = length(dur); % number of ave time bins
cnt1 = 1;
cnt2 = 0;
ptime = zeros(nm,1);      % start time of bin average
nmave = zeros(nm,2);    % number of averages possible and used(ie not filtered out) for each bin
mpwr = ones(nf,nm);     % mean power over time period
mpwrtf = ones(nf,nm);   % mean power with TF applied
% mpwrTF = ones(nf,nm);   % mean power with FIFO interp & TF applied
eltsam = 1; Cltsa = 0;
for m = 1:nm    % loop over time average bins
    I = [];
    I = find(mnum == dur(m));
    nrfM = length(I);
    pwrM = [];
    pwrM = zeros(nf,NA*nrfM);
    fnum = [];
    fnum = unique(H(I,1));
    nfiles = length(fnum);
    
    NBO = 0;    % number of averages (taves) read for mean spectra
    for f = 1:nfiles    % loop over files with same month (probably only 2 max)
        % open ltsa file
        fid = fopen(fn{fnum(f)},'r');
        % samples to skip over in ltsa file
        J = [];
        J = find(H(I,1) == fnum(f));
        nrfRead = length(J);
        skip = H(I(J(1)),3);    % get first byteloc of file for that month
        fseek(fid,skip,-1);    % skip over header + other data
        
        ptime(m) = H(I(J(1)),2);     % save 1st time of this time unit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of time slices to use for calcs need to set 'xxxx*int8' value
        prcsn = [num2str(NA*nf),'*int8'];
        
        % allocate memory
        pwrF = [];
        NB = NA * nrfRead;
        pwrF = zeros(nf,NB);
        
        % number of time slices to skip on each read
        SA = nave - NA;
        skip = nf*(SA);
        % initial time slice to start at
        if fs0 == 2000 || fs0 == 200000 || fs0 == 10000 || fs0 == 64000 || fs0 == 50000
            %  IA = (15+5)*nf;
            IA = (nave + NA)*nf;  % better than commented 200kHz hardwired above
        elseif fs0 == 3200 || fs0 == 320000
            IA = (nave + NA-1)*nf;
        else
            disp('Error: unknown sample rate')
            disp(['fs0 = ',num2str(fs0)])
        end
        fseek(fid,IA,-0);
        
        % read data into File power
        [pwrF,count] = fread(fid,[nf,NB],prcsn,skip);
        if count ~= nf*NB
            disp('error = did not read enough data')
            %             NBO = NB;
            NB = floor(count/nf);
            %             disp(['NBO = ',num2str(NBO)])
            disp(['NB = ',num2str(NB)])
        end
        
        % fill up power
        if nfiles == 1
            pwrM = pwrF;
            NBO = NB;
        else
            if f == 1
                pwrM(1:nf,1:NB) = pwrF;
                NBO = NB;
            else
                pwrM(1:nf,NBO+1:NBO+NB) = pwrF;
                NBO = NBO+NB;
            end
        end
        %disp(['Number of 5s time bins in this Average Time Bin = ',num2str(NBO)]);
        fclose(fid);
        [Cltsa, ialtsa, ibltsa] =  intersect(I,eltsa);
        if Cltsa > 0
            disp([num2str(m),' end of ltsa ',fn{fnum(f)}]);
            eltsam = [eltsam; m];
            Cltsa = 0;
        end
    end  % end for f
    
    nmave(m,1) = NBO;
    %     nmave(m,2) = NBF;
    
    %     cnt2 = cnt2 + size(pwrN,2);
    cnt2 = cnt2 + size(pwrM,2);
    %     pwrA(1:nf,cnt1:cnt2) = pwrM;
    cnt1 = cnt2 + 1;
    
    % mean - these are smooth (floating point pwr values)
    %     mpwr(1:nf,m) = mean(pwrN,2);
    mpwr(1:nf,m) = mean(pwrM,2);
    mpwrtf(1:nf,m) = mpwr(1:nf,m) + Ptf';   % add transfer function
    
end
%REMOVE TIMES outside of 2006 - 2020
d2006 = datenum([2005 0 0 0 0 0]);
d2019 = datenum([2021 0 0 0 0 0]);
gotime = find(ptime < d2019 & ptime > d2006);
ptime = ptime(gotime);
mpwr = mpwr(:,gotime);
mpwrtf = mpwrtf(:,gotime);
%
% Get Wind data
syr(1,:) = datevec(ptime(1));
syr(2,:) = datevec(ptime(end));
%
vec = []; speed = [];
for i = syr(1,1) : syr(2,1)
    if exist([WindFolder,Proj,Site,'windvec.mat'])
        load([WindFolder,Proj,Site,'windvec.mat'],'dnvec','wspeed');
        disp([' Load ',WindFolder,Proj,Site,'windvec.mat'])
    elseif exist([WindFolder,num2str(i),'\',Proj,Site,'windvec.mat'])
        load([WindFolder,num2str(i),'\',Proj,Site,'windvec.mat'],'dnvec','wspeed');
        disp([' Load ',WindFolder,num2str(i),'\',Proj,Site,'windvec.mat'])
    elseif exist([WindFolder,num2str(i),'\',Site,'windvec.mat'])
        load([WindFolder,num2str(i),'\',Site,'windvec.mat'],'dnvec','wspeed');
         disp([' Load ',WindFolder,num2str(i),'\',Site,'windvec.mat'])
    elseif exist([WindFolder,num2str(i),'\',Proj,Short,'windvec.mat'])
        load([WindFolder,num2str(i),'\',Proj,Short,'windvec.mat'],'dnvec','wspeed');
        disp([' Load ',WindFolder,num2str(i),'\',Proj,Short,'windvec.mat'])
    elseif exist([WindFolder,num2str(i),'\',Site,Short,'windvec.mat'])
        load([WindFolder,num2str(i),'\',Site,Short,'windvec.mat'],'dnvec','wspeed');
        disp([' Load ',WindFolder,num2str(i),'\',Site,Short,'windvec.mat'])
    else
        disp('Could Not Find Wind File');
        suggestedWind = fullfile(WindFolder,num2str(i));
        [w_file, w_pathname ] =...
            uigetfile(fullfile(suggestedWind,'*.mat'),'Pick Wind File');
        load([w_pathname,w_file],'dnvec','wspeed');
    end
    if strcmp([Proj,Site],'SOCALB') || strcmp([Proj,Site],'SOCAL_CINMS_B') 
        dnew = dnvec;  % all years in one file
        wsnew = wspeed;
    else
        vec = [vec;dnvec];
        speed = [speed;wspeed];
    end
end
% re-interpolate for hourly from every 6 hr
if ~strcmp([Proj,Site],'SOCALB') && ~strcmp([Proj,Site],'SOCAL_CINMS_B')
    dnew = []; wsnew = [];
    for i = 1 : length(vec)
        dnew((i-1)*6 +1) = vec(i);
        wsnew((i-1)*6 +1) = speed(i);
        if i < length(vec)
            for k = 1 : 5
                dnew((i-1)*6 +1 + k) = (k/6)  * vec(i+1) + ...
                    ((6-k)/6) * vec(i);
                wsnew((i-1)*6 +1 + k) = (k/6)  * speed(i+1) + ...
                    ((6-k)/6) * speed(i);
            end
        end
    end
end
% save results for noise and wind
PARAMS.ltsa.LTSAFolder = fn_pathname; % corrects ltsa folder
save(fullfile(OutFolder,'WindNoiseMat',OutName)) % save all workspace
end