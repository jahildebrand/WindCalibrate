function [ptime,mpwr,mpwrtf,freq,eltsa,dfreq,dnew,wsnew,fs0] = calLTSA49()
global p PARAMS
Proj = p.harp.Proj ;
Site = p.harp.Site ;
Short = p.harp.Short ;
Depl = p.harp.Depl;
WindFolder = p.harp.WindFolder ;
LTSAFolder = p.ltsa.LTSAFolder ;
OutFolder = p.harp.OutFolder ;
OutName = p.harp.OutName1;
NA = p.harp.NA ;     % number of time slices (spectral averages) to read per raw file
tres = p.harp.tres ;
% Establish location of LTSA files
if exist(fullfile(LTSAFolder,[Proj,Depl,Site,'*.ltsa']),'file')
    ltsafile = fullfile(LTSAFolder,[Proj,Depl,Site,'*.ltsa']);
elseif exist(LTSAFolder,'dir')
    ltsafile = LTSAFolder;
else
    ltsafile = LTSAFolder;
end
[fn_files, fn_pathname,~] = uigetfile('*.ltsa','Pick LTSA(s)',LTSAFolder,...
    'MultiSelect','on');
p.ltsa.LTSAFolder = fn_pathname; % corrects ltsa folder
%
% get TF
tf_freq = p.tf.freq;
tf_uppc = p.tf.uppc;
% sort if multiple files selected
sfn = size(fn_files);
if sfn(2) > 1
    fn_files = sort(fn_files);
end

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
         p.ltsa.freq = freq;
        dfreq = PARAMS.ltsa.dfreq;
        p.ltsa.dfreq =dfreq;
        fs0 = PARAMS.ltsa.fs;
        p.ltsa.fs0 = fs0;
    end
end

% Transfer function correction vector
Ptf = interp1(tf_freq,tf_uppc,freq,'linear','extrap');
% fill up header matrix H
PARAMS.ltsa = [];  % clear
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
eltsa = 1; Cltsa = 0;
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
            %disp('error = did not read enough data')
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
        [Cltsa, ~, ~] =  intersect(I,eltsa);
        if Cltsa > 0
            disp([num2str(m),' end of ltsa ',fn{fnum(f)}]);
            eltsa = [eltsa; m];
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
d2021 = datenum([2021 0 0 0 0 0]);
gotime = find(ptime < d2021 & ptime > d2006);
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