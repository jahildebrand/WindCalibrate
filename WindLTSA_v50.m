% Compare Wind to LSTA - SINGLE DEPLOYMENT
% v45 - add mid band  11/2020 - v47 no need to move midband 1/2021 skip48
% v49 add low band 7/2021
% v50 reorganization to allow use of Low, Mid or High independently
% Low = 1 Hz , Mid = 10 Hz, High = 100 Hz bin in LTSA
% JAH 10/2019 begin
% derived from LTSAdailySpectra.m 141103 smw
%% Parameters
clear variables
close all
% p holds var that come from getWindParams
% PARAMS read from the LTSAs
global p
p = getWindParams; % paramter file
% Get TF - allows update without recalculating LTSA averages
p = getTF; % Save tf incase it gets overwritten in recall
Psave = p;

%% calculate or load WindNoise.mat
% Low band assumes 1 Hz LTSA bins
if strcmp(Psave.usel,'y')
    if strcmp(Psave.calavgl,'y') % calculate WindNoise_low.mat file ?
        [ptimel,mpwrl,mpwrtfl,freql,eltsal,dfreql,dnew,wsnew,fs0l] = calLTSAl49;  % reads LTSA files to get averages and gets Wind model
    else
        % if low exists load it
        if exist(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName3))
            disp(['load: ',p.harp.OutName3])
            load(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName3),...
                'ptimel','mpwrl','mpwrtfl','freql','eltsal','dfreql','dnew','wsnew','p','zTDl','fs0l');
            % check if tf is correct
            if p.tf.uppc ~= Psave.tf.uppc
                disp('Need to update low TF')
                return
            end
        else
            disp('No existing Low WindNoise file')
            return
        end
    end
    eltsal = unique(eltsal);
end
% Mid band assumes 10 Hz LTSA bins
if strcmp(Psave.usem,'y')
    if strcmp(Psave.calavgm,'y') % calculate WindNoise_mid.mat file ?
        [ptimem,mpwrm,mpwrtfm,freqm,eltsam,dfreqm,dnew,wsnew,fs0m] = calLTSAm49;  % reads LTSAm files to get averages and gets Wind model
    else
        % if mid exists load it
        if exist(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName2))
            disp(['load: ',p.harp.OutName2])
            load(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName2),...
                'ptimem','mpwrm','mpwrtfm','freqm','eltsam','dfreqm','dnew','wsnew','p','zTDm','fs0m');
            if p.tf.uppc ~= Psave.tf.uppc
                disp('Need to update mid TF')
                return
            end
        else
            disp('No existing Mid WindNoise file')
            return
        end
    end
    eltsam = unique(eltsam);
end
%High band assumes 100 Hz LTSA bins
if strcmp(Psave.use,'y')
    if strcmp(Psave.calavg,'y') % calculate WindNoise.mat file ?
        [ptime,mpwr,mpwrtf,freq,eltsa,dfreq,dnew,wsnew,fs0] = calLTSA49();  % reads LTSA files to get averages and gets Wind model
        load(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName1));
    else
        if exist(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName1))
            disp(['load: ',p.harp.OutName1])
            load(fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName1),...
                'ptime','mpwr','mpwrtf','freq','eltsa','dfreq','dnew','wsnew','p','zTD','fs0');
            if p.tf.uppc ~= Psave.tf.uppc
                disp('Need to update high TF')
                return
            end
        else
            disp('No existing High WindNoise file')
            return
        end
    end
    eltsa = unique(eltsa);
end
Psave.ltsa =  p.ltsa;
p = Psave;
depth = p.tf.depth;
tf_file = p.tf.tffile;

%% Noise model
% Ocean Wind Noise Model
a=2.8; b=600; aof=0; bof=100; cof=12; %Knudsen 8 as in JASA
nfacl = 1000; mfacl = 1000; mfacf = 150; mfaca = 3;%
[kdp,ss,fnm,ms] = NoiseModelnew(depth,a,b,aof,bof,cof,nfacl,mfacl,mfacf,mfaca);  % Knudsen8
fm = .01 : .01 : 5; % 10 Hz - 5 khz in 10 Hz steps
kdm = zeros(11,500);
for i = 1:11
    kdm(i,:) = interp1([fnm(1:10),fnm(12:60)],...
        [kdp(i,1:10),kdp(i,12:60)],fm); %10 Hz -5 kHz
end
kd = kdp(:,11:end);
% plot Knudsen curves
if strcmp(p.NMPlt,'on')
    nmFig = figure(5); clf;
    for i = 1: length(ms)
        semilogx(1000*fnm, kdp(i,:),'k','LineWidth',2);
        if i == 1
            hold on
        end
    end
    i=5;
    semilogx(1000*fnm, kdp(i,:),'r','LineWidth',2); % ss = 4 is log10(ms) = 1
    axis([10,160000,10,95]);
    xlabel('Frequency [Hz]')
    ylabel('dB re uPa^2/Hz')
    ttitle = ['Noise Model ',num2str(depth),' m ',...
        p.harp.dBaseName,' Hyd ',tf_file(1:3)];
    title(ttitle)
    grid on
    hold on
end
%% Make Wind vs Noise plots for cleaning data
% reduce sig figures to make wind and noise times match, accurate to ~ 7 min
xw = round(dnew' .* 100)./100;
% High Frequency = 100 Hz bin
if strcmp(Psave.use,'y')
%     if strcmp(Psave.calavg,'y')
%         fs0 = p.ltsa.fs0;
%     end
    xp = round(ptime .* 100)./100; % reduce sig fig to make match
    [~,inoise,iwind] = intersect(xp,xw);
    iX = 1: length(iwind); % in case there is no mid or low
    if fs0 >= 200000
        if60k = find(freq < 60000 + dfreq/2 & freq > 60000 - dfreq/2); % 60Khz
        if70k = find(freq < 70000 + dfreq/2 & freq > 70000 - dfreq/2); % 70Kh
    end
    if fs0 >= 100000
        if30k = find(freq < 30000 + dfreq/2 & freq > 30000 - dfreq/2); % 30Khz
        if40k = find(freq < 40000 + dfreq/2 & freq > 40000 - dfreq/2); % 40Khz
        if50k = find(freq < 50000 + dfreq/2 & freq > 50000 - dfreq/2); % 50Khz
    end
    if fs0 >= 48000
        if200 = find(freq < 200 + dfreq/2 & freq > 200 - dfreq/2); % 200 Hz
        if300 = find(freq < 300 + dfreq/2 & freq > 300 - dfreq/2); % 300 Hz
        if500 = find(freq < 500 + dfreq/2 & freq > 500 - dfreq/2); % 500 Hz
        if1000 = find(freq < 1000 + dfreq/2 & freq > 1000 - dfreq/2); % 1Khz
        if2000 = find(freq < 2000 + dfreq/2 & freq > 2000 - dfreq/2); % 2Khz
        if10k = find(freq < 10000 + dfreq/2 & freq > 10000 - dfreq/2); % 10Khz
        if20k = find(freq < 20000 + dfreq/2 & freq > 20000 - dfreq/2); % 20Khz
    end
end
% Mid Frequency = 10 Hz bin
if strcmp(Psave.usem,'y')
%     if strcmp(Psave.calavgm,'y')
%         fs0m = p.ltsa.fs0m;
%     end
    xpm = round(ptimem .* 100)./100;
    [~,inoisem,iwindm] = intersect(xpm,xw);
    if strcmp(Psave.use,'y')
        [~,iX,iXm] = intersect(iwind,iwindm);% make wind agree for mid and high
    else
        iXm = 1 : length(iwindm);
    end
    if20m = find(freqm < 20 + dfreqm/2 & freqm > 20 - dfreqm/2); % 20 Hz
    if50m = find(freqm < 50 + dfreqm/2 & freqm > 50 - dfreqm/2); % 50 Hz
    if100m = find(freqm < 100 + dfreqm/2 & freqm > 100 - dfreqm/2); % 100 Hz
    if200m = find(freqm < 200 + dfreqm/2 & freqm > 200 - dfreqm/2); % 200 Hz
    if500m = find(freqm < 500 + dfreqm/2 & freqm > 500 - dfreqm/2); % 500 Hz
    if1000m = find(freqm < 1000 + dfreqm/2 & freqm > 1000 - dfreqm/2); % 1000 Hz
end
% Low frequency = 1 Hz bin
if strcmp(Psave.usel,'y')
%     if strcmp(Psave.calavgl,'y')
%         fs0l = p.ltsa.fs0l;
%     end
    xpl = round(ptimel .* 100)./100;
    [~,inoisel,iwindl] = intersect(xpl,xw);
    if strcmp(Psave.usem,'y')
        [~,~,iXl] = intersect(iwindm,iwindl);% make wind agree for low and mid
    else
        iXl = 1 : length(iwindl);
    end
    if5l = find(freql < 5 + dfreql/2 & freql > 5 - dfreql/2); % 5 Hz
    if10l = find(freql < 10 + dfreql/2 & freql > 10 - dfreql/2); % 10 Hz
    if20l = find(freql < 20 + dfreql/2 & freql > 20 - dfreql/2); % 20 Hz
    if50l = find(freql < 50 + dfreql/2 & freql > 50 - dfreql/2); % 50 Hz
    if100l = find(freql < 100 + dfreql/2 & freql > 100 - dfreql/2); % 100 Hz
    if200l = find(freql < 200 + dfreql/2 & freql > 200 - dfreql/2); % 200 Hz
    if500l = find(freql < 500 + dfreql/2 & freql > 500 - dfreql/2); % 500 Hz
    if1000l = find(freql < 1000 + dfreql/2 & freql > 1000 - dfreql/2); % 1000 Hz
end

%% make vs wind figures with fit to eliminate bad data
isok = cell(1,9); % array to hold edited data
% figure(2)
figure(2); clf; set(2,'name',sprintf('Wind vs Noise'));
set(gcf,'position',[20 500 600 450]);
%
if strcmp(p.usel,'y') && strcmp(p.calavgl,'y') % low cleaning
    h1 = subplot(3,3,1);
    grid on; hold on;
    plot(h1,wsnew(iwindl(iXl)),mpwrtfl(if100l,inoisel(iXl)),'ro'); %use 200 Hz
    [isok{1,1}] = createFit4(wsnew(iwindl(iXl)),mpwrtfl(if100l,inoisel(iXl)),...
        wsnew(iwindl(iXl)),mpwrtfl(if100l,inoisel(iXl)),h1 );
    legend off
    xlabel('Wind Speed m/s');
     ltex = sprintf('%s\n%s','Low Band','100 Hz');
    ylabel(ltex);
    title([p.harp.dBaseName,'Noise vs Wind Speed']);
    %
    h2 = subplot(3,3,2);
    plot(h2,wsnew(iwindl(iXl)),mpwrtfl(if200l,inoisel(iXl)),'ro'); %use 500 Hz
    ltex = '200 Hz';
    grid on; hold on;
    [isok{1,2}] = createFit5(ltex,wsnew(iwindl(iXl)),mpwrtfl(if200l,inoisel(iXl)),...
        wsnew(iwindl(iXl)),mpwrtfl(if200l,inoisel(iXl)),h2 );
    legend off
    %
    h3 = subplot(3,3,3);
    plot(h3,wsnew(iwindl(iXl)),mpwrtfl(if500l,inoisel(iXl)),'ro'); %use 500 Hz
    ltex='500 Hz';
    grid on; hold on;
    [isok{1,3}] = createFit5(ltex,wsnew(iwindl(iXl)),mpwrtfl(if500l,inoisel(iXl)),...
        wsnew(iwindl(iXl)),mpwrtfl(if500l,inoisel(iXl)),h3 );
    legend off
    if ~exist('zTDl','var')
        zTDl = mintersect(isok{1,1}, isok{1,2}, isok{1,3});
        wnfilel = fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName3);
        save(wnfilel) % save all workspace
    else
        disp('using existing zTDl')
    end
end
if strcmp(p.usem,'y') && strcmp(p.calavgm,'y') % mid cleaning
    %
    h4 = subplot(3,3,4);
    plot(h4,wsnew(iwindm(iXm)),mpwrtfm(if200m,inoisem(iXm)),'ro'); %use 500 Hz
    grid on; hold on;
    [isok{1,4}] = createFit4(wsnew(iwindm(iXm)),mpwrtfm(if200m,inoisem(iXm)),...
        wsnew(iwindm(iXm)),mpwrtfm(if200m,inoisem(iXm)),h4 );
    legend off
     ltex = sprintf('%s\n%s','Mid Band','200 Hz');
    ylabel(ltex);
    %
    h5 = subplot(3,3,5);
    plot(h5,wsnew(iwindm(iXm)),mpwrtfm(if500m,inoisem(iXm)),'ro'); %use 500 Hz
    ltex='500 Hz';
    grid on; hold on;
    [isok{1,5}] = createFit5(ltex,wsnew(iwindm(iXm)),mpwrtfm(if500m,inoisem(iXm)),...
        wsnew(iwindm(iXm)),mpwrtfm(if500m,inoisem(iXm)),h5 );
    legend off
    %
    h6 = subplot(3,3,6);
    plot(h6,wsnew(iwindm(iXm)),mpwrtfm(if1000m,inoisem(iXm)),'ro'); %use 500 Hz
    ltex='1000 Hz';
    grid on; hold on;
    [isok{1,6}] = createFit5(ltex,wsnew(iwindm(iXm)),mpwrtfm(if1000m,inoisem(iXm)),...
        wsnew(iwindm(iXm)),mpwrtfm(if1000m,inoisem(iXm)),h6 );
    legend off
    if ~exist('zTDm','var')
        zTDm = mintersect(isok{1,4}, isok{1,5}, isok{1,6});
        wnfilem = fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName2);
        save(wnfilem) % save all workspace
    else
        disp('using existing zTDm')
    end
end
if strcmp(p.use,'y') && strcmp(p.calavg,'y') % high cleaning
    h7 = subplot(3,3,7);
    plot(h7,wsnew(iwind(iX)),mpwrtf(if500,inoise(iX)),'ro'); %use 500 Hz
    grid on; hold on;
    [isok{1,7}] = createFit4(wsnew(iwind(iX)),mpwrtf(if500,inoise(iX)),...
        wsnew(iwind(iX)),mpwrtf(if500,inoise(iX)),h7 );
    legend off
     ltex = sprintf('%s\n%s','High Band','500 Hz');
    ylabel(ltex);
    %
    h8 = subplot(3,3,8);
    plot(h8,wsnew(iwind(iX)),mpwrtf(if1000,inoise(iX)),'ro'); %use 500 Hz
    ltex='1000 Hz';
    grid on; hold on;
    [isok{1,8}] = createFit5(ltex,wsnew(iwind(iX)),mpwrtf(if1000,inoise(iX)),...
        wsnew(iwind(iX)),mpwrtf(if1000,inoise(iX)),h8 );
    legend off
    %
    h9 = subplot(3,3,9);
    plot(h9,wsnew(iwind(iX)),mpwrtf(if2000,inoise(iX)),'ro'); %use 500 Hz
    ltex='2000 Hz';
    grid on; hold on;
    [isok{1,9}] = createFit5(ltex,wsnew(iwind(iX)),mpwrtf(if2000,inoise(iX)),...
        wsnew(iwind(iX)),mpwrtf(if2000,inoise(iX)),h9 );
    legend off
    if ~exist('zTD','var')
        zTD = mintersect(isok{1,7}, isok{1,8}, isok{1,9});
        wnfile = fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName1);
        save(wnfile) % save all workspace
    else
        disp('using existing zTD')
    end
end
%
% wnfignam2 = fullfile(p.harp.OutWinFig,[p.harp.dBaseName,'WindNoisefit']);
% figure(2);
% savefig(wnfignam2)
% disp('Done with cleaning');

%% Loop to allow editing
revise = 'd';
while strcmp(revise,'d')
    if strcmp(p.use,'y') % figure 4
        wnfile = fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName1);
        save(wnfile) % save all workspace
        wsfinal = wsnew(iwind(iX(zTD)));
        mpwfinal = mpwrtf(:,inoise(iX(zTD)));
        smpw = size(mpwfinal);
        nf = smpw(1);
        figure(4); clf;
        plot(wsfinal',mpwfinal(if500,:),'ko'); %use 1 kHz noise
        hold on
        plot(wsfinal',mpwfinal(if1000,:),'bo'); %use 1 kHz noise
        plot(wsfinal,mpwfinal(if10k,:),'ro'); %use 1 %use 10 kH
        plot(wsfinal,mpwfinal(if20k,:),'go'); %use 1 %use 20 kH
        legend('500 Hz','1 kHz','10 kHz','20 kHz','Location','northwest');
        title([p.harp.dBaseName,' High']);
        wnfignam4 = fullfile(p.harp.OutWinFig,[p.harp.dBaseName,'WindNoise']);
        figure(4);
        savefig(wnfignam4)
        %sort into speed bins
        [MPTF,lfor] = WindSort49(wsfinal,mpwfinal); % Mean Pressure TF corrected
        % Make TF Correction
        TFCorr = NaN(8,nf-1);
    end
    if strcmp(p.usem,'y') %figure 40
        wnfilem = fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName2);
        save(wnfilem) % save all workspace
        wsfinalm = wsnew(iwindm(iXm(zTDm)));
        mpwfinalm = mpwrtfm(:,inoisem(iXm(zTDm)));
        smpwm = size(mpwfinalm);
        nfm = smpwm(1);
        figure(40); clf;
        plot(wsfinalm',mpwfinalm(if50m,:),'ko'); %use 50 Hz noise
        hold on
        plot(wsfinalm',mpwfinalm(if100m,:),'bo'); %use 100Hz noise
         plot(wsfinalm',mpwfinalm(if200m,:),'ro'); %use 100Hz noise
          plot(wsfinalm',mpwfinalm(if500m,:),'go'); %use 100Hz noise
        legend('50 Hz','100 Hz','200 Hz','500 Hz','Location','southeast');
        title([p.harp.dBaseName,' Mid']);
        wnfignam40 = fullfile(p.harp.OutWinFig,[p.harp.dBaseName,'WindNoise_mid']);
        figure(40);
        savefig(wnfignam40)
        [MPTFm,lfor] = WindSort49(wsfinalm,mpwfinalm); % sort into speed
        %TF Correction
        TFCorrm = NaN(8,nfm-1); %
    end
    if strcmp(p.usel,'y') %figure 400
        wnfilel = fullfile(p.harp.OutFolder,'WindNoiseMat',p.harp.OutName3);
        save(wnfilel) % save all workspace
        wsfinall = wsnew(iwindl(iXl(zTDl)));
        mpwfinall = mpwrtfl(:,inoisel(iXl(zTDl)));
        smpwl = size(mpwfinall);
        nfl = smpwl(1);
        figure(400); clf;
        plot(wsfinall',mpwfinall(if20l,:),'ko'); %use 20 Hz noise
        hold on
        plot(wsfinall',mpwfinall(if50l,:),'bo'); %use 50Hz noise
        plot(wsfinall',mpwfinall(if100l,:),'ro'); %use 100Hz noise
        plot(wsfinall',mpwfinall(if200l,:),'go'); %use 200Hz noise
        legend('20 Hz','50 Hz','100 Hz','200 Hz','Location','southeast');
        title([p.harp.dBaseName,' Low']);
        wnfignam400 = fullfile(p.harp.OutWinFig,[p.harp.dBaseName,'WindNoise_mid']);
        figure(400);
        savefig(wnfignam400)
        [MPTFl,lfor] = WindSort49(wsfinall,mpwfinall); % sort into speed
        %TF Corr low
        TFCorrl10 = NaN(8,(nfl-1)/10); % reduce length by 10
        freql10 = freql(1:10:nfl); % reduce sampling by 10
        nfl10 = length(freql10);
    end
    %%
    if (~exist('depth') || isnan(depth))
        prompt = ' Please Enter Depth in m: ';
        depth = input(prompt);
    end
    %% Make TF Correction
    % Theory is fnm starts with .01 kHz (need to x 1000)
    % Data is freq and freqm starts with 0 then 100 (Hz)
    % fnm is frequency for the noise model
    % nf = number freq samples; nfm for mid; nfl for low
    itf = 1; % increases with wind speed
    % counts plot: mid-high, low-high, low only, mid only, high only
    iplotmh = 1; iplotlh = 1; iplotl = 1; iplotm = 1; iplot = 1; 
    % lfor is the highest wind force bin found by WindSort.m
    for  i = 2 : lfor    % start at i = 2 ss1 end i=9 ss8
        % Average MPTF: AMl AMm and AM start at 1Hz 10 Hz and 100 Hz
        if strcmp(p.use,'y') %
            smptf = size((MPTF{i}(2:nf,:)')); % 200 Hz to 30 kHz
            if smptf(1) > 1
                AM = mean(MPTF{i}(2:nf,:)');
            elseif smptf(1) == 1
                AM = MPTF{i}(2:nf,:)';
            end
            TFCorr(itf,1:nf-1) = kd(i,1:nf-1) - AM(1:nf-1); % first value is 100 Hz last highest freq
        end
        if strcmp(p.usem,'y') %
            smptfm = size((MPTFm{i}(2:nfm,:)')); %
            if smptfm(1) > 1
                AMm = mean(MPTFm{i}(2:nfm,:)'); % 20 Hz to 5 kHz
            elseif smptfm(1) == 1
                AMm = MPTFm{i}(2:nfm,:)'; %
            end
            TFCorrm(itf,1:nfm-1) = kdm(i,1:nfm-1) - AMm(1:nfm-1); % first value is 10 Hz last 1000 Hz
        end
        if strcmp(p.usel,'y') %
            smptfl = size((MPTFl{i}(2:nfl,:)')); %
            if smptfl(1) > 1
                AMl = mean(MPTFl{i}(2:nfl,:)'); % 2 Hz to 1 kHz
            elseif smptfl(1) == 1
                AMl = MPTFl{i}(2:nfl,:)'; %
            end
            % correct for 1 Hz sampling of AMl and 10 Hz sampling of kdm
            TFCorrl10(itf,:) = kdm(i,1:nfl10-1) - AMl(10:10:nfl-1); % first is 10 Hz last 1000 Hz
        end
        itf = itf + 1;

        %% Freq Plot
        if strcmp(p.FrePlt,'on')
            if strcmp(p.usem,'y') && strcmp(p.use,'y') % mid + high
                if iplotmh ==1
                    ssFigmh = figure(66); clf;
                    hold on
                end
                figure(ssFigmh);
                subplot(4,2,iplotmh)  % 8 subplots = 4 x 2
                semilogx(freq(5:nf),MPTF{i}(5:nf,:)); % from 200 Hz to 30 kHz
                hold on
                semilogx(freqm(3:41),MPTFm{i}(3:41,:)); % from 20 Hz to 200 Hz
                semilogx(1000*fnm,kdp(i,:),'r','LineWidth',3); % theory as a line
                semilogx(freq(5:nf),AM(4:nf-1),'k','Linestyle',':','LineWidth',3); % from 200 Hz to 30 kHz
                semilogx(freqm(3:41),AMm(2:40),'k','Linestyle',':','LineWidth',3); % from 20 Hz to 100 Hz
                grid on
                ftxt =['Beaufort Force' ,num2str(ss(i))];
                text(1000,94,ftxt)
                v = [20 (nf-1)*dfreq 25 110];
                xticks([10 100 1000 10000 100000])
                axis(v)
                if iplotmh > 6
                    xlabel('Frequency [Hz]')
                end
                if any(iplotmh == [1,3,5,7])
                    ylabel('dB re uPa^2/Hz')
                end
                if iplotmh == 1
                    title([p.harp.dBaseName,' Mid+High']);
                end
                iplotmh = iplotmh + 1;
            end
            %
            if strcmp(p.usel,'y') && strcmp(p.use,'y') % low + high
                % low freq plot
                if iplotlh ==1
                    ssFiglh = figure(606); clf;
                end
                figure(ssFiglh);
                subplot(4,2,iplotlh)  % 8 subplots = 4 x 2
                semilogx(freq(5:nf),MPTF{i}(5:nf,:)); % from 400 Hz to 100 kHz
                hold on
                semilogx(freql(2:401),MPTFl{i}(2:401,:)); % from 10 Hz to 400 Hz
                semilogx(1000*fnm,kdp(i,:),'r','LineWidth',3); % theory as a line
                semilogx(freq(5:nf),AM(4:nf-1),'k','Linestyle',':','LineWidth',3); % from 400 Hz to end
                semilogx(freql(2:401),AMl(1:400),'k','Linestyle',':','LineWidth',3); % from 10 Hz to 400 Hz
                grid on
                ftxt =['Beaufort Force' ,num2str(ss(i))];
                text(700,94,ftxt)
                v = [10 (nf-1)*dfreq 25 110];
                xticks([10 100 1000 10000 100000])
                axis(v)
                if iplotlh > 6
                    xlabel('Frequency [Hz]')
                end
                if any(iplotlh == [1,3,5,7])
                    ylabel('dB re uPa^2/Hz')
                end
                if iplotlh == 1
                    title([p.harp.dBaseName,' Low+High']);
                end
                iplotlh = iplotlh + 1;
            end
            %
            if strcmp(p.usel,'y') % low only freq plot
                if iplotl ==1
                    ssFigl = figure(600); clf;
                end
                figure(ssFigl);
                subplot(4,2,iplotl)  % 8 subplots = 4 x 2
                semilogx(freql(2:nfl),MPTFl{i}(2:nfl,:)); % from 10 Hz to nfl Hz
                hold on
                semilogx(1000*fnm,kdp(i,:),'r','LineWidth',3); % theory as a line
                semilogx(freql(2:nfl),AMl(1:nfl-1),'k','Linestyle',':','LineWidth',3); % from 10 Hz to 400 Hz
                grid on
                ftxt =['Beaufort Force' ,num2str(ss(i))];
                text(70,94,ftxt)
                v = [10 (nfl-1)*dfreql 50 110];
                xticks([10 100 1000])
                axis(v)
                if iplotl > 6
                    xlabel('Frequency [Hz]')
                end
                if any(iplotl == [1,3,5,7])
                    ylabel('dB re uPa^2/Hz')
                end
                if iplotl == 1
                    title([p.harp.dBaseName,' Low']);
                end
                iplotl = iplotl + 1;
            end
            %
            if strcmp(p.usem,'y') % mid only freq plot
                if iplotm ==1
                    ssFigm = figure(60); clf;
                end
                figure(ssFigm);
                subplot(4,2,iplotm)  % 8 subplots = 4 x 2
                semilogx(freqm(2:nfm),MPTFm{i}(2:nfm,:)); % from 10 Hz to nfm Hz
                hold on
                semilogx(1000*fnm,kdp(i,:),'r','LineWidth',3); % theory as a line
                semilogx(freqm(2:nfm),AMm(1:nfm-1),'k','Linestyle',':','LineWidth',3); % from 10 Hz to 400 Hz
                grid on
                ftxt =['Beaufort Force' ,num2str(ss(i))];
                text(70,94,ftxt)
                v = [10 (nfm-1)*dfreqm 50 110];
                xticks([10 100 1000])
                axis(v)
                if iplotm > 6
                    xlabel('Frequency [Hz]')
                end
                if any(iplotm == [1,3,5,7])
                    ylabel('dB re uPa^2/Hz')
                end
                if iplotm == 1
                    title([p.harp.dBaseName,' Mid']);
                end
                iplotm = iplotm + 1;
            end
            if strcmp(p.use,'y') % mid only freq plot
                if iplot ==1
                    ssFig = figure(6); clf;
                end
                figure(ssFig);
                subplot(4,2,iplot)  % 8 subplots = 4 x 2
                semilogx(freq(3:nf),MPTF{i}(3:nf,:)); % from 200 Hz to 30 kHz
                hold on
                semilogx(freq(2:nf),kd(i,1:nf-1),'r','LineWidth',3); % theory as a line
                semilogx(freq(3:nf),AM(2:nf-1),'k','Linestyle',':','LineWidth',3); % from 200 Hz to 30 kHz
                grid on
                ftxt =['Beaufort Force' ,num2str(ss(i))];
                text(1000,94,ftxt)
                v = [200 (nf-1)*dfreq 25 100];
                xticks([1000 10000 100000 ])
                axis(v)
                if iplot > 6
                    xlabel('Frequency [Hz]')
                end
                if any(iplot == [1,3,5,7])
                    ylabel('dB re uPa^2/Hz')
                end
                if iplot == 1
                    title([p.harp.dBaseName,' High']);
                end
                iplot = iplot + 1;
            end
        end
    end
    % Save SS Figures
    if strcmp(p.use,'y')
        ssfile = fullfile(p.harp.OutFolder,'VsFreqVsWind',p.harp.OutVsFreq);
        savefig(ssFig,ssfile) % high only
        if strcmp(p.usem,'y')
            ssfilemh = fullfile(p.harp.OutFolder,'VsFreqVsWind',...
                [p.harp.OutVsFreq(1:end-4),'mh.fig']);
            savefig(ssFigmh,ssfilemh) % mid + high
        end
        if strcmp(p.usel,'y')
            ssfilelh = fullfile(p.harp.OutFolder,'VsFreqVsWind',...
                [p.harp.OutVsFreq(1:end-4),'lh.fig']);
            savefig(ssFiglh,ssfilelh) % low + high
        end
    end
    if strcmp(p.usem,'y') %
        ssfilem = fullfile(p.harp.OutFolder,'VsFreqVsWind',p.harp.OutVsFreqm);
        savefig(ssFigm,ssfilem)
    end
    if strcmp(p.usel,'y') %
        ssfilel = fullfile(p.harp.OutFolder,'VsFreqVsWind',p.harp.OutVsFreql);
        savefig(ssFigl,ssfilel)
    end
    % REgression 
    if strcmp(p.use,'y')%  REgress high
        if (fs0 == 200000 || fs0 == 320000)% Frequencies for plotting
            ifr = [.2, .5, 1, 2, 5, 10, 20, 30, 40, 50, 75, 100]; % freq in kHz
        elseif (fs0 == 64000 || fs0 == 96000)
            ifr = [.2, .5, 1, 2, 5, 10, 20, 30, 40]; % freq in kHz
        elseif fs0 == 48000 || fs0 == 50000
            ifr = [.2, .5, 1, 2, 5, 10, 20]; % freq in kHz
        end
        [SlopeLR,OffSetLR,LRr2,SlopeTS,OffSetTS,SlopeTSu,OffSetTSu,fTS] = WNRegress1(...
            p.harp.OutFolder,p.harp.OutVsWind1,p.harp.Proj,p.harp.Site,p.harp.Depl,ttitle,depth,...
            wsfinal,mpwfinal,ptime,dnew,kd,ms,fs0,ifr,freq,...
            dfreq,p.RegPlt,8);
    end
    if strcmp(p.usem,'y') % REgress mid
        ifrm = [.02,.05,.1,.2,.3,.4,.5,.6,1];% Frequencies for plotting
        [SlopeLRm,OffSetLRm,LRr2m,SlopeTSm,OffSetTSm,SlopeTSum,OffSetTSum,fTSm] = WNRegress1(...
            p.harp.OutFolder,p.harp.OutVsWind2,p.harp.Proj,p.harp.Site,p.harp.Depl,ttitle,depth,...
            wsfinalm,mpwfinalm,ptimem,dnew,kdm,ms,fs0m,ifrm,freqm,...
            dfreqm,p.RegPlt,80);
    end
    if strcmp(p.usel,'y') % REgress low
        ifrl = [.005,.01,.02,.05,.1,.2,.3,.4,.5,.6,.8,1];% Frequencies for plotting
        [SlopeLRl,OffSetLRl,LRr2l,SlopeTSl,OffSetTSl,SlopeTSul,OffSetTSul,fTSl] = WNRegress1(...
            p.harp.OutFolder,p.harp.OutVsWind3,p.harp.Proj,p.harp.Site,p.harp.Depl,ttitle,depth,...
            wsfinall,mpwfinall,ptimel,dnew,kdm,ms,fs0l,ifrl,freql,...
            dfreql,p.RegPlt,800);
    end
    %% TF correction
    imaxSS = min(8,lfor) -1; % max SS to use in TFCorr, ifor or imaxSS =8 (ss7)
    if strcmp(p.usel,'y')
        % for low freq use only ss3 - ss7 (i = 4 - 8)
        MTFCorrl10 = mean(TFCorrl10(3:imaxSS,:),'omitnan');%starts with freq = 10 Hz
        figure(900); clf;
        semilogx(freql10(2:end),MTFCorrl10,'r-','LineWidth',3); %
        v = [10 (nfl-1) -10 10];
        axis(v)
        grid on
        title([p.harp.dBaseName,' Hydrophone ',tf_file(1:3)])
        xlabel('Frequency [Hz]')
        ylabel('dB re uPa//counts');
        % make new TF
        % Default TF CORRECTION > 300 Hz < 20 kHz
        col = floor(300/10);  % 300 Hz start
        coh = floor(nfl10-1-15);  % mod tf cutoff frequencies in Hz
        % Transfer function correction vector
        Ptf = interp1(p.tf.freq,p.tf.uppc,freql10,'linear','extrap');
        [TFoldl, TFnewl] = tfmakel49(10*dfreql,nfl10,col, coh, freql10, Ptf, MTFCorrl10);
        TFFigl = figure(1000); clf;
        [TFFigl] = tffigmakel(TFFigl,TFnewl,TFoldl,tf_file,p.harp.dBaseName,col,coh);
    end
    if strcmp(p.usem,'y')
        % for mid freq use only ss3 - ss7 (i = 4 - 8)
        MTFCorrm = mean(TFCorrm(3:imaxSS,:),'omitnan');%starts with freq = 10 Hz
        nMTm = isnan(MTFCorrm); % correct Nan
        inMTm = find(nMTm > 0);
        if (~isempty(inMTm))
            for i = 1 : length(inMTm)
                MTFCorrm(inMTm(i))=(MTFCorrm(inMTm(i)-1)+MTFCorrm(inMTm(i)+1))/2 ;
            end
        end
        figure(90); clf;
        semilogx(freqm(2:nfm),MTFCorrm(1:(nfm-1)),'r-','LineWidth',3); %
        v = [10 (nfm-1)*dfreqm -10 5];
        axis(v)
        grid on
        title([p.harp.dBaseName,' Hydrophone ',tf_file(1:3)])
        xlabel('Frequency [Hz]')
        ylabel('dB re uPa//counts');
        % make new TF
        % Default TF CORRECTION > 300 Hz < 20 kHz
        col = floor(300/10);  % 300 Hz start
        coh = floor((nfm-1-100));  % mod tf cutoff frequencies in Hz
%         freqmten = freqm(11:10:nfm);
        % Transfer function correction vector
        Ptf = interp1(p.tf.freq,p.tf.uppc,freqm,'linear','extrap');
        [TFoldm, TFnewm] = tfmakel49(dfreqm,nfm,col, coh, freqm, Ptf, MTFCorrm);
        TFFigm = figure(100); clf;
        [TFFigm] = tffigmakel(TFFigm,TFnewm,TFoldm,tf_file,p.harp.dBaseName,col,coh);
    end
    if strcmp(p.use,'y') %
        % for high freq use only ss1 - ss7
        MTFCorr = mean(TFCorr(2:imaxSS,:),'omitnan');  %starts with freq = 100 Hz
        nMT = isnan(MTFCorr); % correct for NaN
        inMT = find(nMT > 0);
        if (~isempty(inMT))
            for i = 1 : length(inMT)
                MTFCorr(inMT(i))=(MTFCorr(inMT(i)-1)+MTFCorr(inMT(i)+1))/2 ;
            end
        end
        MTFCorra = MTFCorr; %in case no mid
        figure(9); clf;
        semilogx(freq(2:nf),MTFCorr(1:nf-1),'r','LineWidth',3); %
        hold on
        if strcmp(p.usem,'y')
            MTFCorra = [MTFCorrm(10:10:40),MTFCorr(5:end)]; %replace point at 100 Hz
            semilogx(freq(2:nf),MTFCorra(1:nf-1),'k:','LineWidth',3); %
            semilogx(freqm(2:101),MTFCorrm(1:100),'r--','LineWidth',3); %
        end
        if strcmp(p.usel,'y')
            semilogx(freql10(2:end),MTFCorrl10,'r-.','LineWidth',3); %
            if ~strcmp(p.usem,'y')
                MTFCorra = [MTFCorrl10(1:4),MTFCorr(5:end)]; %replace point at 100 Hz
                semilogx(freq(2:nf),MTFCorra(1:nf-1),'k:','LineWidth',3); %
            end
        end
        v = [10 (nf-1)*dfreq -10 5];
        axis(v)
        grid on
        title([p.harp.dBaseName,' Hydrophone ',tf_file(1:3)])
        xlabel('Frequency [Hz]')
        ylabel('dB re uPa//counts')
        % Make New TF
        % Default TF CORRECTION > 300 Hz < 20 kHz
        col = floor(300/100);  % 300 Hz start
        coh = floor(20000/100);  % mod tf cutoff frequencies in Hz
        % Transfer function correction vector
        Ptf = interp1(p.tf.freq,p.tf.uppc,freq,'linear','extrap');
        [TFold, TFnew] = tfmake49(dfreq,nf,col, coh, freq, Ptf, MTFCorra,p.tf.tfn);
        %Make TF Figure
        TFFig = figure(10); clf;
        [TFFig] = tffigmake(TFFig,TFnew,TFold,tf_file,p.harp.dBaseName,col,coh);
    end
    
    %% Revise or Save data
    revise = input('Revise Data: d ; Cutoff: c; Bad x; Quit q; ','s');
    if strcmp(revise,'d')
        disp('Revise Data')
        celnums = inputdlg({'Enter Freq', 'Add=a Subtact=s','edit l m h'}, ...
            'Data Edit', [1 20; 1 20; 1 20]);
        efreq = str2double(celnums{1});
        addsub = (celnums{2});
        lmh = celnums{3};
        % make figure to edit
        figure(1); clf; set(1,'name',sprintf('Wind vs Noise')); h1 = gca;
        grid on; hold on;
        legend(h1,[num2str(efreq),' Hz'],'Location','southeast');
        if strcmp(lmh,'l') && strcmp(p.usel,'y')
            if efreq == 50
                ifx = if50l;
            elseif efreq == 100
                ifx = if100l;
            elseif efreq == 500
                ifx = if500l;
            elseif efreq == 1000
                ifx = if1000l;
            else
                ifx = if500l;
                disp([num2str(efreq),' not available for edit ... using 500'])
            end
            plot(h1,wsnew(iwindl),mpwrtfl(ifx,inoisel),'o'); %use low Hz noise
            [zTDl] = addsub49(addsub,zTDl);

        elseif strcmp(lmh,'m') && strcmp(p.usem,'y')
            if efreq == 50
                ifx = iffif;
            elseif efreq == 100
                ifx = ifhun;
            elseif efreq == 500
                ifx = if500;
            elseif efreq == 1000
                ifx = if1000;
            else
                ifx = if1000;
                disp([num2str(efreq),' not available for edit ... using 1000'])
            end
            plot(h1,wsnew(iwindm),mpwrtfm(ifx,inoisem),'o'); %use mid select
            [zTDm] = addsub49(addsub,zTDm);

        elseif strcmp(lmh,'h') && strcmp(p.use,'y')
            if efreq == 500
                ifx = if500;
            elseif efreq == 1000
                ifx = if1000;
            elseif efreq == 2000
                ifx = if2000;
            elseif efreq == 10000
                ifx = if10k;
            elseif efreq == 20000
                ifx = if20k;
            else
                ifx = if2000;
                disp([num2str(efreq),' not available for edit ... using 2000'])
            end
            plot(h1,wsnew(iwind),mpwrtf(ifx,inoise),'o'); %use selected noise
            [zTD] = addsub49(addsub,zTD);
        end
        %
    elseif strcmp(revise,'c')
        while (revise == 'c')
            cutlow = input('Low cutoff Hz: ');
            cuthigh = input('High cuttoff Hz: ');
            if strcmp(p.use,'y') %high freq
                col = floor(cutlow/100);
                coh = floor(cuthigh/100);
                Ptf = interp1(p.tf.freq,p.tf.uppc,freq,'linear','extrap');
                [TFold, TFnew] = tfmake49(dfreq,nf,col, coh, freq, Ptf, MTFCorra,p.tf.tfn);
                figure(TFFig); clf;
                [TFFig] = tffigmake(TFFig,TFnew,TFold,tf_file,p.harp.dBaseName,col,coh);
            end
            if  strcmp(p.usem,'y') % mid freq
                col = floor(cutlow/10);  % 300 Hz start
                coh = floor(cuthigh/10);  % mod tf cutoff frequencies in Hz
                if col > nfm-1-100
                    col = nfm-1-100;
                end
                if coh > nfm-1 -100
                    coh = nfm-1 - 100;
                end
                Ptf = interp1(p.tf.freq,p.tf.uppc,freqm,'linear','extrap');
                [TFoldm, TFnewm] = tfmakel49(dfreqm,nfm,col, coh, freqm, Ptf, MTFCorrm);
                figure(TFFigm); clf;
                [TFFigm] = tffigmakel(TFFigm,TFnewm,TFoldm,tf_file,p.harp.dBaseName,col,coh);
            end
            if  strcmp(p.usel,'y') % low freq
                col = floor(cutlow/10);  % 300 Hz start
                coh = floor(cuthigh/10);  % mod tf cutoff frequencies in Hz
                if col > nfl10-1
                    col = nfl10-1;
                end
                if coh > nfl10-1 -15
                    coh = nfl10-1 -15;
                end
                Ptf = interp1(p.tf.freq,p.tf.uppc,freql10,'linear','extrap');
                [TFold, TFnewl] = tfmakel49(10*dfreql,nfl10,col, coh, freql10, Ptf, MTFCorrl10);
                figure(TFFigl); clf;
                [TFFigl] = tffigmakel(TFFigl,TFnewl,TFold,tf_file,p.harp.dBaseName,col,coh);
            end
            revise = input('Cutoff: c; Quit q; ','s');
        end
    elseif strcmp(revise,'x') % bad result move files to "bad" folder
        SaveTF = 'no';
        wnfilebad = fullfile(p.harp.OutFolder,'BAD',p.harp.OutName1);
        status1 = movefile(wnfile,wnfilebad);
        wnfignam2bad = fullfile(p.harp.OutBad,[p.harp.dBaseName,'WindNoisefit']);
        status2 = movefile([wnfignam2,'.fig'],wnfignam2bad);
        wnfignam3bad = fullfile(p.harp.OutBad,[p.harp.dBaseName,'WindNoise']);
        status3 = movefile([wnfignam3,'.fig'],wnfignam3bad);
        ssfilebad = fullfile(p.harp.OutBad,p.harp.OutVsFreq);
        status4 = movefile(ssfile,ssfilebad);
        VsWindfile = fullfile(p.harp.OutFolder,'VsFreqVsWind',OutVsWind);
        VsWindfilebad = fullfile(p.harp.OutBad,OutVsWind);
        status5 = movefile(VsWindfile,VsWindfilebad);
        if (status1 && status2 && status3 && status4 && status5)
            disp('Successful Move to BAD folder')
        else
            disp(' Failed move to BAD folder')
        end
    elseif strcmp(revise,'q')
        % keep revise ~= 'd'
    else
        revise = 'd';
    end
end
% save in inverse sensitivity tf format in original TF Folder
if strcmp(p.SaveTF,'yes')
    tfn = str2num(tf_file(1:3));
     if strcmp(p.use,'y') %high fre
        tfcorrfile = fullfile(p.harp.OutFolder,'TFCorr',p.harp.OutTFCorr);
        freq = freq(2:nf); % remove 0 freq
        save(tfcorrfile,'freq','MTFCorra','MTFCorr','TFCorr','depth','tfn',...
            'SlopeLR','OffSetLR','LRr2','SlopeTS','OffSetTS','fTS','SlopeTSu','OffSetTSu',...
            'col','coh');
        % save TF_Wind in Folder in Output
        tfnewfile = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnew.tf']);
        save(tfnewfile,'TFnew','-ascii','-tabs');
        tfnewfig = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewfig']);
        savefig(TFFig,tfnewfig)
        tfnewfigpdf = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewfig.pdf']);
        saveas(TFFig,tfnewfigpdf)
    end
    if strcmp(p.usem,'y') %mid fre
        tfcorrfilem = fullfile(p.harp.OutFolder,'TFCorr',...
            [p.harp.OutTFCorr(1:end-4),'m.mat']);
         freqm = freqm(2:nfm);
        save(tfcorrfilem,'freqm','MTFCorrm','TFCorrm','depth','tfn',...
            'SlopeLRm','OffSetLRm','LRr2m','SlopeTSm','OffSetTSm','fTSm','SlopeTSum','OffSetTSum',...
            'col','coh');
        % save TF_Wind in Folder in Output
        tfnewfilem = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewm.tf']);
        save(tfnewfilem,'TFnewm','-ascii','-tabs');
        tfnewfigm = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewfigm']);
        savefig(TFFigm,tfnewfigm)
        tfnewfigmpdf = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewfigm.pdf']);
        saveas(TFFigm,tfnewfigmpdf)
    end
    if strcmp(p.usel,'y')
        tfcorrfilel = fullfile(p.harp.OutFolder,'TFCorr',...
            [p.harp.OutTFCorr(1:end-4),'l.mat']);
         freql = freql(2:nfl);
        save(tfcorrfilel,'freql','MTFCorrl10','TFCorrl10','depth','tfn',...
            'SlopeLRl','OffSetLRl','LRr2l','SlopeTSl','OffSetTSl','fTSl','SlopeTSul','OffSetTSul',...
            'col','coh');
        % save TF_Wind in Folder in Output
        tfnewfilel = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewl.tf']);
        save(tfnewfilel,'TFnewl','-ascii','-tabs');
        tfnewfigl = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewfigl']);
        savefig(TFFigl,tfnewfigl)
        tfnewfigpdfl = fullfile(p.harp.OutFolder,'TF_Wind',...
            [tf_file(1:3),'_',p.harp.Proj,p.harp.Site,p.harp.Depl,'_TFnewfigl.pdf']);
        saveas(TFFigl,tfnewfigpdfl)

    end
end
%