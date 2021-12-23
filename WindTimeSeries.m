%Jah Extract Wind time series from CCMP data
% Oct 2021 added ERA5 option
% Oct 8, 2019
%
clear;
wmod = 'CCMP';
pre = 'H:\Wind_Model\CCMP\';
mid = '\CCMP_Wind_Analysis_';
post = '_V02.0_L3.0_RSS.nc';
%
% wmod = 'ERA5';
% pre = 'H:\Wind_Model\ERA5\';
% mid = '\';
% post = '.nc';
for yr = 2016 : 2018
    year = num2str(yr);
    % deployment site
    % site
    site = 'NRS16';
    outdir = 'H:\Wind_Data\GofMX';  % Output file directory
    disp([site,' ',year]);
    % average of all deployment locations for each site
    %G of AK
    %      slat = 58.66302424; slon =	148.0451015;  %CB
    %  slat = 56.24340417;	slon = 142.7572458;  %PT
    %      slat = 56.34023333;	slon = 145.1855944; %QN
    %      slat = 57.51366667;	slon = 146.5008333; %AB
    %      slat = 59.00897667; slon =	148.90422; %CA
    %      slat = 57.33508889;	slon = 150.6878222; %KO
    %      slat = 52.076000000; slon = 184.3601667; %BD
%     slat = 57-13.440;	slon = 150-31.700; %KOA

    %     %SOCAL
    %      slat = 32.89312;	slon = 117.3790167; %P
    %        slat = 33.25184697;	slon = 118.2506545; %A
    %       slat = 33.22731667;	slon = 118.2755229; %A2
    %       slat = 34.274958; slon = 	120.0244633; %B
    %      slat = 34.31946708; slon =	120.8053204; %C
    %       slat = 33.51263095; slon =	119.2500944; %M
    %       slat = 32.36977843; slon =	118.5645078; %N
    %       slat = 32.65695833; slon =	119.4771906; %E
    %       slat = 32.84412432; slon =	119.8441243; %H
    %       slat = 33.820175; slon = 	118.6290333; %Q
    % %      slat = 33.16013889; slon =	120.0092806; %R
    %      slat = 32.484875; slon =	118.2724111; %S
    %      slat = 32.91499;    slon = 120.37539; %SN
    %       slat = 32.88024; slon =	117.5742; %T
    %       slat = 32.9263; slon = 118.6170; %G
    %      slat = 33.47808333; slon = 122.5496667;  % CCE!
    %     slat = 35.400000000; slon=  121.562500000; % DCPPC
    %      slat = 31.747000000;  slon = 121.378000000; %CORC
    % slat = 36.4351;	slon = 121.9608; % Granite Canyon
    % slat =31.85163333;	slon =	118.4845167; % U
    % slat =29.14103333	;	slon =118.2609667; %BajaGI
    % slat = 33+(03.023/60); slon =	118+(40.893/60); %FLIP 08N1

    % GOFMX
    %   slat = 29.0490; slon = 86.0975; %GofMx DC
    %   slat = 25.53484394; slon = 84.6340803; %GOM DT
    %   slat = 27.55671212; slon = 91.1689303; %GOM GC
    %   slat = 28.84642381;  slon = 88.46614286; %GOM MC1
    %   slat = 28.97951;  slon = 88.46810667; %GOM MC2
    %   slat = 29.25399697; slon = 88.29345455; % GOM MP
    %   slat = 25.02466; slon =	84.39344333; % HH
    % slat = 27.7330833; slon = 92.46221806; %GOM EF
    % slat =27.88453333;	slon =92.40936667;%GOM EI
    % slat =	27.65405;	slon =	93.39408333;% GOM WF
    % slat = 28.6292;	slon = 	90.04048333; %GOM GI
    % slat =29.28108333;	slon =	87.85831667; %GOM EP
   slat =  28.00;  slon = 86.994; %GOM NRS06 2014
    slat = 28.2502; slon = 86.8327; %GOM NRS06 2016
 
    % %WAT
    % slat = 30.75156667 ; slon =	 79.8587375; %JAXD
    % slat = 37.16671667;  slon =	27.996; %NFCA
    % slat = 35.58413333; slon = 	74.74985; %HATB
    % slat = 39.83248333;	slon= 69.98205; %WC
    % slat = 40.26331667;	slon=	67.98626667;% OC
    % slat =	40+ (01.967/60); slon = 67+ (59.301/60); %BR
    % slat = 39.190775; slon = 72.22791667; % BC
    % slat = 39.83248333; slon =  69.98205; % NC
    % slat = 30.58340833; slon = 77.390575; % BS
    % slat = 33.666325; slon = 76.00875833;  % GS
    % slat = 32.10649167;  slon =	77.09220833;% BP

    %ADRIA
    % slat = 42+(20.651/60); slon = 360-(17+(49.235/60)); %JD east long

    %ANTAR
    %ARCTIC
    % slat = 72.72497; slon = 74.5; % Pond Inlet
    % slat = 74.34078; slon = 94.53752; % Barrow Strait

    %HAWAII
    % slat = 17.96308; slon = 145.481117; % Pagan in East Lon
    % slat =	27+ (43.774/60); slon = 175+ (33.221/60); %PHA12

    %OCNMS
    %OOI
    % slat = 45.81681; slon =	129.75421; % AXBS
    % slat = 44.51519; slon =	125.38987; % AXBS
    %SanctSound
    % slat = 36.798; slon =	121.976; % MB01
    % slat = 48.3938; slon =	124.654; % OC01

    % convert to East Lon
    slon = 360 - slon;
    %pre = '/Volumes/CzikoAudio1/WIND/CCMP/';
    %pre = '/Users/jah 1 2/Ddirectory/Whale/WIND/CCMP/';

    % get latitude and longitude loction from first file

    las = [mid,year,'01','01',post];

    f = fullfile(pre,year,las);
    lat = ncread(f,'latitude');
    lon = ncread(f,'longitude');
    % find point near deployment site

    slat25 = (round((slat+.125) * 4))/4 - .125;
    slon25 = (round((slon+.125) * 4))/4 - .125;

    % Difference with wind site?
    disp([' Del-Lat = ',num2str(slat - slat25),' deg']);
    disp([' Del-Lon = ',num2str(slon - slon25),' deg']);
    ilat = find(lat == slat25);
    ilon = find(lon == slon25);
    dnvec = zeros(4*365,1);  ivec = 0;
    uwvec = dnvec;  vwvec = dnvec;
    for month = 1:12
        if month < 10
            ml = num2str(month);
            m = ['0',ml];
        else
            m = num2str(month);
        end
        for day = 1:31
            if day < 10
                dl = num2str(day);
                d = ['0',dl];
            else
                d = num2str(day);
            end
            las = [mid,year,m,d,post];
            f = fullfile(pre,year,las);
            isfile = exist(f,'file');
            if isfile == 2
                %n = ncinfo(f); % read data for one day
                tihr = ncread(f,'time');
                uwnd = ncread(f,'uwnd');
                vwnd = ncread(f,'vwnd');
                for i = 1:4
                    ivec = ivec + 1;
                    dnvec(ivec) = datenum([1987,0,1,tihr(i),0,0]);
                    %test date
                    if i == 1
                        dntest = datenum([yr,month,day,0,0,0]);
                        if dntest ~= dnvec(ivec)
                            disp(['Date Prob ',dntest, dnvec(ivec)]);
                        end
                    end
                    uwvec(ivec) = uwnd(ilon,ilat,i);
                    vwvec(ivec) = vwnd(ilon,ilat,i);
                end
            else
                disp(['Not Valid Date:  ',year,m,d]);
            end
        end
    end
    % plot wind speed
    figure
    subplot(2,1,1)
    wspeed = sqrt(uwvec.^2 + vwvec.^2);
    % take out bad data
    isok = find(wspeed > 0);
    wspeed = wspeed(isok);
    dnvec = dnvec(isok);
    uwvec = uwvec(isok);
    vwvec = vwvec(isok);
    %
    plot(dnvec,wspeed);
    datetick('x','mmm')
    title([site,' Wind Speed ',num2str(year),'  lat= ',num2str(slat25),...
        '  lon= ',num2str(slon25)]);
    ylabel('Wind Speed in m/s');
    xlabel(' Month ');
    %plot wind direction
    subplot(2,1,2)
    wdir = 180*atan2(vwvec,uwvec)/pi; % v is north u is east
    wmetdir = wdir;
    idir = find (wdir < 0); % meterological convention 0-360 deg
    wmetdir(idir) = 360 + wdir(idir);
    plot(dnvec,wmetdir);
    datetick('x','mmm')
    title([site,' Wind Direction ',num2str(year),'  lat= ',num2str(slat25),...
        '  lon= ',num2str(slon25)]);
    ylabel('Meterological Wind Direction re: North');
    xlabel(' Month ');
    %save results
    fs = fullfile(outdir,year,[site,'windvec']);
    save(fs,'site','year','slat25','slon25','dnvec','uwvec','vwvec',...
        'wspeed','wmetdir');
    saveas(gcf,fs,'pdf');
    %
end