function [kdp,ss,fp,ms] = NoiseModelnew(depth,a,b,aof,bof,cof,nfacl,mfacl,mfacf,mfaca)
%Ambient Noise Curve
% includes mid and high frequency
%With Acknowldegement to Donald Ross and Gordon Wenz
%J Hildebrand 1-26-2008, 1-2020, 10-2020
% Kundsen8 a=2.8; b=100; aof=0; bof=100; cof=12; %new model 
% nfacl = 1000; 
nfacld = 400; % nfacl = nfac linear depth, nfacld non linear
nfacf(1) = 201; nfacf(2) = 101; nfacf(3) = 62; nfacf(4) = 41; %nfacf = freq for non-linear
% 
off = 50 + a*exp(-depth/b);
% frequency in kHz
fm = .01 : .01 : .1; % 10 Hz steps 10 Hz - 100 Hz
f = .1 : .1 : 160; % 100 Hz steps 100 Hz - 100 kHz
fp = [fm,f]; % frequency for plots
%
ms = [1, 2.5, 4.5, 6.7, 9.4, 12.3, 15.5, 19, 22.6, 26.5, 30.5];
lms = log10(ms);
ss = [.5, 1,2,3,4,5,6,7,8,9, 10];
beau = [1,2,3,4,5,6,7,8,9,10,11];
fWi = [100, 500, 1000, 5000, 10000, 20000, 30000, 100000, 160000];
nl = zeros(length(ms),length(f));
mfac = zeros(1,length(f));
% Noise Model Here
%%
% Model for Knudson curves 10 Hz - 100 Hz
n = 0.5; m= -1;  off = off -8.;
for iw = 1 % wind less than 1 m/s
    for iifm = 1:10 % 
        ofac = exp(-depth/bof)*aof*(cof-iifm/10)/cof; % depth and frequency
        nlm(iw,iifm) = ofac + off+6.4 + (n*20* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ;
        osavem(iw,iifm) = ofac + off+6.4;
        nsavem(iw,iifm) = n;
        msavem(iw,iifm) = m * 7 / 10;
    end
end
n = 1; 
for iw = 2 % wind 1 - 2.5 m/s
    for iifm = 1:10 % freq 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+1.5 + (n*20* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ;
        osavem(iw,iifm) = ofac + off+1.5 ;
        nsavem(iw,iifm) = n;
        msavem(iw,iifm) = m * 7 / 10;
    end
   
end
n = 1; 
for iw = 3 % wind 2.5 - 4.5 m/s
    for iifm = 1:10 % freq 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off-0.5 + (n*20* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; 
        osavem(iw,iifm) = ofac + off-0.5 ;
        nsavem(iw,iifm) = n;
        msavem(iw,iifm) = m * 7 / 10;
    end
   
end
n = 1.25;  off = off - 10;
for iw = 4 % wind 4.5 - 6.7 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+3.2 + (n *22* log10(ms(iw))) +4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+3.2 ;
        nsavem(iw,iifm) = n * 22 /20;
        msavem(iw,iifm) = m * 7 / 10;
    end
   
end
n=1.25;
for iw = 5 % wind 6.7 - 9.4 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+1.83 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+1.83 ;
        nsavem(iw,iifm) = n * 23 /20;
        msavem(iw,iifm) = m * 7 / 10;
    end
    
end
n=1.25; 
for iw = 6 % wind 9.4 - 12.3 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+2.4 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+2.4;
        nsavem(iw,iifm) = n* 23/20;
        msavem(iw,iifm) = m * 7 / 10;
    end
   
end
n=1.25; 
for iw = 7 % wind > than 12 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+2.18 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+2.18 ;
        nsavem(iw,iifm) = n* 23/20;
        msavem(iw,iifm) = m * 7 / 10;
    end
    
end
n=1.25; 
for iw = 8 % wind > than 12 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+2.27 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+2.27 ;
        nsavem(iw,iifm) = n* 23/20;
        msavem(iw,iifm) = m * 7 / 10;
    end
   
end
n=1.25; 
for iw = 9 % wind > than 12 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+2.54 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+2.54 ;
        nsavem(iw,iifm) = n * 23/20;
        msavem(iw,iifm) = m * 7 / 10;
    end
   
end
n=1.25; 
for iw = 10 % wind > than 12 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof- iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+2.72 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+2.72 ;
        nsavem(iw,iifm) = n * 23/20;
        msavem(iw,iifm) = m * 7 / 10;
    end
    
end
n=1.25; 
for iw = 11 % wind > than 12 m/s
    for iifm = 1:10 % 100-1000 Hz
        ofac = exp(-depth/bof)*aof*(cof-iifm/10)/cof;
        nlm(iw,iifm) = ofac + off+2.8 + (n *23* log10(ms(iw))) + 4 *m* log10(fm(iifm)) ; %m=1
        osavem(iw,iifm) = ofac + off+2.8 ;
        nsavem(iw,iifm) = n * 23/20;
        msavem(iw,iifm) = m * 7 / 10;
    end
end
%%
% Model for Knudson curves 100 Hz - 100 kHz
n = 0.5; m= 1;  off = off +18.;
for iw = 1 % wind less than 1 m/s
    for iif = 1:4 % freq 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof; % depth and frequency
        nl(iw,iif) = ofac + off+5.4 + (n*20* log10(ms(iw))) + 3 *(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off+5.4;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off+0.3 + (n*20* log10(ms(iw))) -  10.5*(m - mfac(iif))* log10(f(iif)) ; %9.5
        osave(iw,iif+iifm) = ofac + off+0.3;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 10.5 / 10;
    end
    for iif = 10:40 %  1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac+ off+0.1 + (n*20* log10(ms(iw))) - 10.4 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off+0.1;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 10.4 / 10;
    end
    for iif = 41:100 %  4100 Hz - 10000 Hz
        nl(iw,iif) = off+2 + (n*20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+2;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
    for iif = 101:length(f) % > 10000 Hz
        mfac(iif) = (1 - depth/mfacl)* mfaca * m * (mfacf - iif)/length(f);
        if iif > 400
            mfac(iif) = mfac(400);
        elseif iif < mfacf
            mfac(iif) = 0;
        end
        nl(iw,iif) = off+2 + (n*20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+2;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
end
n = 1; m= 1; 
for iw = 2 % wind 1 - 2.5 m/s
    for iif = 1:4 % freq 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-0.14 + (n*20* log10(ms(iw))) + 2.3 *(m - mfac(iif))* log10(f(iif)) ; 
        osave(iw,iif+iifm) = ofac + off-0.14 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-4.7 + (n*20* log10(ms(iw))) - 10.7*(m - mfac(iif))* log10(f(iif)) ;%11.5
        osave(iw,iif+iifm) = ofac + off - 4.7 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 10.7 / 10;
    end
    for iif = 10:40 %  1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-5 + (n*20* log10(ms(iw))) - 10.2 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-5 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 10.2 / 10;
    end
    for iif = 41:100 %  4100 Hz - 10000 Hz
        nl(iw,iif) = off-3 + (n*20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off-3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
    for iif = 101:length(f) % > 10000 Hz
        nl(iw,iif) = off-3 + (n*20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off-3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
end
n = 1; m= 1; 
for iw = 3 % wind 2.5 - 4.5 m/s
    for iif = 1:4 % freq 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-2.81 + (n*20* log10(ms(iw))) + 1.6 *(m - mfac(iif))* log10(f(iif)) ; 
                osave(iw,iif+iifm) = ofac + off-2.81 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-6.8 + (n*20* log10(ms(iw))) - 10*(m - mfac(iif))* log10(f(iif)) ;%9.3
        osave(iw,iif+iifm) = ofac + off-6.8 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 10 / 10;
    end
    for iif = 10:40 %  1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-7 + (n*20* log10(ms(iw))) - 8.7 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-7 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 8.7 / 10;
    end
    for iif = 41:100 %  4100 Hz - 10000 Hz
        nl(iw,iif) = off-4 + (n*20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off-4 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
    for iif = 101:length(f) % > 10000 Hz
        nl(iw,iif) = off-4. + (n*20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off-4 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
end
n = 1.25; m= 1; off = off - 10;
for iw = 4 % wind 4.5 - 6.7 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off+0.18 + (n *22* log10(ms(iw))) + 0.9 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off+0.18 ;
        nsave(iw,iif+iifm) = n * 22 /20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-0.1 + (n *20* log10(ms(iw))) -  7*(m - mfac(iif))* log10(f(iif)) ;%7.5
        osave(iw,iif+iifm) = ofac + off-0.1 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 7 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-0.4 + (n *20* log10(ms(iw))) - 10.2 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.4 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 10.2 / 10;
    end
    n = 1.0;
    for iif = 41:200 % 4100 Hz - 20000 Hz
        nl(iw,iif) = off+5.75 + (n *20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+5.75 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
    for iif = 201:length(f) % 20000 Hz - 100000 Hz
        nfac = (depth/nfacl+ 0.04*exp(-depth/nfacld))*0.3*n*(iif-nfacf(1))/length(f);
        nl(iw,iif) = off+5.7 + ((n - nfac) *20* log10(ms(iw))) - 13.5 *(m - mfac(iif))* log10(f(iif));%m = 1.6
        osave(iw,iif+iifm) = off+5.7 ;
        nsave(iw,iif+iifm) = n-nfac;
        msave(iw,iif+iifm) = - m * 13.5 / 10;
    end
end
n=1.25; m= 1; 
for iw = 5 % wind 6.7 - 9.4 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-1.9 + (n *23* log10(ms(iw))) + 0.2 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-1.9 ;
        nsave(iw,iif+iifm) = n * 23 /20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off+0 + (n *20* log10(ms(iw))) - 6*(m - mfac(iif))* log10(f(iif)) ;%5
        osave(iw,iif+iifm) = ofac + off ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 6 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-0.3 + (n *20* log10(ms(iw))) - 12.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 12.5 / 10;
    end
    n = 1.0;
    for iif = 41:200 % 4100 Hz - 20000 Hz
        nl(iw,iif) = off+5.5 + (n *20* log10(ms(iw))) - 14 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+5.5 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14 / 10;
    end
    for iif = 201:length(f) % 20000 Hz - 100000 Hz
        nfac = (depth/nfacl+ 0.08*exp(-depth/nfacld))*0.3*n*(iif-nfacf(1))/length(f);
        nl(iw,iif) = off+5.5 + ((n - nfac) *20* log10(ms(iw))) - 14 *(m - mfac(iif))* log10(f(iif));%m = 1.6
        osave(iw,iif+iifm) = off+5.5 ;
        nsave(iw,iif+iifm) = n-nfac;
        msave(iw,iif+iifm) = - m * 14 / 10;
    end
end
n=1.25; m= 1; 
for iw = 6 % wind 9.4 - 12.3 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-2.10 + (n *23* log10(ms(iw))) - 0.5 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-2.10;
        nsave(iw,iif+iifm) = n* 23/20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-0.12 + (n *20* log10(ms(iw))) - 6*(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off-0.12 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 6 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-0.3 + (n *20* log10(ms(iw))) - 13.7 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 13.7 / 10;
    end
    n = 1.0;
    for iif = 41:100 % 4100 Hz - 10000 Hz
        nl(iw,iif) = off+5.6 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) =  off+5.6 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    for iif = 101:length(f) % 10000 Hz - 100000 Hz
        nfac = (depth/nfacl+ 0.17*exp(-depth/nfacld))*0.7*n*(iif-nfacf(2))/length(f);
        nl(iw,iif) = off+5.65 + ((n - nfac) *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)); %m = 1.6
        osave(iw,iif+iifm) = off+5.65 ;
        nsave(iw,iif+iifm) = n-nfac;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
        %  no wind dependence > 10 kHz for > 15 m/s
    end
end
n=1.25; m= 1; 
for iw = 7 % wind > than 12 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-3.02 + (n *23* log10(ms(iw))) - 1.2 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-3.02 ;
        nsave(iw,iif+iifm) = n* 23/20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-.3 + (n *20* log10(ms(iw))) - 5*(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 5 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-.3 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    n = 1.0;
    for iif = 41:61 % 4100 Hz - 6000 Hz
        nl(iw,iif) = off+5.75 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+5.75 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    for iif = 62:length(f) % 6000 Hz - 100000 Hz
        nfac(iif) = (depth/nfacl+ 0.24*exp(-depth/nfacld))*1.2*n*(iif-nfacf(3))/(length(f));
        nl(iw,iif) = off+5.85 + ((n - nfac(iif)) *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+5.85 ;
        nsave(iw,iif+iifm) = n-nfac(iif);
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
end
n=1.25; m= 1; 
for iw = 8 % wind > than 12 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-3.54 + (n *23* log10(ms(iw))) - 1.9 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-3.54 ;
        nsave(iw,iif+iifm) = n* 23/20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-.3 + (n *20* log10(ms(iw))) - 6*(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 6 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-.3 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    n = 1.0;
    for iif = 41:length(f) % 4100 Hz - 100000 Hz
        nfac(iif) = (0.6*exp(-depth/nfacld) + depth/nfacl)*1.8*n*(iif-nfacf(4))/(length(f));
        nl(iw,iif) = off+6.1 + ((n - nfac(iif)) *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) =  off+6.1 ;
        nsave(iw,iif+iifm) = n - nfac(iif);
        msave(iw,iif+iifm) = - m * 14.5 / 10;
        %  no wind dependence > 10 kHz for > 15 m/s
    end
end
n=1.25; m= 1; 
for iw = 9 % wind > than 12 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-4.06 + (n *23* log10(ms(iw))) - 2.6 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-4.06 ;
        nsave(iw,iif+iifm) = n * 23/20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-.3 + (n *20* log10(ms(iw))) - 6*(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 6 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-.3 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    n= 1.0;
    for iif = 41:length(f) % 4100 Hz - 100000 Hz
        nfac(iif) = (depth/nfacl+ 1*exp(-depth/nfacld))*2.2*n*(iif-nfacf(4))/(length(f));
        nl(iw,iif) = off+6.4 + ((n - nfac(iif)) *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+6.4 ;
        nsave(iw,iif+iifm) = n - nfac(iif);
        msave(iw,iif+iifm) = - m * 14.5 / 10;
        %  no wind dependence > 10 kHz for > 15 m/s
    end
end
n=1.25; m= 1; 
for iw = 10 % wind > than 12 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-4.58 + (n *23* log10(ms(iw))) - 3.3 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-4.58 ;
        nsave(iw,iif+iifm) = n * 23/20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-.3 + (n *20* log10(ms(iw))) - 6*(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 6 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-.3 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    n=1.0;
    for iif = 41:length(f) % 4100 Hz - 100000 Hz
        nfac(iif) = (depth/nfacl+ 1.6*exp(-depth/nfacld))*2.5*n*(iif-nfacf(4))/(length(f));
        nl(iw,iif) = off+6.85 + ((n - nfac(iif)) *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+6.85 ;
        nsave(iw,iif+iifm) = n - nfac(iif);
        msave(iw,iif+iifm) = - m * 14.5 / 10;
        %  no wind dependence > 10 kHz for > 15 m/s
    end
end
n=1.25; m= 1; 
for iw = 11 % wind > than 12 m/s
    for iif = 1:4 % 100-400 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-5.2 + (n *23* log10(ms(iw)))  - 4 *(m - mfac(iif))* log10(f(iif)) ; %m=1
        osave(iw,iif+iifm) = ofac + off-5.2 ;
        nsave(iw,iif+iifm) = n * 23/20;
        msave(iw,iif+iifm) = m * 7 / 10;
    end
    for iif = 5:9 % 500-900 Hz
        ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        nl(iw,iif) = ofac + off-.3 + (n *20* log10(ms(iw))) - 6*(m - mfac(iif))* log10(f(iif)) ;
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 6 / 10;
    end
    for iif = 10:40 % 1000 Hz - 4000 Hz
        if iif <= cof
            ofac = exp(-depth/bof)*aof*(cof-iif)/cof;
        end
        nl(iw,iif) = ofac +off-.3 + (n *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = ofac + off-0.3 ;
        nsave(iw,iif+iifm) = n;
        msave(iw,iif+iifm) = - m * 14.5 / 10;
    end
    n = 1.0;
    for iif = 41:length(f) % 4100 Hz - 100000 Hz
        nfac(iif) = (depth/nfacl+ 1.8*exp(-depth/nfacld))*2.7*n*(iif-nfacf(4))/(length(f));
        nl(iw,iif) = off+7.15 + ((n - nfac(iif)) *20* log10(ms(iw))) - 14.5 *(m - mfac(iif))* log10(f(iif)) ; %m = 1.6
        osave(iw,iif+iifm) = off+7.15 ;
        nsave(iw,iif+iifm) = n - nfac(iif);
        msave(iw,iif+iifm) = - m * 14.5 / 10;
        %  no wind dependence > 10 kHz for > 15 m/s
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% depth dependance correction to Knudsen spectra
% numerically solve eqn (3)->(4) from Kurahshi and Gratta 2007
% or similarly eqn (10)->(12) Short 2005 IEEE
% depth of hydrophone:
% assume dipole sources so n=2 in cos^(n-1)(theta)
h = depth;   % [meters
% start,end step angle [rad]
ti = 0;
to = pi/2;
dt = to/90;
% alpha -> sound absorbtion coefficient
T = 10; %           T is temperature in deg-C  
S = 35; %           S is salinity in PSU      
pH = 8; %           pH (default is 8)
[alpha] = AMCAtten(fp,h,T,S,pH); % Ainslie and McColm
alphah = alpha * h/1000;
% e^-ah
eah = 10.^(-alphah/10);
% loop over angle (theta)
Joah = 0;
Joo= 0;
for t = ti:dt:to
    ct = cos(t);
    sct = sec(t);
    st = sin(t);
    po = ct * st;
    pc = po * eah.^sct;
    Joah = Joah + pc;
    Joo = Joo + po;
end
% depth dependance correction to Knudsen spectra
ddc = 10 .* log10(Joah ./ Joo);
kdp = [nlm,nl];% theoretical noise level with ss
kdp = kdp + ones(length(ms),1)*ddc; % depth dependent correction
% Thermal Noise curve
nt = zeros(1,length(fp));
for iif = 1:length(fp)
    nt(iif) = -15 + 20 * log10(fp(iif));
end
% %combine kd and nt
for j = 1 : length(kdp(:,1))
    for i = 1 : length(nt)
        if kdp(j,i) < nt(i)
            kdp(j,i) = nt(i);
        end
    end
end
%
end