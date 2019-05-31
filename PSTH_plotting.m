%% Script for computing FFT for STN and GPe

% Arguments
%Spike times (linear_S=[times,number ID])

% Output
%FFT plots

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%

% Number of neurons
%STN -> Mandali2015
nSTN=32; % (nSTNxnSTN network size)
Mstn=nSTN;
Nstn=nSTN;
Pstn=Mstn*Nstn;

%GPe -> Mandali2015
nGPe=32; % (nGPexnGPe network size)
Mgpe=nGPe;
Ngpe=nGPe;
Pgpe=Mgpe*Ngpe;

%SNc -> Cullen2015
nSNc=8; % (nSNcxnSNc network size)
Msnc=nSNc;
Nsnc=nSNc;
Psnc=Msnc*Nsnc;

% Time parameters
dt=0.1; % Timestep
tspan=0:dt:dur;
Ttime=numel(tspan);

sec=0.001;
binsize=100; % msec
fs=10000;
temp1=ST_psth(stn_firings);
temp2=ST_psth(gpe_firings);
temp3=ST_psth(snc_firings);
% Trade-off b/w binsize and desired results
[ph1,ed1,r1]=psth(temp1,binsize,fs,Pstn,Ttime,0); %225-binsize
[ph2,ed2,r2]=psth(temp2,binsize,fs,Pgpe,Ttime,0);
[ph3,ed3,r3]=psth(temp3,binsize,fs,Psnc*2,Ttime,0);

fig81=figure(81);
set(fig81, 'Position', [5, 50, 1920, 955]);
subplot(311)
plot(ed1,r1,'b');
% xlim([0 numel(ed1)]);
title('STN')
refline([0 mean(r1)]);
I11=strcat(num2str(mean(r1)));legend(I11);
subplot(312)
plot(ed2,r2,'r')
% xlim([0 numel(ed2)]);
title('GPe')
refline([0 mean(r2)]);
I11=strcat(num2str(mean(r2)));legend(I11);
subplot(313)
plot(ed3,r3,'k')
% xlim([0 numel(ed2)]);
title('SNc')
refline([0 mean(r3)]);
I11=strcat(num2str(mean(r3)));legend(I11);
f811=strcat('PSTH_',filename0);
saveas(fig81,f811,'png');