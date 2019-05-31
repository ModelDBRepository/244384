%% Spiking Excitotoxicity Model (SEM)
%Script for running the simulation

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%References
%Muddapu VR, Mandali A, Chakravarthy VS, Ramaswamy S (2019) A computational model of loss of dopaminergic cells in Parkinson’s disease
%due to glutamate-induced excitotoxicity. Front Neural Circuits 13:11
%Available at: https://www.frontiersin.org/articles/10.3389/fncir.2019.00011/abstract [Accessed February 25, 2019].

%%
clc;clear;close all;
tic
time=clock;curdate=time(3);curmonth=time(2); % Tracking time

dur=50000; % Duration of simulation in milliseconds

randinit=1; % 0-No random initializatiion, 1-Random initializatiion
nolat=1; % 0-Lateral connections off, 1-Lateral connections on
pd=1; % 0-PD condition on, 1-PD condition off
lstsn=1; % 0-STN-SNc connection off, 1-STN-SNc connection on

wsg=1; % Connection strength from STN to GPe
wgs=20; % Connection strength from GPe to STN

sthrsnc=10.6; % Stress thershold
sthrsnc1=deci2str(sthrsnc);

wstsn=1; % Connection strength from STN to SNc
wstsn1=deci2str(wstsn);

GPUon=0; % 0-GPU usage off, 1-GPU usage on

filename0=strcat('SEM_Wstsn',num2str(wstsn1),'_Sthr',num2str(sthrsnc1),'_',num2str(dur),'msec_',num2str(curdate),'_',num2str(curmonth));

if(GPUon==1)
    [snc_firings,stn_firings,gpe_firings,DA,srnd,simtime]=SEM_GPU(dur,filename0,randinit,nolat,pd,lstsn,wsg,wgs,sthrsnc,wstsn);
else
    [snc_firings,stn_firings,gpe_firings,DA,srnd,simtime]=SEM(dur,filename0,randinit,nolat,pd,lstsn,wsg,wgs,sthrsnc,wstsn);
end

toc