function [wlatsnc]= weightcal_snc(DA)

%% Computing lateral weight matrix of SNc

% Arguments
%DA: Dopamine signal

% Output
%wlatsnc: Dopamine-dependent lateral connection weights for SNc

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
%----------------SNc_lateral connections-----------%
smax = 0.1;  %Strength of lateral connections in snc
rs = 1.6; %Radius of lateral connections in snc
nlatsnc = 15; % number of laterals in snc

ssmax = (smax*(1+35*DA))/(1-(0.1*DA*DA));
wlatsnc = calclatwts(nlatsnc,ssmax, rs);

end