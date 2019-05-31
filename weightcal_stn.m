function [wlatstn]= weightcal_stn(DA)

%% Computing lateral weight matrix of STN

% Arguments
%DA: Dopamine signal

% Output
%wlatgpe: Dopamine-dependent lateral connection weights for STN

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
%----------------STN_lateral connections-----------%
smax = 1.3;  %Strength of lateral connections in stn
rs = 1.4; %Radius of lateral connections in stn
nlatstn = 11; % number of laterals in stn

ssmax = smax.*exp(-4.87.*DA);
wlatstn = calclatwts(nlatstn,ssmax, rs);

end