function [wlatgpe]= weightcal_gpe(DA)

%% Computing lateral weight matrix of GPe

% Arguments
%DA: Dopamine signal

% Output
%wlatgpe: Dopamine-dependent lateral connection weights for GPe

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
%----------------GPe_lateral connections-----------%
smax = 0.1002;  %Strength of lateral connections in gpe
rs = 1.6; %Radius of lateral connections in gpe
nlatgpe = 15; % number of laterals in gpe

ssmax = (smax.*(1+43.*DA))./(1-(0.1.*DA.*DA));
wlatgpe = calclatwts(nlatgpe,ssmax, rs);

end