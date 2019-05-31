function [output]=sub_sampling(input,dt)

%% Subsampling of time series data

% Arguments
%input: time series data

% Output
%output: Subsampled input data

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
[mm,nn]=size(input);
Ntime=(nn)*dt;

output=zeros(mm,Ntime);
for tr=1:mm
    start=1;stop=1/dt;
    for te=1:Ntime
        output(tr,te)=mean(input(tr,start:stop));
        start=start+1/dt;stop=stop+1/dt;
    end
end

end