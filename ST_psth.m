function [temp1]=ST_psth(matrix)

%% Extracting spike times from neuronal firing data

% Arguments
% matrix = firings (column1->Time point; column2->Neuron ID)

% Output
% temp1 = spike times for each neuron (row->neuron ID; column->timepoints)

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
neur=unique(matrix(:,2));
temp1=[];
for i=1:numel(neur)
    temp=matrix((matrix(:,2)==neur(i)));
    temp1=[temp1;temp];
end

end