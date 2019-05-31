function [Rvalue,Ravg]=mrcalculate(linear_S,Nneur,Ntime)

%% Computing synchrony value across time

% Arguments
%linear_S: Spike times (linear_S=[times,number ID])
%Nneur: Number of neurons
%Ntime: Simulation time

% Output
%Rvalue: Synchrony value across time
%Ravg: Average synchrony value

% References
%Pinsky PF, Rinzel J (1995) Synchrony measures for biological neural networks. Biol Cybern 73:129–137.

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
Rvalue=[];phi=[];

phi=3000*ones(Nneur,Ntime-1);
for neur=1:Nneur
    temptime=linear_S((linear_S(:,2)==neur));
    % temptime =[4    12    21    30    60    78   100   117   126   163   503   652   797   857   940   943];
    j=1;
    
    while j<numel(temptime)
        for i=temptime(j):1:temptime(j+1)-1
            phi(neur,i)=(2*pi*(i-temptime(j)))/(temptime(j+1)-temptime(j));
        end
        j=j+1;
    end
end
a=sqrt(-1);
tempM=sum(phi)/numel(phi);
M=exp(a*tempM);
Rvalue=((sum(exp(a*phi))/neur))./M;
Ravg=sum(abs(Rvalue))/numel(Rvalue);

end