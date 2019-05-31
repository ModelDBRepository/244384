function [Rvalue,Ravg]=mrintercalculate(linear_Stn,linear_Gpe,Nneur,Ttime)

%% Computing synchrony value between STN and GPe across time

% Arguments
%linear_Stn: STN Spike times (linear_S=[times,number ID])
%linear_Gpe: GPe Spike times (linear_S=[times,number ID])
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
Rvalue=[];phistn=zeros(Nneur,Ttime);
phigpe=zeros(Nneur,Ttime);
for neur=1:Nneur
    temptime1=linear_Stn((linear_Stn(:,2)==neur));
    temptime2=linear_Gpe((linear_Gpe(:,2)==neur));
    j=1;k=1;
    for i=1:Ttime % initial transients
        if j<numel(temptime1)
            if i>=temptime1(j) && i<=temptime1(j+1)
                phistn(neur,i)=((2*pi*(i-temptime1(j)))/(temptime1(j+1)-temptime1(j)));
                if i==temptime1(j+1)
                    j=j+1;
                end
            end
        end
        
        if k<numel(temptime2)
            if i>=temptime2(k) && i<=temptime2(k+1)
                phigpe(neur,i)=((2*pi*(i-temptime2(k)))/(temptime2(k+1)-temptime2(k)));
                if i==temptime2(k+1)
                    k=k+1;
                end
            end
        end
        
    end
end

phistn1=phistn(:,:);
phigpe1=phigpe(:,:);

a=sqrt(-1);
tempMstn=sum(phistn1)/numel(phistn1);
tempMgpe=sum(phigpe1)/numel(phigpe1);
M=exp(a*((tempMstn+tempMgpe)/2));
M=M';

sumphistn=(sum(exp(a*phistn1))/neur);
sumphigpe=(sum(exp(a*phigpe1))/neur);
sumphistn=sumphistn';
sumphigpe=sumphigpe';
Rvalue=((sumphistn+sumphigpe)/2)./M;
Ravg=sum(abs(Rvalue))/numel(Rvalue);

end