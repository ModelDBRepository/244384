function ff=deci2str(x)

%% Converting decimal number into string

% Arguments
%x: input (decimal number)

% Output
%ff: string

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
f=strcat(num2str(x));
ff=strrep(f,'.','-');

end