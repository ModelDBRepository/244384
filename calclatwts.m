function w = calclatwts(nlat, sig, rad)

%% Computing lateral weight matrix

% Arguments
%nlat is the size of the window. It must be an odd number.
%sig: strength of the connections
%rad: radius of neighborhood

% Output
%w: weight matrix with size (nlat x nlat)

%%
%Created on 2016
%@author: Vignayanandam R. Muddapu (CNS@IIT-Madras)

%%
w = zeros(nlat, nlat);
ic = (nlat+1)/2;
jc = (nlat+1)/2;


for i = 1:nlat,
   for j = 1:nlat,
        dis = (i-ic)*(i-ic) + (j-jc)*(j-jc);
        w(i, j) = sig*exp(-dis/(rad*rad)) ; 
        if(i==j)
            w(i,j)=0;
        end
   end
end