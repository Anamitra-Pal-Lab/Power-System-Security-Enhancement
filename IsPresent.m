function [ flag, pos ] = IsPresent( Arr, Val )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: Checks if a given number is contained    
% in a specific array
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag = 0;
pos = 0;
for i = 1:length(Arr)
    if Val==Arr(i)
        flag = 1;
        pos = i;
        break;
    end
end


end

