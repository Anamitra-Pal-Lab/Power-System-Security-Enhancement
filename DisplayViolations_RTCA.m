function [out] = DisplayViolations_RTCA(CL_Sp_vio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program displays the violations 
% detected by the RTCA. RTCA identifies critical contingencies 
% that create post-contingency branch overloads.
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    fprintf('-------------------------------------------- \n');
    fprintf('Contingencies that create post-contingency branch overloads are as follows: \n');
    fprintf('-------------------------------------------- \n');
    [row, ~] = size(CL_Sp_vio);
    for i = 1:row
       fprintf('%f-%f \n',CL_Sp_vio(i,2),CL_Sp_vio(i,3));        
    end    
    out = 1;
    fprintf('-------------------------------------------- \n');
end

