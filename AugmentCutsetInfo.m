function [ K, Tm, Cutset_FT, time ] = AugmentCutsetInfo(K_new,Tm_new,Cutset_FT_new,K,Tm,Cutset_FT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: The new violations identified by the FT, are 
% augmented with the violations detected in a previous iteration,  
% such that the corrective action can be initiated with respect 
% to all the violations
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    tic;

    % Initialization:
    NumOfCritCutset = length(K_new(:,1));
    MaxNumBranchCritCutset = length(K_new(1,:)); 

    % Make the transfer margins of the saturated cut-sets 
    % addressed in previous iteration as zero
    
    [ row_K_old, ~] = size(K);
    for i = 1:row_K_old
        Tm(i,1) = 0;
    end

    % Augment the new cutsets with their respective transfer margins
    ncutset = length(Tm)+1;
    [row_K, col_K] = size(K_new);
    for r = 1:length(K_new(:,1))   
        K(ncutset,[1:col_K]) = K_new(r,[1:col_K]);    
        Tm(ncutset,1) = Tm_new(r);    
        [row_set, col_set] = size(Cutset_FT_new(:,:,r));
        Cutset_FT([1:row_set],[1:col_set],ncutset) = Cutset_FT_new([1:row_set],[1:col_set],r); 
        ncutset = ncutset+1;
    end 

    time = toc;

end

