function [ PTDF_true, PTDF_approx, LODF, B_full, H_full, time ] = Update_PTDF_LODF_B_H( B_full, H_full, Bus, Branch, BranchOut,RoundOffFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: This program updates the system matrices
% following a branch outage in the system. 
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
%% Updating the H matrix:
    [row_H, col_H] = size(H_full);    
    for i = 1:length(BranchOut)
        BranchNum = BranchOut(i);
        F = Branch(BranchNum,1);
        T = Branch(BranchNum,2);        
        H_full(BranchNum,:) = zeros(1, col_H);        
    end
    
%% Updating the B matrix:    
    % Updating the susceptance (B) matrix changes only four entries of the
    % matrix and hence saves computation time
    for i = 1:length(BranchOut)
        BranchNum = BranchOut(i);
        F = Branch(BranchNum,1);
        T = Branch(BranchNum,2);
        B_full(F,F) = B_full(F,F)-abs(B_full(F,T));
        B_full(T,T) = B_full(T,T)-abs(B_full(T,F));
        B_full(F,T) = 0; 
        B_full(T,F) = 0; 
    end
    
%% Finding the new PTDF matrix:
    noofbus = length(Bus);
    B = B_full([1:noofbus-1],[1:noofbus-1]);    
    H = H_full(:,1:noofbus-1);     
    
% Perform matrix operation to obtain the PTDF matrix:
    X = inv(B);
    PTDF = H*X; 
    PTDF_true = PTDF;
% For all PTDF values lesser than 0.02, round them down to zero;
    if (RoundOffFlag==1)
       [r,c] = find(abs(PTDF)<0.02);
        for i = 1:length(r)
            PTDF(r(i),c(i)) = 0;
        end
    end
    PTDF_approx = PTDF;
    
% From the PTDF matrix, we now create the LODF matrix:        
    PTDF_full = horzcat(PTDF,zeros(length(Branch(:,1)),1));
    [nl, nb] = size(PTDF_full);
    f = Branch(:, 1);
    t = Branch(:, 2);
    Cft =  sparse([f; t], [1:nl 1:nl]', [ones(nl, 1); -ones(nl, 1)], nb, nl);
    H = PTDF_full * Cft;
    h = diag(H, 0);
    LODF = H ./ (ones(nl, nl) - ones(nl, 1) * h');
    h_diff = abs(ones(length(h),1)-h);
    [ pos_ar ] = find(h_diff<0.00001);
    LODF = LODF - diag(diag(LODF)) - eye(nl, nl);
    for  i = 1:length(pos_ar) 
       pos_val = pos_ar(i); 
       LODF([1:nl],pos_val) = zeros(nl,1);
       LODF(pos_val,[1:nl]) = zeros(1,nl);
       LODF(pos_val,pos_val) = -1;
    end        
time = toc;
        
end

