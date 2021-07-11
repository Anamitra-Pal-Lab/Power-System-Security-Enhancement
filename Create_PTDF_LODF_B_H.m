function [ PTDF_true, PTDF_approx, LODF, B_full, H_full, time ] = Create_PTDF_LODF_B_H(Bus,Branch,RoundOffFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: This function creates the power transfer
% distribution factor (PTDF), line outage distribution factor (LODF), the
% susceptance matrix (B) and the branch-bus matrix (H)
% 
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    tic;        
    %% Create the H (branch-bus matrix):
    H = [];
    for i = 1:length(Branch(:,1))
        Status = Branch(i,8);
        if Status==1
           FromBus = Branch(i,1);
           ToBus = Branch(i,2);
           H(i,FromBus) = 1/Branch(i,4);
           H(i,ToBus) = (-1)*1/Branch(i,4);
        else
            FromBus = Branch(i,1);
            ToBus = Branch(i,2);    
            H(i,FromBus) = 0;
            H(i,ToBus) = 0;
        end
    end        
    
    %% Build the B matrix after monitoring the branch statuses:
    B = zeros(length(Bus(:,1)));
    BusOld = Bus(:,1);
    for i = 1:length(Branch(:,1))
        Status = Branch(i,8);
        if Status==1
            FromBus = Branch(i,1);
            ToBus = Branch(i,2);
            xline = Branch(i,4);           
            B(FromBus,ToBus) = B(FromBus,ToBus)-1/(xline);
            B(ToBus,FromBus) = B(ToBus,FromBus)-1/(xline);
            B(FromBus,FromBus) = B(FromBus,FromBus)+1/(xline);
            B(ToBus,ToBus) = B(ToBus,ToBus)+1/(xline);                
        end    
    end
    
    %% Adjust the B and H matrices to account for the reference bus
    % B matrix: Remove the entire row and column for the reference bus 
    % H matrix: Remove the reference bus column from the H matrix
    noofbus = length(B(:,1));
    B_full = B;
    H_full = H;
    B_temp = B([1:noofbus-1],[1:noofbus-1]);
    B = B_temp;
    H_temp = H(:,1:noofbus-1); 
    H = H_temp;
    
   %% Perform matrix operation to obtain the PTDF matrix
    X = inv(B);
    PTDF = H*X;
    PTDF_true = PTDF;
    
    if (RoundOffFlag==1)        
        [r,c] = find(abs(PTDF)<0.02);
        for i = 1:length(r)
            PTDF(r(i),c(i)) = 0;
        end
    end
    PTDF_approx = PTDF;
    
    %% From the PTDF matrix, we now create the LODF matrix        
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