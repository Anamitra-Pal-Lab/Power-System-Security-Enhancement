function [ Shortlist, time ] = ModifiedShortlistAssets( BranchFlowChange, EdgeList, Branch )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: This program shortlists the transmission 
% assets that must be re-evaluated by the feasibility test (FT)
% algorithm following a generation redispatch in the system.
% This logic for this program is based on the M-SA algorithm.
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
len = length(BranchFlowChange);
[row_lim,col_lim] = size(EdgeList);
countF = 1;
Shortlist = [];

for r = 1:row_lim   
    Common = 0;
    if (EdgeList(r,1)~=0)
        for i = 1:2:len
            F = BranchFlowChange(1,i);
            T = BranchFlowChange(1,i+1);        
        
            col_lim = length(find(EdgeList(r,:)~=0)); % This line is newly added
            for c = 1:2:col_lim-1
                if (F==EdgeList(r,c) && T==EdgeList(r,c+1)) || (F==EdgeList(r,c+1) && T==EdgeList(r,c))
                    Common = 1;
                    break;           
                end
            end
            if (F==Branch(r,1) && T==Branch(r,2)) || (F==Branch(r,2) && T==Branch(r,1))
                Common = 1;               
            end
            if Common==1
                break;
            end        
        end    
        if (Common==1)
            Shortlist(countF,1) = r;
            Shortlist(countF,2) = Branch(r,1);
            Shortlist(countF,3) = Branch(r,2);
            countF = countF + 1;
        end
    end
end
time = toc;

end

