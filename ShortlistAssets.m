function [ Shortlist, time ] = ShortlistAssets( Branch, EdgeList, LineOutNumber )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program shortlists the transmission  
% assets that must be evaluated by the FT following a branch 
% outage. The logic for this program is based on the SA 
% algorithm
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
    [row_renum,col_renum] = size(Branch);
    [row_list, col_list] = size(EdgeList);

    vec = zeros(row_list,1);
    EdgeList = horzcat(EdgeList,vec);


    [row, col] = size(EdgeList);
    count = 1;
    l_Col = [];
    l_Col = find(EdgeList(LineOutNumber,:)==0);
    count = 1;

    Shortlist = [];
    
    for eno = 1:length(EdgeList(:,1))
        flag = 0;   
        e_Col = [];
        e_Col = find(EdgeList(eno,:)==0);      
        for e_C = 1:2:e_Col(1)-2        
            for l_C = 1:2:l_Col(1)-2            
                if (EdgeList(LineOutNumber,l_C)==EdgeList(eno,e_C) && EdgeList(LineOutNumber,l_C+1)==EdgeList(eno,e_C+1))             
                    Shortlist(count,1) = eno;
                    Shortlist(count,2) = Branch(eno,1);
                    Shortlist(count,3) = Branch(eno,2);
                    count = count+1;
                    flag = 1;
                    break;                
                end
                if (EdgeList(LineOutNumber,l_C)==EdgeList(eno,e_C+1) && EdgeList(LineOutNumber,l_C+1)==EdgeList(eno,e_C))             
                    Shortlist(count,1) = eno;
                    Shortlist(count,2) = Branch(eno,1);
                    Shortlist(count,3) = Branch(eno,2);
                    count = count+1;
                    flag = 1;
                    break;                
                end
            end        
            if flag==1
                break;
            end        
        end    
    end
time = toc;
end

