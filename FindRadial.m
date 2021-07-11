function [ Radial, time ] = FindRadial( Branch, A )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program finds the list of radial
% branches in the system
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
    Radial = [];
    count = 1;
    for i = 1:length(Branch(:,1))
        Fbus = Branch(i,1);
        Tbus = Branch(i,2);
        if (Branch(i,8)==1)
            A(Fbus,Tbus) = 0;
            A(Tbus,Fbus) = 0;        
            [S,path]=graphshortestpath(A,Fbus,Tbus,'Method','BFS','Directed','true');  
            if (S==Inf)
                Radial(count,1) = i;
                Radial(count,2) = Fbus;
                Radial(count,3) = Tbus;
                count = count + 1;
            end 
            A(Fbus,Tbus) = 1;
            A(Tbus,Fbus) = 1;  
        end
    end
time = toc;

end

