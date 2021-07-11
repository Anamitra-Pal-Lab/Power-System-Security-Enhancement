function [ out ] = DisplayViolations_FT( CL_Sp_vio,CutsetStack_vio )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: This program displays the violations 
% (post-contingency cut-set saturation) identified by the FT 
% algorithm
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ row_Cl,~ ] = size(CL_Sp_vio);
if row_Cl>0
    fprintf('-------------------------------------------- \n');
    fprintf('Contingencies that create saturated cut-sets: \n');
    fprintf('-------------------------------------------- \n');
    for i = 1:length(CL_Sp_vio(:,1))
        fprintf('Case %d :',i);
        fprintf('Outage of %d-%d saturates cut-set K%d by %f MW, where K%d={',CL_Sp_vio(i,2),CL_Sp_vio(i,3),i,CL_Sp_vio(i,4),i);     
        [row,col] = size(CutsetStack_vio(:,:,1));
        for j = 1:row
            F = CutsetStack_vio(j,1,i);
            T = CutsetStack_vio(j,2,i);
            LastFlag = 0;
            if j==row
                LastFlag = 1;
            else     
                if CutsetStack_vio(j+1,1,i)==0
                    LastFlag = 1;
                end
            end
        
            if F>0
                if (LastFlag==0)
                    fprintf('%d-%d,',F,T);
                else
                    fprintf('%d-%d',F,T);
                end
            end
        end
        fprintf('} \n');
    end
    fprintf('-------------------------------------------- \n');    
else
   fprintf('-------------------------------------------- \n');
   fprintf('No contingencies create saturated cut-sets: \n');    
   fprintf('-------------------------------------------- \n');
end
out = 1;

end

