function [CL_Sp_vio, CutsetStack_vio, EdgeList, flag_vio, time] = FeasibilityTestBasecase( Flow, Capacity, A, Branch )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program implements the feasibility 
% test (FT) algorithm in the base-case scenario for all 
% transmission assets
% 
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
% Find the lines which have a latent capacity of zero, 
% and increase its capacity by a small margin
for Line=1:length(Branch(:,1))    
    if Capacity(Branch(Line,1),Branch(Line,2))==0
        Capacity(Branch(Line,1),Branch(Line,2)) = 0.0001;
    elseif Capacity(Branch(Line,2),Branch(Line,1))==0
        Capacity(Branch(Line,2),Branch(Line,1)) = 0.0001;
    end
end

%% Analyzing different transmission assets by the FT algorithm:
CL_Sp = [];
count = 1;
CL_Sp_vio = [];
CutsetStack_vio = [];
EdgeList = zeros(length(Branch(:,1)),1);
count_radial = 1;
for Line=1:length(Branch(:,1))         
        FlagPresBefore = 0;        
        if FlagPresBefore==0    
            [ LoseFlag, PathAr, CurrentFlow, FlowCap, FlowInjAr, flag_Radial, EdgeSat, Cutset ] = CheckIfLose_Cutset(Branch, Line, Flow, Capacity, A);
            [ row, col ] = size(PathAr);
            EdgeCount = 1;
            EnterLoop = 0;
            for R = 1:row
                for C = 1:col-1
                    if PathAr(R,C+1)>0
                        PresentFlag = 0;
                        if EnterLoop==1
                            Col_list = length(EdgeList(Line,:));                
                            for k = 1:Col_list-1
                                if EdgeList(Line,k)==PathAr(R,C) && EdgeList(Line,k+1)==PathAr(R,C+1)
                                    PresentFlag=1; 
                                end                    
                            end
                        end
                
                        if PresentFlag==0
                            EdgeList(Line,EdgeCount) = PathAr(R,C);
                            EdgeList(Line,EdgeCount+1) = PathAr(R,C+1);
                            EdgeCount = EdgeCount+2;
                        end
                        EnterLoop = 1;
                    end
                end
            end
            
            if LoseFlag==0          
                CL_Sp(count,1) = Line;
                CL_Sp(count,2) = Branch(Line,1);
                CL_Sp(count,3) = Branch(Line,2);
                CL_Sp(count,4) = FlowCap-CurrentFlow;
                CL_Sp(count,5) = flag_Radial;
                [row,col] = size(PathAr);
                PathStack([1:row],[1:col],count) = PathAr;
                [row_e,col_e] = size(EdgeSat);
                EdgeSatStack([1:row_e],[1:col_e],count) = EdgeSat;
                [row_K,col_K] = size(Cutset);
                CutsetStack([1:row_K],[1:col_K],count) = Cutset;
                count = count + 1;
            end 
            NoOfPaths(Line,1) = size(PathAr,1);
        end            
end

%% Check if there are non-radial special assets detected by the FT algorithm
flag_vio = 0;
count = 1;

if (isempty(CL_Sp)==1)
    CL_Sp_vio = [];
    CutsetStack_vio = [];
else
    for i = 1:length(CL_Sp(:,1))
        value = CL_Sp(i,5);
        if value==0
            flag_vio = 1;
            CL_Sp_vio(count,:) = CL_Sp(i,:);
            [r,c] = size(CutsetStack(:,:,i));
            CutsetStack_vio([1:r],[1:c],count) = CutsetStack([1:r],[1:c],i);
            count = count + 1;
        end
    end
end    

fprintf('\n Total number of paths traversed = %d \n',sum(NoOfPaths));
num_non_zero = length(find(NoOfPaths~=0));
fprintf('Average number of paths traversed = %d \n',sum(NoOfPaths)/num_non_zero);

time = toc;

end
