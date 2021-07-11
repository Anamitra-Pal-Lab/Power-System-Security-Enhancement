function [ CL_Sp_vio, CutsetStack_vio, EdgeList, flag_vio, time ] = FeasibilityTestOnShortlist( Flow, Capacity, A, Branch, ShortlistedEdges, EdgeList )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: This program performs the feasibility test
% (FT) algorithm on the shortlisted assets following a change in 
% generation redispatch in the system
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization:
CL_Sp = [];
PathStack = [];
EdgeSatStack = [];
CutsetStack = [];
tic;
count = 1;
[rowF, colF] = size(ShortlistedEdges);

% Evaluate the shortlisted branches by the FT algorithm:
for i=1:rowF    
    Line = ShortlistedEdges(i,1);    
    [ LoseFlag, PathAr, CurrentFlow, FlowCap, FlowInjAr, flag_Radial, EdgeSat, Cutset ] = CheckIfLose3_Break_Cutset( Branch, Line, Flow, Capacity, A );    
    [ row, col ] = size(PathAr);
    EdgeCount = 1;
    EnterLoop = 0;
    EdgeList(Line,:) = zeros(1,length(EdgeList(Line,:)));
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
    [row, col] = size(PathAr);
    [row_e, col_e] = size(EdgeSat);
    [row_K,col_K] = size(Cutset);
    if LoseFlag==0
        PathInterest([1:row],[1:col],count) = PathAr;
        CL_Sp(count,1) = Line;
        CL_Sp(count,2) = Branch(Line,1);
        CL_Sp(count,3) = Branch(Line,2);
        CL_Sp(count,4) = FlowCap-CurrentFlow;
        CL_Sp(count,5) = flag_Radial;
        PathStack(1:row,1:col,count) = PathAr;        
        EdgeSatStack(1:row_e,1:col_e,count) = EdgeSat;
        CutsetStack([1:row_K],[1:col_K],count) = Cutset;
        count = count + 1;
    end
end

% Check if there are non-singleton violations:
flag_vio = 0;
count = 1;
CL_Sp_vio = [];
CutsetStack_vio = [];
[row,col] = size(CL_Sp);
for i = 1:row
    value = CL_Sp(i,5);
    if value==0
        flag_vio = 1;
        CL_Sp_vio(count,:) = CL_Sp(i,:);
        [r,c] = size(CutsetStack(:,:,i));
        CutsetStack_vio([1:r],[1:c],count) = CutsetStack([1:r],[1:c],i);
        count = count + 1;
    end
end
time = toc;

end
