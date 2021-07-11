function [ Flow, Capacity, A CL_Sp, EdgeList, PathStack, EdgeSatStack, CutsetStack, time ] = OutageAnalysis( Branch, Flow, Capacity, LineOutNumber, EdgeList, A )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program performs feasibility test (FT) 
% following a branch outage in the system. Therefore, it involves 
% the following:
% (a) The Update Scheme (UPS) for updating the weighted graphs
%     after the outage
% (b) The Shortlisting Assets (SA) algorithm to determine the 
%     assets which should be evaluated by FT
% (c) The feasibility test (FT) on the shortlisted set of assets 
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
    [ LoseFlag, PathAr, CurrentFlow, FlowCap, FlowInjAr ] = CheckIfLose_Cutset( Branch, LineOutNumber, Flow, Capacity, A );
%% The Update Scheme (UPS):
    if LoseFlag==1    
        [  Flow, Capacity ] = UpdateScheme( Branch, LineOutNumber, Flow, Capacity );    
        A(Branch(LineOutNumber,1),Branch(LineOutNumber,2)) = 0;
        A(Branch(LineOutNumber,2),Branch(LineOutNumber,1)) = 0;
    else 
        fprintf('\n Warning! Outage of the branch saturates a cut-set. \n');    
    end

%% The Shortlisting Assets (SA):
    Shortlist = ShortlistAssets( Branch, EdgeList, LineOutNumber );

% The branches with zero latent capacities are indentified, and 
% the latent capacities are increased by a small margin to ensure that all
% cut-sets are identified properly
    for Line = 1:length(Branch(:,1))
        if Capacity(Branch(Line,1),Branch(Line,2))==0 && Branch(Line,4)==1
            Capacity(Branch(Line,1),Branch(Line,2)) = 0.0001;
        elseif Capacity(Branch(Line,2),Branch(Line,1))==0 && Branch(Line,4)==1
            Capacity(Branch(Line,2),Branch(Line,1)) = 0.0001;
        end
    end

%% Feasibility Test (FT) on shortlisted assets:
    CL_Sp = [];
    PathStack = [];
    EdgeSatStack = [];
    CutsetStack = [];

    count = 1;
    [rowF, colF] = size(Shortlist);
    
    for i=1:rowF
    
        Line = Shortlist(i,1);
        FlagPresBefore = 0;
        
        if FlagPresBefore==0    
            [ LoseFlag, PathAr, CurrentFlow, FlowCap, FlowInjAr, flag_Radial, EdgeSat, Cutset ] = CheckIfLose_Cutset( Branch, Line, Flow, Capacity, A );
    
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
            if (LoseFlag==0) && (flag_Radial==0)
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
    end
time = toc;

end

