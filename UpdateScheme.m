function [  Flow, Capacity, time ] = UpdateScheme( mpcNewbranch, Line, Flow, Capacity )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program updates the flow and 
% latent capacity graphs based upon the UPS algorithm. The 
% logic for this program is based on the UPS algorithm.
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;    
% Find the flow through the branch
    BusA = mpcNewbranch(Line,1);
    BusB = mpcNewbranch(Line,2);
    NewFlowSheet = Flow;
    NewFlowSheet(NewFlowSheet<0) = 0;
    [Bus1, Bus2, flow] = find(NewFlowSheet);
    found = 0;
    for i = 1:length(Bus1)
        if (Bus1(i)==BusA && Bus2(i)==BusB) || (Bus1(i)==BusB && Bus2(i)==BusA)
            FromBus = Bus1(i);
            ToBus = Bus2(i);
            CurrentFlow = flow(i);
            found = 1;
        end
    end
    if found==0
        FromBus = BusA;
        ToBus = BusB;
        CurrentFlow = 0;
    end

% Remove the branch from the flow and latent capacity graphs
    Flow(FromBus, ToBus) = 0;
    Flow(ToBus, FromBus) = 0;
    Capacity(FromBus,ToBus) = 0;
    Capacity(ToBus, FromBus) = 0;

% Re=route the flow through the set of indirect paths:
    FlowCap = 0; 
    LoseFlag = 0;
    EdgeTouch = [];
    TouchCount = 1;
    if CurrentFlow==0
        LoseFlag = 1;
    else
        countP = 1;
        while (1<2)       
            [S,path]=graphshortestpath(Capacity,FromBus,ToBus,'Method','BFS','Directed','true');
      
            if S==Inf
                break;
            end 
            MaxCap = 9999;
            for k=1:S
                From = path(k);To = path(k+1);          
                if MaxCap>Capacity(From,To)
                    MaxCap = Capacity(From,To);
                end
            end
       
            FlowInj = MaxCap;
            if FlowInj>CurrentFlow
                FlowInj = CurrentFlow;
            end
                  
            for k=1:S
                From = path(k);To = path(k+1);                                      
                Flow(From,To) = Flow(From,To) + FlowInj; 
                Flow(To,From) = Flow(To,From) - FlowInj;
                Capacity(From,To) = Capacity(From, To) - FlowInj;
                Capacity(To,From) = Capacity(To, From) + FlowInj;
            end
            countP = countP + 1;
       
            CurrentFlow = CurrentFlow-FlowInj;
            if CurrentFlow==0
                break;
            end
        end
    end
time = toc;
end

