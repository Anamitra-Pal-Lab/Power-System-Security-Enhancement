function [ LoseFlag, PathAr, CurrentFlow, FlowCap, FlowInjAr, flag_Radial, EdgeSat, Cutset] = CheckIfLose_Cutset( LinesArray, Line, Flow, Capacity, A )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program evaluates if a specific 
% transmission outage will create post-contingency cut-set 
% saturation, based upon the FT algorithm
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BusA = LinesArray(Line,1);
    BusB = LinesArray(Line,2);

    NewFlowSheet = Flow;
    NewFlowSheet(NewFlowSheet<0) = 0;

    [Bus1, Bus2, flow] = find(NewFlowSheet);
    found = 0;

% Obtain the current flow of the specified branch:
    for i = 1:length(Bus1)
        if ((Bus1(i)==BusA && Bus2(i)==BusB) || (Bus1(i)==BusB && Bus2(i)==BusA))
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

% Updates the incidence matrix accordingly:
    A(FromBus,ToBus) = 0;
    A(ToBus,FromBus) = 0;
    
% Checks if the outage branch is a radial branch or not
    [S,path]=graphshortestpath(A,FromBus,ToBus,'Method','BFS','Directed','true'); 
    flag_Radial = 0;
    if S==Inf    
        flag_Radial = 1; 
    end

% Remove the line from the latent capacity graph
    Capacity(FromBus,ToBus) = 0;
    Capacity(ToBus, FromBus) = 0;

%% Find the maximum power that can be transferred along the indirect paths 
    FlowCap = 0; 
    LoseFlag = 0;
    FlowInjAr = [];
    countP = 1;
    PathAr = [];
    EdgeSat = [];
    countS = 1;
% If the power flow through the "direct path" is more than the "maximum
% power that can be transferred through the "indirect paths", then it
% creates a saturated cut-set.
    while (1<2)       
      [S,path]=graphshortestpath(Capacity,FromBus,ToBus,'Method','BFS','Directed','true');  
      if S==Inf
            break;
      end 
      
      if S<Inf     
            PathAr(countP,[1:length(path)]) = path;            
            MaxCap = 9999;
            for k=1:S
                From = path(k);To = path(k+1);          
                if MaxCap>Capacity(From,To)
                    MaxCap = Capacity(From,To);
                end
            end
            FlowInj = MaxCap;
            for k=1:S
                From = path(k);To = path(k+1);
                Flow(From,To) = Flow(From,To) + FlowInj; 
                Flow(To,From) = Flow(To,From) - FlowInj;
                Capacity(From,To) = Capacity(From, To) - FlowInj;
                Capacity(To,From) = Capacity(To, From) + FlowInj;
                if Capacity(From,To)<0.0001
                    % Finding the saturated edges after flow injection:
                    EdgeSat(countS,1) = From;
                    EdgeSat(countS,2) = To;
                    countS = countS + 1;
                end               
            end
            FlowCap = FlowCap + FlowInj;
            FlowInjAr(countP,1) = FlowInj;
            countP = countP + 1;
            if FlowCap>=CurrentFlow  
                LoseFlag = 1;
                break;
            end
      end
    end
% Saturated cut-sets with a transfer margin lesser than 0.001 are ignored
    if (abs(CurrentFlow-FlowCap)<0.001)
        LoseFlag = 1;
    end

    Cutset = [];

%% Find the saturated cut-set:    
    if LoseFlag==0        
        V_insub = [];
        if LoseFlag==0 && flag_Radial==0
            % Group the vertices:
            [row_P, col_P] = size(PathAr);
            for i = 1:row_P
                V_insub = horzcat(V_insub,PathAr(i,:));        
            end    
            V_insub = unique(V_insub);
            V_insub(V_insub==0) = [];
        end

        V_reach_F = []; V_reach_T = []; 
        countF = 1; countT = 1;
        [row_V, col_V] = size(V_insub);
        for v = 1:col_V     
            [S,path]=graphshortestpath(Capacity,FromBus,V_insub(v),'Method','BFS','Directed','true');
            if S<Inf
                V_reach_F(countF,1) = V_insub(v); countF = countF+1;        
            else
                V_reach_T(countT,1) = V_insub(v); countT = countT+1;
            end    
        end

        K = 1;
        Cutset(K,1) = FromBus;Cutset(K,2) = ToBus; K = K+1;
        [row_E, col_E] = size(EdgeSat);

        for i = 1:row_E    
            F = EdgeSat(i,1); T = EdgeSat(i,2);    
            [ flag_F, pos ] = IsPresent( V_reach_F, F );    
            [ flag_T, pos ] = IsPresent( V_reach_T, T );   
            if flag_F==1 && flag_T==1 
                Cutset(K,1) = F;
                Cutset(K,2) = T;
                K = K+1;
            end    
        end
    end

end