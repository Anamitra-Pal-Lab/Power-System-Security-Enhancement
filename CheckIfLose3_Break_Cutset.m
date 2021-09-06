function [ LoseFlag, PathAr, CurrentFlow, FlowCap, FlowInjAr, flag_Radial, EdgeSat, Cutset] = CheckIfLose3_Break_Cutset( LinesArray, Line, Flow, Capacity, A )
%% Algorithm to Determine if you can loose a line

%% Determine the "From Bus" and the "To Bus" of the specific line.

% BusA = mpcNew.branch(Line,1);
% BusB = mpcNew.branch(Line,2);

BusA = LinesArray(Line,1);
BusB = LinesArray(Line,2);

NewFlowSheet = Flow;
NewFlowSheet(NewFlowSheet<0) = 0;

[Bus1, Bus2, flow] = find(NewFlowSheet);
found = 0;


for i = 1:length(Bus1)
    if ((Bus1(i)==BusA && Bus2(i)==BusB) || (Bus1(i)==BusB && Bus2(i)==BusA))
        FromBus = Bus1(i);
        ToBus = Bus2(i);
        CurrentFlow = flow(i);
        found = 1;
    end
end

%%Special Case:
% FromBus_n = ToBus;
% ToBus_n = FromBus;
% CurrentFlow_n =-CurrentFlow; 
% 
% FromBus = FromBus_n;
% ToBus = ToBus_n;
% CurrentFlow =CurrentFlow_n; 


if found==0
   FromBus = BusA;
   ToBus = BusB;
   CurrentFlow = 0;
end

%% Find if Radial Branch or Not:
A(FromBus,ToBus) = 0;
A(ToBus,FromBus) = 0;
%% 
[S,path]=graphshortestpath(A,FromBus,ToBus,'Method','BFS','Directed','true');
%% 
flag_Radial = 0;
if S==Inf    
   flag_Radial = 1; 
end

%% Remove the line from the Capacity graph
Capacity(FromBus,ToBus) = 0;
Capacity(ToBus, FromBus) = 0;

%% Check how much power you can transfer from the "FromBus" to "ToBus" ?
FlowCap = 0; 
LoseFlag = 0;
FlowInjAr = [];
countP = 1;
% if CurrentFlow==0
%     LoseFlag = 1;
%     PathAr = [];
%     FlowInjAr = [];
% else

PathAr = [];
EdgeSat = [];
countS = 1;
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
              % Finding saturated edges during the path search itself:
%               EdgeSat(countP,1) = From;
%               EdgeSat(countP,2) = To;
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
if (abs(CurrentFlow-FlowCap)<0.001)
    LoseFlag = 1;
end

Cutset = [];

if LoseFlag==0
%% Logic 3: Find the Cut-set

V_insub = [];
if LoseFlag==0 && flag_Radial==0
    % Group vertices:
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
    
% The above is a better way to group the vertices rather than what is given
% below.
    
%     [S,path]=graphshortestpath(Capacity,V_insub(v),ToBus,'Method','BFS','Directed','true');
%     if S<Inf
%         V_reach_T(countT,1) = V_insub(v); countT = countT+1;        
%     end
end

% Screen out the saturated edges among which: "whose one end belongs to "V_reach_F" and the other end belongs to "V_reach_T"
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