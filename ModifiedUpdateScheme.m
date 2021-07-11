function [ Flow, Capacity, BranchFlowChange, time ] = ModifiedUpdateScheme( Flow, Capacity, GeneratorPositiveChange, GeneratorNegativeChange, Branch )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: This program creates an updated "flow" 
% and "latent capacity graph" after change in generation 
% in the system. The logic for this program is based on the
% M-UPS algorithm.
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
% Initialization:
    DontSelect = [];
    countD = 1;
    problem = 0;
    NoOfBFS = 0;
    EdgesFlowChange = []; 
    count = 1;
    countBT = 1;
    BranchFlowChange = [];
    countChange = 1;
    BranchTouch = [];
    FlowOrg = Flow;
    CapacityOrg = Capacity;
    GenIncData = GeneratorPositiveChange;
    GenDecData = GeneratorNegativeChange;
    GenPosInc = GenIncData(:,1);GenInc = GenIncData(:,2);
    GenPosDec = GenDecData(:,1);GenDec = abs(GenDecData(:,2));
    
% Select Source-Sink pairs from generator increase and decrease pairs to update the flow and capacity graphs:
    while (1<2)
      FF = CheckZeros(GenDec);
      GG = CheckZeros(GenInc);
      if (Is_i_Present(0,FF)==1) || (Is_i_Present(0,GG)==1)
      else
        break;
      end
      
      if (sum(GenDec)<0.01 && sum(GenInc)<0.01)
        break;
      end
      
      for i = 1:length(GenDec)
          if GenDec(i)~=0
              Sink = GenPosDec(i);
              break;
          end
      end
      
     %% Selection of the source:
     % Select a "source" depending upon the position of the "sink"    
        for j = 1:length(GenInc)
            if GenInc(j)~=0 && problem==0
                Source = GenPosInc(j);
                break;
            else
                if GenInc(j)~=0 && IsPresent(DontSelect,j)==0
                    Source = GenPosInc(j);                    
                    break;
                end
            end
        end 
        
% Finding the shortest path from the Source to the Sink and finding out the maximum capacity of the path.        
        while (1<2)       
            [S,path]=graphshortestpath(Capacity,Source,Sink,'Method','BFS','Directed','true');NoOfBFS = NoOfBFS + 1;
            if S<Inf
                for ii = 1:S
                    F = path(ii);T=path(ii+1);
                    len = length(BranchFlowChange);
                    if len>0
                        Present = 0;
                        for kk = 1:2:len
                            Fbr = BranchFlowChange(1,kk);
                            Tbr = BranchFlowChange(1,kk+1);                                                        
                            if (F==Fbr && T==Tbr) || (F==Tbr && T==Fbr)
                                Present = 1;
                                break;                            
                            end
                        end
                        if Present==0
                            BranchFlowChange(1,countChange) = F;
                            BranchFlowChange(1,countChange+1) = T;
                            countChange = countChange+2;
                        end
                    else
                        BranchFlowChange(1,countChange) = F;
                        BranchFlowChange(1,countChange+1) = T;
                        countChange = countChange+2;
                    end
                end                                    
            end
                        
            if ((S==Inf) && (IfCloseToZero(GenDec(i))==0) && (IfCloseToZero(GenInc(j))==0))            
                problem = 1;
                DontSelect(countD) = j; countD = countD+1;
            else
                problem = 0;
                DontSelect = [];
                countD = countD+1;
            end
            if ((S==Inf) || (GenDec(i)==0) || (GenInc(j)==0))
                break;
            end 
            MaxCap = 999999;
            for k=1:S
                From = path(k);To = path(k+1);          
                if MaxCap>Capacity(From,To)
                    MaxCap = Capacity(From,To);
                end
            end 
% Determine the flow injection along a given path
            if ((GenDec(i)<=MaxCap) && (GenDec(i)<=GenInc(j)))
                FlowInj = GenDec(i);
            elseif ((GenInc(j)<=MaxCap) && (GenInc(j)<=GenDec(i)))
                FlowInj = GenInc(j);
            elseif ((MaxCap<=GenDec(i)) && (MaxCap<=GenInc(j)))
                FlowInj = MaxCap;            
            end            
% Update the load and generation values
            GenDec(i) = GenDec(i)-FlowInj;
            GenInc(j) = GenInc(j)-FlowInj;
% Update the "flow" and "latent capacity" graphs for
% power injection along the given path
            for k=1:S
                From = path(k);To = path(k+1);
                Flow(From,To) = Flow(From,To) + FlowInj; 
                Flow(To,From) = Flow(To,From) - FlowInj;
                Capacity(From,To) = Capacity(From, To) - FlowInj;
                Capacity(To,From) = Capacity(To, From) + FlowInj;
            end   
        end              
    end    
time = toc;    
end

