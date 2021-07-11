function [ Flow, Capacity, A, time ] = NetworkFlowAlgorithm( Bus, Branch, Gen, BusLoad )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Description: The network flow algorithm (NFA) creates the 
% "flow" and "latent capacity graphs" based upon the 
% conservation of energy.
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
%% Initialize arrays and variables
NoOfBus = length(Bus);
NoOfBranch = length(Branch);
BranchSt = ones(NoOfBranch,1);

GenPos = Gen(:,1);
Generation = Gen(:,2);
LoadPos = BusLoad(:,1);
Load = BusLoad(:,2);

NoOfBFS = 0;
problem = 0;

%% Initialize the "flow" and "latent capacity" graphs
Capacity = sparse(NoOfBus,NoOfBus); % Latent capacity graph
Flow = sparse(NoOfBus,NoOfBus); % Flow graph
A = sparse(NoOfBus,NoOfBus); % Incidence matrix
for k=1:NoOfBranch 
     if BranchSt(k)==1
        Capacity(Branch(k,1),Branch(k,2)) = Capacity(Branch(k,1),Branch(k,2)) + Branch(k,3);
        Capacity(Branch(k,2),Branch(k,1)) = Capacity(Branch(k,2),Branch(k,1)) + Branch(k,3);
        A(Branch(k,1),Branch(k,2)) = 1; A(Branch(k,1),Branch(k,1)) = 1; 
        A(Branch(k,2),Branch(k,1)) = 1; A(Branch(k,2),Branch(k,2)) = 1; 
     end
end
DontSelect = [];
countD = 1;

%% Create the "flow" and "latent capacity" graphs iteratively 
while (1<2)
      FF = CheckZeros(Load);
      GG = CheckZeros(Generation);
      if (Is_i_Present(0,FF)==1) || (Is_i_Present(0,GG)==1)
      else
        break;
      end
      
      for i = 1:length(Load)
          if Load(i)~=0
              Sink = LoadPos(i);
              break;
          end
      end
      
     % Selection of the source:         
        for j = 1:length(Generation)
            if Generation(j)~=0 && problem==0
                Source = GenPos(j);
                break;
            else
                if Generation(j)~=0 && IsPresent(DontSelect,j)==0
                    Source = GenPos(j);                    
                    break;
                end
            end
        end 
        
        % Finding the maximum power that can be injected along
        % the shortest path from the source to the sink        
        while (1<2)       
            [S,path]=graphshortestpath(Capacity,Source,Sink,'Method','BFS','Directed','true');NoOfBFS = NoOfBFS + 1;
            if S==Inf && IfCloseToZero(Load(i))==0 && IfCloseToZero(Generation(j))==0            
                problem = 1;
                DontSelect(countD) = j; countD = countD+1;
            else
                problem = 0;
                DontSelect = [];
                countD = countD+1;
            end
            if S==Inf || Load(i)==0 || Generation(j)==0
                break;
            end 
            MaxCap = 999999;
            for k=1:S
                From = path(k);To = path(k+1);          
                if MaxCap>Capacity(From,To)
                    MaxCap = Capacity(From,To);
                end
            end 
            
        % Determining the flow that will be injected along the path     
            if Load(i)<=MaxCap && Load(i)<=Generation(j)
                FlowInj = Load(i);
            elseif Generation(j)<=MaxCap && Generation(j)<=Load(i)
                FlowInj = Generation(j);
            elseif MaxCap<=Load(i) && MaxCap<=Generation(j)
                FlowInj = MaxCap;            
            end
            
        % Updating the source and sink values:
            Load(i) = Load(i)-FlowInj;
            Generation(j) = Generation(j)-FlowInj;
            
        % Updating the "flow" and "latent capacity" graph based upon the power
        % transferred along different paths
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

