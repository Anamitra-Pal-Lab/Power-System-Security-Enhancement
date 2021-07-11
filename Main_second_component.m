%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Name: Main_program_second_component
% 
% Program Description: This program implements the proposed second 
% component (FT-rCA) during successive outages in a power network
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear all
close all

%% Load the input data
mpc = loadcase('case118_J2.m');
load Data_118bus_J2.mat;

%% Initialize different matrices:
BusGraph = Bus;
BranchGraph = Branch(:,[1:2]);BranchGraph(:,3) = Branch(:,7);BranchGraph(:,4) = Branch(:,8);
GeneratorGraph = Generator(:,[1:2]);
LoadGraph = Load(:,[1:2]);
GenBusNumAr = Generator(:,1);
PgenOldAr = Generator(:,2);  

%% Build the "flow" and "latent capacity" graphs:
% The graphs are built based on the graph-theory based 
% network flow algorithm
[ Flow, Capacity, A, ~ ] = NetworkFlowAlgorithm(BusGraph,BranchGraph,GeneratorGraph,LoadGraph);

% An alternate way of building the "flow" and "latent capacity" 
% graphs is to use a DC power flow solution in the base-case scenario
% Flow = sparse(NoOfBus,NoOfBus);
% Capacity = sparse(NoOfBus,NoOfBus);
% for i = 1:length(Branch(:,1))       
%     Flow(Branch(i,1),Branch(i,2)) = Branch(i,6); % dc power flows 
%     Flow(Branch(i,2),Branch(i,1)) = (-1)*Branch(i,6);
%     Capacity(Branch(i,1),Branch(i,2)) = Branch(i,7)-Flow(Branch(i,1),Branch(i,2));    
%     Capacity(Branch(i,2),Branch(i,1)) = Branch(i,7)-Flow(Branch(i,2),Branch(i,1));    
% end


%% Initialize different variables:
fprintf('\n------ System condition: Base-case (No outage) ------\n');
BranchOut = [];
count = 1;
flag_vio_out = 1;
Net_change_cost = 0;
continue_flag = 1;
K = [];
Cutset_FT = [];
Total_load_base = sum(Load(:,2));

%% Create the PTDF, LODF, B and H matrices: 
fprintf('Creating the PTDF matrix \n');
baseMVA = 100;
NoOfBus = length(Bus);
RoundOffFlag = 0;
[PTDF_true, PTDF, LODF, B_full, H_full, ~] = Create_PTDF_LODF_B_H(Bus,Branch,RoundOffFlag);

%% Perform feasibility test (FT) for all branches in the base-case scenario:
[ CL_Sp_vio, CutsetStack_vio, EdgeList, flag_vio, ~ ] = FeasibilityTestBasecase( Flow, Capacity, A, BranchGraph );

fprintf('-------------------------------------------- \n');

initial = 1;

while (continue_flag==1)
    
if (flag_vio==1)
    %% Display the violations detected by the FT algorithm:
    DisplayViolations_FT( CL_Sp_vio,CutsetStack_vio );               
              
    %% Create inputs for the relaxed corrective action (rCA):
    if ((initial==1) || (isempty(K)==1))
        [ K, Tm, Cutset_FT, ~ ] = CreateInput_cutset( CL_Sp_vio, CutsetStack_vio, Branch);             
    else            
        [ K, Tm, Cutset_FT, ~] = AugmentCutsetInfo(K_new,Tm_new,Cutset_FT_new,K,Tm,Cutset_FT);
    end
    
    while (flag_vio==1)
%% Perform the relaxed corrective action (rCA):                     
           [ GeneratorNegativeChange, GeneratorPositiveChange, LoadNegativeChange, LoadPositiveChange, Branch, Load, Generator, Soln_Flag, Flow_dc, Rate_dc, tot_change_cost, ~ ] = RelaxedCorrectiveAction( K, Tm, Cutset_FT, PTDF, Bus, Branch, Generator, Load );           
           count = count + 1;    
           % Group all the injection increase and injection decrease 
           % together in separate matrices
           InjectionPositiveChange = [];
           Temp = []; Temp = LoadNegativeChange; Temp(:,2) = (-1)*Temp(:,2); 
           InjectionPositiveChange = vertcat(GeneratorPositiveChange,Temp);            
           Temp = [];Temp = LoadPositiveChange; 
           InjectionNegativeChange = vertcat(GeneratorNegativeChange,Temp);           
           if Soln_Flag==0
               % This implies that the optimization problem in the rCA has 
               % not converged and there is a problem
                break;
           end
%% Update the system based upon the redispatch solution:           
           mpc.gen(:,2) = Generator(:,2);
           for nload = 1:length(Load(:,1))
                LoadBus = Load(nload,1);
                loc = find(mpc.bus(:,1)==LoadBus);
                mpc.bus(loc,3) = Load(nload,2);
           end
           Res = rundcpf(mpc);
           Branch(:,6) = Res.branch(:,14);           
           
%% Update the "flow" and "latent capacity graphs" following redispatch based on M-UPS:
           [Flow,Capacity,BranchFlowChange,~] = ModifiedUpdateScheme(Flow, Capacity, InjectionPositiveChange, InjectionNegativeChange, BranchGraph);           
            
%% Perform shortlisting assets following redispatch based on M-SA algorithm:
           [ShortlistedEdges, ~] = ModifiedShortlistAssets(BranchFlowChange, EdgeList, BranchGraph);          
                                
%% Perform Feasibility Test (FT) on shortlisted assets:
           [CL_Sp_vio, CutsetStack_vio, EdgeList, flag_vio, ~] = FeasibilityTestOnShortlist( Flow, Capacity, A, BranchGraph, ShortlistedEdges, EdgeList );
                
           if (flag_vio==1)
%% Display additional violations (if any) due to the redispatch:
               DisplayViolations_FT( CL_Sp_vio,CutsetStack_vio );               
%% Create the inputs for rCA to mitigate the combined violations:
               [K_new, Tm_new, Cutset_FT_new, ~ ] = CreateInput_cutset( CL_Sp_vio, CutsetStack_vio, Branch); 
               [K, Tm, Cutset_FT, ~] = AugmentCutsetInfo(K_new,Tm_new,Cutset_FT_new,K,Tm,Cutset_FT);
           end    
    end        
    if (flag_vio==0)
%% Display the results of the corrective action:       
         fprintf('-----------------------------------------------------------------\n');
         fprintf('The second component has alleviated all post-contingency cut-set saturation \n');           
         fprintf('-----------------------------------------------------------------\n');
         GeneratorCost_Ar = Generator(:,6).*Generator(:,2).^2+Generator(:,5).*Generator(:,2); TotalGeneratorCost = sum(GeneratorCost_Ar);
         LoadCost_Ar = (Load(:,2)-Load(:,4)).*Load(:,3); TotalLoadCost = sum(LoadCost_Ar);            
         fprintf('Production cost = $ %f \n',TotalGeneratorCost);                 
         Net_load_shed = Total_load_base-sum(Load(:,2));
         fprintf('Total amount of load shed = %f MW \n',Net_load_shed);         
    else
         fprintf('-----------------------------------------------------------------\n');   
         fprintf('Warning: All post-contingency cut-set saturation cannot be mitigated \n');         
         fprintf('-----------------------------------------------------------------\n');                           
    end
    
else
    fprintf('-----------------------------------------------------------------\n');
    fprintf('No violations are detected by the FT algorithm \n');
    fprintf('-----------------------------------------------------------------\n');
end

%% Check if there are successive branch outages in the system:
   LineOutNumber = input('\n Enter the branch number which is out (Press 0 and enter if you do not want to continue) ?');
   if LineOutNumber==0
      continue_flag = 0;
      break;
        
   else
      fprintf('\n------ System condition: Outage of branch (%d-%d)------\n',Branch(LineOutNumber,1),Branch(LineOutNumber,2));            

%% Update the system matrices following the branch outage:      
        mpc.branch(LineOutNumber,11) = 0;        
% Update the system matrices, instead of re-building 
% the matrices from scratch
        [ PTDF_true, PTDF, LODF, B_full, H_full, ~ ] = Update_PTDF_LODF_B_H( B_full, H_full, Bus, Branch, LineOutNumber,RoundOffFlag);                                 
        Res = rundcpf(mpc);                
        Branch(:,6) = Res.branch(:,14);
        Branch(:,8) = Res.branch(:,11);
        BranchGraph(:,4) = Branch(:,8);
        
        
%% A successive FT has to be performed following the branch outage
% However, the successive FT must involve the UPS algorithm and the 
% SA algorithm for fast computation
      [ Flow, Capacity, A, CL_Sp_vio, EdgeList, PathStack, EdgeSatStack, CutsetStack_vio, ~ ] = OutageAnalysis( BranchGraph, Flow, Capacity, LineOutNumber, EdgeList, A );     
      [ row_vio, col_vio ] = size(CL_Sp_vio);
      if (row_vio>=1)
          flag_vio = 1;          
          initial = 0;                   
% Creates input for the next relaxed corrective action (rCA):
          [ K_new, Tm_new, Cutset_FT_new, ~ ] = CreateInput_cutset( CL_Sp_vio, CutsetStack_vio, Branch);          
      else
          flag_vio = 0;  
      end                  
   end    
end
