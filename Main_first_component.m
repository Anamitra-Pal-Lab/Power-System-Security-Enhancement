%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Program Name: Main_program_first_component
% 
% Program Description: This program implements the proposed first 
% component (FT-RTCA-iCA) during successive outages in a power network
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
%% Load the input data:
mpc = loadcase('case118_J2.m');
load Data_118bus_J2.mat;

%% Initialize different matrices:
Generator(:,4) = zeros(length(Generator),1); 
loc_negative = find(Generator(:,2)<0);
Generator(loc_negative,3) = Generator(loc_negative,2);
Generator(loc_negative,4) = Generator(loc_negative,2);
initial = 1;
continue_flag = 1;
K = [];
baseMVA = 100;
NoOfBus = length(Bus);
BusGraph = Bus;
BranchGraph = Branch(:,[1:2]);BranchGraph(:,3) = Branch(:,7);BranchGraph(:,4) = Branch(:,8);
GeneratorGraph = Generator(:,[1:2]);
LoadGraph = Load(:,[1:2]);

%% Settings: 
% Rank_limit controls the size of the contingency list 
% used in RTCA
Rank_limit = 54; 
% RoundOffFlag determines if the PTDF values will be 
% approximated below a specified threshold.
% 0: no approximation; 1: approximation  
RoundOffFlag = 0; 

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

%% Find the list of radial branches in the system:
[ Radial, ~ ] = FindRadial( Branch, A );

%% Create the PTDF, LODF, B and H matrices: 
[ PTDF_true, PTDF, LODF, B_full, H_full, ~ ] = Create_PTDF_LODF_B_H(Bus,Branch,RoundOffFlag);

%% Perform Contingency Ranking:
[Bompard_rank, ~] = ContingencyRanking(Bus, Branch, Load, Generator, PTDF_true);

%% Perform RTCA using DC power flows using the results of contingency ranking: 
[ Vio, flag_vio_rtca, ~] = DC_RTCA_Ranking(mpc,Bompard_rank,Rank_limit,Radial);
if (isempty(Vio)==0)
    Pre_Vio_L = Vio(:,1);
else
    Pre_Vio_L = [];
end

%% Perform feasibility test (FT) for all branches in the base-case scenario:
[ CL_Sp_vio, CutsetStack_vio, EdgeList, flag_vio_ft, ~ ] = FeasibilityTestBasecase( Flow, Capacity, A, BranchGraph );
K_rtca = [];
K_ft = [];
outage_number = 1;
Total_load_base = sum(Load(:,2));

while (continue_flag==1)    
    if ((flag_vio_rtca==1) || (flag_vio_ft==1))  
%% Display violations detected by RTCA:
        if (flag_vio_rtca==1)
            DisplayViolations_RTCA(Vio);
        end
%% Display violations detected by FT:
        if (flag_vio_ft==1)
            DisplayViolations_FT( CL_Sp_vio,CutsetStack_vio );
        end
        
%% Create inputs from RTCA for the iCA:        
        if (isempty(Vio)==0)
            if ((initial==1) || (isempty(K_rtca)==1))
                K_rtca = Vio(:,1);
            else
                K_rtca = vertcat(K_rtca,K_rtca_new);
                K_rtca = unique(K_rtca);
            end
        end        

%% Create inputs from FT for the iCA:        
        if ((initial==1) || (isempty(K_ft)==1))
                [ K_ft, Tm, Cutset_FT, ~ ] = CreateInput_ODC( CL_Sp_vio, CutsetStack_vio, Branch);          
        else            
                [ K_ft, Tm, Cutset_FT, ~] = AugmentCutsetInfo(K_ft_new,Tm_new,Cutset_FT_new,K_ft,Tm,Cutset_FT);
        end                   
        count_unique = 1;
        K_ft_unique = [];
        Tm_unique = [];
        Cutset_FT_unique = [];
        [row_K_ft, ~] = size(K_ft);
        
        for i = 1:row_K_ft
            branch_num = K_ft(i);
            flag = IsPresent(K_rtca,branch_num);
            if (flag==0)
                K_ft_unique(count_unique,:) = K_ft(i,:);
                Tm_unique(count_unique,1) = Tm(i,1);
                Cutset_FT_unique(:,:,count_unique) = Cutset_FT(:,:,i);
                count_unique = count_unique + 1;
            end
        end
        
%% Perform the Integrated Corrective Action (iCA):
        [ GeneratorNegativeChange, GeneratorPositiveChange, LoadNegativeChange, LoadPositiveChange, Branch, Load, Generator, Soln_Flag, tot_change_cost, ~ ] = IntegratedCorrectiveAction(K_rtca, PTDF, LODF, Bus, Branch, Generator, Load, Radial, K_ft_unique, Tm_unique, Cutset_FT_unique);
                
%% Update the system based upon the redispatch solution:
        if (Soln_Flag==0)
            break;
        end                
        mpc.gen(:,2) = Generator(:,2);
        for nload = 1:length(Load(:,1))
             LoadBus = Load(nload,1);
             loc = find(mpc.bus(:,1)==LoadBus);
             mpc.bus(loc,3) = Load(nload,2);
        end
        Res = rundcpf(mpc);
        Branch(:,6) = Res.branch(:,14);   
                
%% Perform RTCA following redispatch:
        [ Vio, flag_vio_rtca, ~ ] = DC_RTCA_Ranking( mpc, Bompard_rank, Rank_limit, Radial);        
        if (isempty(Vio)==0)                    
            Vio_L = Vio(:,1);
            Intersect_Pre_Vio_L = intersect(Pre_Vio_L,Vio_L);
            if (length(Intersect_Pre_Vio_L)==length(Vio_L))
                flag_vio_rtca = 0;
            else
                Pre_Vio_L = Vio_L;
            end
        end        
        if (flag_vio_rtca==1)
             K_rtca_new = Vio(:,1);
             initial = 0;
        end
        
        % Group all the injection increase and injection decrease
        % together in separate matrices
        InjectionPositiveChange = [];
        Temp = []; Temp = LoadNegativeChange; Temp(:,2) = (-1)*Temp(:,2); 
        InjectionPositiveChange = vertcat(GeneratorPositiveChange,Temp);            
        Temp = [];Temp = LoadPositiveChange; 
        InjectionNegativeChange = vertcat(GeneratorNegativeChange,Temp);
        
 %% Update the "flow" and "latent capacity graphs" following redispatch based on M-UPS algorithm:
        [Flow,Capacity,BranchFlowChange,~] = ModifiedUpdateScheme(Flow, Capacity, InjectionPositiveChange, InjectionNegativeChange, BranchGraph);
                           
%% Perform shortlisting assets following redispatch based on M-SA algorithm:
        [ShortlistedEdges, ~] = ModifiedShortlistAssets(BranchFlowChange, EdgeList, BranchGraph);        
                               
%% Perform Feasibility Test (FT) on shortlisted assets:
        [CL_Sp_vio, CutsetStack_vio, EdgeList, flag_vio_ft, ~] = FeasibilityTestOnShortlist( Flow, Capacity, A, BranchGraph, ShortlistedEdges, EdgeList );
        
        if (flag_vio_ft==1)               
%% Create inputs for Optimal Dispatch Change (ODC) for new cutsets:
            [K_ft_new, Tm_new, Cutset_FT_new, ~ ] = CreateInput_cutset( CL_Sp_vio, CutsetStack_vio, Branch);               
        end                       
                        
    else
%% Display the results of the corrective action:       
         fprintf('-----------------------------------------------------------------\n');
         fprintf('The first component has alleviated all post-contingency cut-set saturation and critical branch overloads \n');           
         fprintf('-----------------------------------------------------------------\n');
         GeneratorCost_Ar = Generator(:,6).*Generator(:,2).^2+Generator(:,5).*Generator(:,2); TotalGeneratorCost = sum(GeneratorCost_Ar);
         LoadCost_Ar = (Load(:,2)-Load(:,4)).*Load(:,3); TotalLoadCost = sum(LoadCost_Ar);            
         fprintf('Production cost = $ %f \n',TotalGeneratorCost);                 
         Net_load_shed = Total_load_base-sum(Load(:,2));
         fprintf('Total amount of load shed = %f MW \n',Net_load_shed);              
        
%% Check if there are successive branch outages in the system:
        LineOutNumber = input('\n Enter the branch number which is out (Press 0 and enter if you do not want to continue) ?');        
        fprintf('\n');
        fprintf('\n ******** New Outage: *************\n');  

        if (LineOutNumber==0)
            continue_flag = 0;
            break;
            
        else
%% Update the system matrices following the branch outage: 
            mpc.branch(LineOutNumber,11) = 0;
            Res = rundcpf(mpc);
            Branch(:,6) = Res.branch(:,14);
            Branch(:,8) = Res.branch(:,11);
            BranchGraph(:,4) = Branch(:,8);
            
            A(Branch(LineOutNumber,1),Branch(LineOutNumber,2)) = 0;
            A(Branch(LineOutNumber,2),Branch(LineOutNumber,1)) = 0;  
                                                                       
% Update the system matrices, instead of re-building 
% the matrices from scratch
            [ PTDF_true, PTDF, LODF, B_full, H_full, ~ ] = Update_PTDF_LODF_B_H( B_full, H_full, Bus, Branch, LineOutNumber,RoundOffFlag);     
            
%% Find the radial branches for the new system:
            [ Radial, ~ ] = FindRadial( Branch, A ); 
            
%% Perform contingency ranking:            
            [Bompard_rank, ~] = ContingencyRanking(Bus, Branch, Load, Generator, PTDF_true);

%% Perform RTCA using DC power flows:
            [ Vio, flag_vio_rtca, ~ ] = DC_RTCA_Ranking( mpc, Bompard_rank, Rank_limit, Radial);            
             if (isempty(Vio)==0)                    
                    Vio_L = Vio(:,1);
                    Pre_Vio_L = Vio_L;
             end                                    
             if (flag_vio_rtca==1)
                    K_rtca_new = Vio(:,1);                
             end
                         
%% A successive FT has to be performed following the branch outage
% However, the successive FT must involve the UPS algorithm and the 
% SA algorithm for fast computation
            [ Flow, Capacity, A, CL_Sp_vio, EdgeList, PathStack, EdgeSatStack, CutsetStack_vio, ~ ] = OutageAnalysis( BranchGraph, Flow, Capacity, LineOutNumber, EdgeList, A );     
            
            [ row_vio, col_vio ] = size(CL_Sp_vio);
            flag_vio_ft = 0;
            if (row_vio>=1)
                flag_vio_ft = 1;          
                initial = 0;                                 
                [ K_ft_new, Tm_new, Cutset_FT_new, ~ ] = CreateInput_cutset( CL_Sp_vio, CutsetStack_vio, Branch);                          
            end              
            
            if ((flag_vio_rtca==0) && (flag_vio_ft==0))
                fprintf('\n There are no violations detected by the Feasibility Test (FT) and DC-RTCA \n');
            end
                 
            outage_number = outage_number + 1;
        end     
    end        
end