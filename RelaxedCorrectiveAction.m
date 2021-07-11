function [ GeneratorNegativeChange, GeneratorPositiveChange, LoadNegativeChange, LoadPositiveChange, Branch, Load, Generator, Soln_Flag, Flow_dc, Rate_dc, tot_change_cost, time ] = RelaxedCorrectiveAction( K, Tm, Cutset_FT, PTDF, Bus, Branch, Generator, Load)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program solves the optimization 
% problem for the relaxed corrective action (rCA) used 
% in the second component
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    

    tic;        
    % Find the length of different arrays related to the cut-set violations
    NumOfCritCutset = length(K(:,1));
    MaxNumBranchCritCutset = length(K(1,:));
    [numrow_K, numcol_K, numsheet] = size(Cutset_FT);
    
    %% Initialize flow across different cutsets:
    flow_cutset = [];
    for snum = 1:numsheet
        flow_total = 0;
        for rnum = 1:numrow_K
            if Cutset_FT(rnum,1,snum)~=0  
                FromBus = Cutset_FT(rnum,1,snum); 
                ToBus = Cutset_FT(rnum,2,snum);           
                for i = 1:length(Branch(:,1)) 
                    if (FromBus==Branch(i,1) && ToBus==Branch(i,2))
                        flow_val = Branch(i,6);
                        flow_total = flow_total+flow_val;
                    elseif (FromBus==Branch(i,2) && ToBus==Branch(i,1))
                        flow_val = (-1)*Branch(i,6);
                        flow_total = flow_total+flow_val;                   
                    end
                end           
            end
            flow_cutset(snum,1) = flow_total;
        end
    end  
 %% Set-up the objective function:     
    b = Generator(:,5); % The linear cost coefficient
    c = Generator(:,6); % The quadratic cost coefficient    
    Pg_old_ar = Generator(:,2); % Old power generation    
    f_gen_lin = (2*(Pg_old_ar.*c) + b);     
    f_load = Load(:,3); 
    f = vertcat(f_gen_lin,f_load); 
    
    noofline = length(Branch(:,1));
    noofgen = length(Generator(:,1));
    noofload = length(Load(:,1));
    
%% Constraints for the conservation of energy:
    Pivot = 1;
    count_Sa = 1;
    for i = 1:noofgen
       Xa(count_Sa,1) = Pivot;Ya(count_Sa,1) = i;Va(count_Sa,1) = 1;
       count_Sa = count_Sa + 1;
    end
    for i = noofgen+1:(noofgen+noofload)
       Xa(count_Sa,1) = Pivot;Ya(count_Sa,1) = i;Va(count_Sa,1) = -1;
       count_Sa = count_Sa + 1;
    end
    Rhs_conserve = [0];
    Sign_conserve = [ '=' ];
    
     %% Constraints for the injection limits:
    Constraint_pinj = eye(noofgen+noofload,noofgen+noofload);
    count_Sb = 1;
    for i = 1:(noofgen+noofload)         
         Xb(count_Sb,1) = Pivot+i;Yb(count_Sb,1) = i;Vb(count_Sb,1) = 1;
        count_Sb = count_Sb+1;
    end
    Pivot = Pivot+noofgen+noofload;
    X_val = Pivot+1:Pivot+noofgen+noofload; X_val = X_val'; 
    Xb = vertcat(Xb,X_val); 
    Yb = repmat(Yb,2,1); 
    Vb = repmat(Vb,2,1);
    Pivot = Pivot+noofgen+noofload;       
    for ngen = 1:noofgen
        GenBusNum = Generator(ngen,1);
        Pgen_old = Generator(ngen,2);
        Pgen_max = Generator(ngen,3);
        Pgen_min = Generator(ngen,4);
        Rhs_pinj_max(ngen,1) = Pgen_max-Pgen_old;               
        Rhs_pinj_min(ngen,1) = Pgen_min-Pgen_old;        
        Sign_pinj_max(ngen,1) = '<';
        Sign_pinj_min(ngen,1) = '>';                        
    end
    
    %% LHS and RHS for the injection limits:
    for nload = 1:noofload
        LoadBusNum = Load(nload,1);
        Pload_old = Load(nload,2);        
        Rhs_pinj_max(noofgen+nload,1) = 0;
        Rhs_pinj_min(noofgen+nload,1) = -Pload_old;
        Sign_pinj_max(noofgen+nload,1) = '<';
        Sign_pinj_min(noofgen+nload,1) = '>';
    end
    
    %% Constraints for pre-contingency power flow in each branch:    
    Constraint_flow = zeros(noofline,noofgen+noofload);
    count_Sc = 1;
    noofflow_cstr = 0;
    for nline = 1:noofline
        flag_flow = 0;
        if (Branch(nline,8)==1)
            for ngen = 1:noofgen
                GenBusNum = Generator(ngen,1);
                if (GenBusNum==length(Bus))
                    PTDF_val = 0;
                else        
                    PTDF_val = PTDF(nline,GenBusNum);
                end
                if abs(PTDF_val)>10^-5
                    Constraint_flow(nline,ngen) = PTDF_val; 
                    Xc(count_Sc,1) = Pivot+noofflow_cstr+1; Yc(count_Sc,1) = ngen; Vc(count_Sc,1) = PTDF_val;
                    count_Sc = count_Sc+1;
                    flag_flow = 1;
                end
            end
            for nload = 1:noofload
                LoadBusNum = Load(nload,1);
                if (LoadBusNum==length(Bus))
                    PTDF_val = 0;
                else
                    PTDF_val = PTDF(nline,LoadBusNum);
                end
%                 if (PTDF_val~=0)
                if abs(PTDF_val)>10^-5
                    Constraint_flow(nline,noofgen+nload) = (-1)*PTDF_val;
                    Xc(count_Sc,1) = Pivot+noofflow_cstr+1; Yc(count_Sc,1) = noofgen+nload; Vc(count_Sc,1) = (-1)*PTDF_val;
                    count_Sc = count_Sc+1;
                    flag_flow = 1;
                end
            end     
            if (flag_flow==1)
                flow_old = Branch(nline,6);
                flow_max = Branch(nline,7);
                flow_min = (-1)*Branch(nline,7);
                Rhs_MaxFlow(noofflow_cstr+1,1) = flow_max-flow_old; 
                Rhs_MinFlow(noofflow_cstr+1,1) = flow_min-flow_old; 
                Sign_Maxflow(noofflow_cstr+1,1) = '<';
                Sign_Minflow(noofflow_cstr+1,1) = '>';
                noofflow_cstr = noofflow_cstr + 1;
            end
        end
    end
    Pivot = Pivot + noofflow_cstr;
    X_val = Xc+noofflow_cstr*ones(length(Xc),1);
    Xc = vertcat(Xc,X_val);
    Yc = repmat(Yc,2,1);
    Vc = repmat(Vc,2,1);
    Pivot = Pivot + noofflow_cstr;
    
    %% Constraints for cut-set power transfer limit:
      %% Constraints for the power transfer across the cut-set:
    [row_K, col_K] = size(K);   
    count_Sd = 1;
    for ncutset = 1:row_K        
        for ngen = 1:noofgen
            GenBusNum = Generator(ngen,1);   
            PTDF_cutset = 0;
            for nbranch = 1:col_K
                if K(ncutset,nbranch)~=0
                    BranchNum = K(ncutset,nbranch);
                    % Check if the direction of the branch is same the direction
                    % of the cut-set.
                    F_Branch = Branch(BranchNum,1);
                    T_Branch = Branch(BranchNum,2);
                    Sign = 0;
                    if IsPresent(Cutset_FT(:,1,ncutset),F_Branch)==1 && IsPresent(Cutset_FT(:,2,ncutset),T_Branch)==1
                        Sign = 1;
                    elseif IsPresent(Cutset_FT(:,2,ncutset),F_Branch)==1 && IsPresent(Cutset_FT(:,1,ncutset),T_Branch)==1
                        Sign = -1;
                    else
                        Sign = 0;
                    end
                    if (GenBusNum < length(Bus(:,1)))
                        PTDF_val = Sign*PTDF(BranchNum,GenBusNum);
                    else
                        PTDF_val = 0;
                    end
                    PTDF_cutset = PTDF_cutset+PTDF_val; 
                end
            end
           if abs(PTDF_cutset)>10^-5
               Xd(count_Sd,1) = Pivot+ncutset;
               Yd(count_Sd,1) = ngen;
               Vd(count_Sd,1) = PTDF_cutset;
               count_Sd = count_Sd + 1;
           end
        end
        for nload = 1:noofload
            LoadBusNum = Load(nload,1);   
            PTDF_cutset = 0;
            for nbranch = 1:col_K
                if K(ncutset,nbranch)~=0
                    BranchNum = K(ncutset,nbranch);
                    % Check if the direction of the branch is same the direction
                    % of the cut-set.
                    F_Branch = Branch(BranchNum,1);
                    T_Branch = Branch(BranchNum,2);
                    Sign = 0;
                    if IsPresent(Cutset_FT(:,1,ncutset),F_Branch)==1 && IsPresent(Cutset_FT(:,2,ncutset),T_Branch)==1
                        Sign = 1;
                    elseif IsPresent(Cutset_FT(:,2,ncutset),F_Branch)==1 && IsPresent(Cutset_FT(:,1,ncutset),T_Branch)==1
                        Sign = -1;
                    else
                        Sign = 0;
                    end  
                    if (LoadBusNum < length(Bus(:,1)))
                        PTDF_val = Sign*PTDF(BranchNum,LoadBusNum);
                    else
                        PTDF_val = 0;
                    end                        
                    PTDF_cutset = PTDF_cutset+(-1)*PTDF_val; 
                end
            end
           if abs(PTDF_cutset)>10^-5
               Xd(count_Sd,1) = Pivot+ncutset;
               Yd(count_Sd,1) = noofgen+nload;
               Vd(count_Sd,1) = PTDF_cutset;
               count_Sd = count_Sd + 1;
           end
        end     
        
        Tot_rate = 0;
        Tot_flow = 0;
        for i = 1:length(K(ncutset,:))
            if K(ncutset,i)==0
                break;
            end                
            if (i>1)
                Tot_rate = Tot_rate+Branch(K(ncutset,i),7);
            end
            A_Branch = Branch(K(ncutset,i),1);
            B_Branch = Branch(K(ncutset,i),2);
            if (IsPresent(Cutset_FT(:,1,ncutset),A_Branch)==1) && (IsPresent(Cutset_FT(:,2,ncutset),B_Branch)==1)
                Tot_flow = Tot_flow + Branch(K(ncutset,i),6);
            else
                Tot_flow = Tot_flow + (-1)*Branch(K(ncutset,i),6);
            end                                           
        end  
        Rhs_cutset(ncutset,1) = Tot_rate-Tot_flow;
        Sign_cutset(ncutset,1) = '<';
    end
    Pivot = Pivot + row_K;
    
    %% Concatenate all constraints:
    X = vertcat(Xa,Xb,Xc,Xd);
    Y = vertcat(Ya,Yb,Yc,Yd);
    V = vertcat(Va,Vb,Vc,Vd);    
    Constraint_SP = sparse(X,Y,V);
    
    RHS = vertcat(Rhs_conserve,Rhs_pinj_max,Rhs_pinj_min,Rhs_MaxFlow,Rhs_MinFlow,Rhs_cutset);
    SIGN = vertcat(Sign_conserve,Sign_pinj_max,Sign_pinj_min,Sign_Maxflow,Sign_Minflow,Sign_cutset);
    
     %% Use the quadratic cost coefficients:
    f_quad_gen = zeros(noofgen+noofload);
    for i = 1:noofgen
        c_quad = c(i,1);    
        f_quad_gen(i,i) = c_quad;
    end    
    
    %% Set the model parameters:
    model.obj = f;
    model.Q = sparse(f_quad_gen); 
    model.A = sparse(Constraint_SP); 
    model.sense = SIGN;
    model.rhs = RHS;
    model.lb = Rhs_pinj_min; 

    clear params;
    params.outputflag = 0;
    result = gurobi(model, params);
    
    if length(result.status)==7
        Soln_Flag = 1;
        xf = result.x;
        %% Compute all measurement values after solving the optmization:
        % New branch flows:
        flow_old = zeros(length(Branch(:,1)),4);
        flow_new = zeros(length(Branch(:,1)),4);
        flow_old(:,1) = Branch(:,1);
        flow_old(:,2) = Branch(:,2);
        flow_old(:,3) = Branch(:,6);
        flow_old(:,4) = Branch(:,7);
        delta_flow = Constraint_flow*xf;
        flow_new(:,1) = Branch(:,1);
        flow_new(:,2) = Branch(:,2);
        flow_new(:,3) = flow_old(:,3)+delta_flow;
        flow_new(:,4) = Branch(:,7);
        % New dispatch:
        gen_old = zeros(length(Generator(:,1)),2);
        gen_new = zeros(length(Generator(:,1)),2);
        gen_old(:,1) = Generator(:,1);
        gen_old(:,2) = Generator(:,2);
        load_old(:,1) = Load(:,1);
        load_old(:,2) = Load(:,2);        
        delta_inj = Constraint_pinj*xf;        
        delta_pgen = delta_inj([1:noofgen],1);
        delta_pload = delta_inj([noofgen+1:noofgen+noofload],1);        
        Generator_New(:,1) = Generator(:,1);
        Generator_New(:,2) = gen_old(:,2)+delta_pgen;
        Load_New(:,1) = Load(:,1);
        Load_New(:,2) = load_old(:,2)+delta_pload;
        
        %% Finding the actual cost using quadratic and linear cost coefficients:
        cost_linear = transpose(f);
        cost_quad = horzcat(transpose(c),zeros(1,noofload));
        tot_change_cost = cost_linear*xf + cost_quad*(xf.^2);
        
        %% Find the positions where non-zero changes have occurred in Pgen:
        [ indpos, ~ ] = find(delta_pgen>0.001);
        [ indneg, ~ ] = find(delta_pgen<-0.001);
        GeneratorPositiveChange(:,1) = Generator(indpos,1);
        GeneratorPositiveChange(:,2) = delta_pgen(indpos);

        GeneratorNegativeChange(:,1) = Generator(indneg,1);
        GeneratorNegativeChange(:,2) = delta_pgen(indneg);
        
        %% Find the positions where non-zero changes have occurred in Pload:
        [ indpos, ~ ] = find(delta_pload>0.001);
        [ indneg, ~ ] = find(delta_pload<-0.001);
        LoadPositiveChange(:,1) = Load(indpos,1);
        LoadPositiveChange(:,2) = delta_pload(indpos);
        
        LoadNegativeChange(:,1) = Load(indneg,1);
        LoadNegativeChange(:,2) = delta_pload(indneg);
        
        %% Get the data for the next stage:        
        Branch(:,6) = flow_new(:,3);        
        Generator(:,2) = Generator_New(:,2);        
        Load(:,2) = Load_New(:,2);
        
        else
        Soln_Flag = 0;
        GeneratorNegativeChange = [];    
        GeneratorPositiveChange = [];
        LoadNegativeChange = [];
        LoadPositiveChange = [];
        tot_change_cost = 0;
    end
       
    %% Print the change in dispatches on the screen:
   if (Soln_Flag==1) 
        fprintf('-------------------------------------------- \n');
        fprintf('Total decrease in load = %f \n',round(sum(LoadNegativeChange(:,2))));
        fprintf('Total increase in load = %f \n',round(sum(LoadPositiveChange(:,2))));        
        fprintf('Total decrease in dispatch  = %f \n', round(sum(GeneratorNegativeChange(:,2))));
        fprintf('Total increase in dispatch  = %f \n', round(sum(GeneratorPositiveChange(:,2))));         
        fprintf('Total change in cost = $ %f \n',round(tot_change_cost));
        fprintf('-------------------------------------------- \n');
   else
        fprintf('-------------------------------------------- \n');
        fprintf('No feasible solution obtained! \n');
        fprintf('-------------------------------------------- \n');
   end
   
   %% Obtain the dc power flow graph:
   NoOfBus = length(Bus);   
   Flow_dc = sparse(NoOfBus,NoOfBus);
   Rate_dc = sparse(NoOfBus,NoOfBus);
   for i = 1:length(Branch(:,1))       
       Flow_dc(Branch(i,1),Branch(i,2)) = Branch(i,6);
       Flow_dc(Branch(i,2),Branch(i,1)) = (-1)*Branch(i,6);
       Rate_dc(Branch(i,1),Branch(i,2)) = Branch(i,7);
   end
   time = toc;
end