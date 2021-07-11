function [ GeneratorNegativeChange, GeneratorPositiveChange, LoadNegativeChange, LoadPositiveChange, Branch, Load, Generator, Soln_Flag, tot_change_cost, time ] = IntegratedCorrectiveAction(K_rtca, PTDF, LODF, Bus, Branch, Generator, Load, RadialLines,K_ft_unique, Tm, Cutset_FT )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program solves the optimization 
% problem for the integrated corrective action (iCA), 
% used in the proposed first component
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
tic;    

%% Set-up the objective function:
    b = Generator(:,5); 
    c = Generator(:,6);     
    Pg_old_ar = Generator(:,2); 
    f_gen_lin = (2*(Pg_old_ar.*c) + b);     
    f_load = Load(:,3); 
    f = vertcat(f_gen_lin,f_load);  
    
    noofline = length(Branch(:,1));
    noofgen = length(Generator(:,1));
    noofload = length(Load(:,1));   
           
    row_K = size(K_rtca,1);
    ContingencySet = [];
    count = 1;
    for lnum = 1:row_K  
        if (IsPresent(RadialLines,K_rtca(lnum,1))~=1) && (Branch(K_rtca(lnum,1),8)==1)
            ContingencySet(count,1) = K_rtca(lnum,1);             
            count = count+1;
        end
    end
    noofconting = size(ContingencySet,1);
    
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
    
    %% Constraints for the power injection limits:
    % Constraints for the injection maximum limit:
    Constraint_pinj = eye(noofgen+noofload,noofgen+noofload);
    count_Sb = 1;
    for i = 1:(noofgen+noofload)         
         Xb(count_Sb,1) = Pivot+i;Yb(count_Sb,1) = i;Vb(count_Sb,1) = 1;
        count_Sb = count_Sb+1;
    end
    Pivot = Pivot+noofgen+noofload;
    % Constraints for the injection minimum limit:
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
                if (PTDF_val~=0)
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
                if (PTDF_val~=0)
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
    
%% Constraints for post-contingency branch flows:     
    count_post = 1;
    count_Sd = 1; 
    noofpostconting_cstr = 0;
    Xd = []; Yd = []; Vd = [];Rhs_MaxFlow_post = [];Rhs_MinFlow_post = [];
    Sign_Maxflow_post = [];Sign_Minflow_post = []; 
    
    for Cline = 1:size(ContingencySet,1)
        l = ContingencySet(Cline,1);
        flow_old_l = Branch(l,6);
        for k = 1:noofline 
            flag_flow_post = 0;             
             if (Branch(k,8)==1)                   
                LODF_k_l = LODF(k,l);
                for ngen = 1:noofgen
                    GenBusNum = Generator(ngen,1);
                    if (GenBusNum==length(Bus))
                        PTDF_k = 0;
                        PTDF_l = 0; 
                    else                    
                        PTDF_k = PTDF(k,GenBusNum);
                        PTDF_l = PTDF(l,GenBusNum);
                    end
                    Value = PTDF_k+PTDF_l*LODF_k_l;
                    if (Value~=0)
                        Xd(count_Sd,1) = Pivot+count_post; Yd(count_Sd,1) = ngen; Vd(count_Sd,1) = Value;
                        count_Sd = count_Sd + 1;     
                        flag_flow_post = 1;
                    end
                end
                for nload = 1:noofload
                    LoadBusNum = Load(nload,1);
                    if (LoadBusNum==length(Bus))
                        PTDF_k = 0;
                        PTDF_l = 0;
                    else
                        PTDF_k = PTDF(k,LoadBusNum);
                        PTDF_l = PTDF(l,LoadBusNum);
                    end
                    Value = (-1)*PTDF_k + (-1)*PTDF_l*LODF_k_l;
                    if (Value~=0)
                        Xd(count_Sd,1) = Pivot+count_post; Yd(count_Sd,1) = noofgen+nload; Vd(count_Sd,1) = Value;
                        count_Sd = count_Sd + 1; 
                        flag_flow_post = 1;
                    end
                end
                if (flag_flow_post==1)
                    flow_old_k = Branch(k,6);
                    flow_max_k = Branch(k,7);
                    flow_min_k = (-1)*Branch(k,7);                
                    Rhs_MaxFlow_post(noofpostconting_cstr+1,1) = flow_max_k-(flow_old_k+flow_old_l*LODF_k_l); 
                    Rhs_MinFlow_post(noofpostconting_cstr+1,1) = flow_min_k-(flow_old_k+flow_old_l*LODF_k_l); 
                    Sign_Maxflow_post(noofpostconting_cstr+1,1) = '<';
                    Sign_Minflow_post(noofpostconting_cstr+1,1) = '>';
                    noofpostconting_cstr = noofpostconting_cstr+1;                                                               
                    count_post = count_post + 1;
                end
            end
        end
    end            
    
    Pivot = Pivot+noofpostconting_cstr;
    X_val = Xd+noofpostconting_cstr*ones(length(Xd),1);
    
    Xd = vertcat(Xd,X_val);
    Yd = repmat(Yd,2,1);
    Vd = repmat(Vd,2,1);    
    Pivot = Pivot+noofpostconting_cstr;  
    
%% Constraints for cutset power trasnfer:    
    [row_K, col_K] = size(K_ft_unique);   
    count_Se = 1;
    Xe = [];
    Ye = [];
    Ve = [];
    Rhs_cutset = [];
    Sign_cutset = [];
    
    for ncutset = 1:row_K        
        for ngen = 1:noofgen
            GenBusNum = Generator(ngen,1);   
            PTDF_cutset = 0;
            for nbranch = 1:col_K
                if K_ft_unique(ncutset,nbranch)~=0
                    BranchNum = K_ft_unique(ncutset,nbranch);
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
           if (PTDF_cutset~=0)
               Xe(count_Se,1) = Pivot+ncutset;
               Ye(count_Se,1) = ngen;
               Ve(count_Se,1) = PTDF_cutset;
               count_Se = count_Se + 1;
           end
        end
        for nload = 1:noofload
            LoadBusNum = Load(nload,1);   
            PTDF_cutset = 0;
            for nbranch = 1:col_K
                if K_ft_unique(ncutset,nbranch)~=0
                    BranchNum = K_ft_unique(ncutset,nbranch);            
                    F_Branch = Branch(BranchNum,1);
                    T_Branch = Branch(BranchNum,2);
                    Sign = 0;
                    if (IsPresent(Cutset_FT(:,1,ncutset),F_Branch)==1 && IsPresent(Cutset_FT(:,2,ncutset),T_Branch)==1)
                        Sign = 1;
                    elseif (IsPresent(Cutset_FT(:,2,ncutset),F_Branch)==1 && IsPresent(Cutset_FT(:,1,ncutset),T_Branch)==1)
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
           if (PTDF_cutset~=0)
               Xe(count_Se,1) = Pivot+ncutset;
               Ye(count_Se,1) = noofgen+nload;
               Ve(count_Se,1) = PTDF_cutset;
               count_Se = count_Se + 1;
           end
        end
        
        Tot_rate = 0;
        Tot_flow = 0;
        for i = 1:length(K_ft_unique(ncutset,:))
            if K_ft_unique(ncutset,i)==0
                break;
            end                
            if (i>1)
                Tot_rate = Tot_rate+Branch(K_ft_unique(ncutset,i),7);
            end
            A_Branch = Branch(K_ft_unique(ncutset,i),1);
            B_Branch = Branch(K_ft_unique(ncutset,i),2);
            if (IsPresent(Cutset_FT(:,1,ncutset),A_Branch)==1) && (IsPresent(Cutset_FT(:,2,ncutset),B_Branch)==1)
                Tot_flow = Tot_flow + Branch(K_ft_unique(ncutset,i),6);
            else
                Tot_flow = Tot_flow + (-1)*Branch(K_ft_unique(ncutset,i),6);
            end                                           
        end  
        Rhs_cutset(ncutset,1) = Tot_rate-Tot_flow;
        Sign_cutset(ncutset,1) = '<';
    end
    Pivot = Pivot + row_K;
   
    X = vertcat(Xa,Xb,Xc,Xd,Xe);
    Y = vertcat(Ya,Yb,Yc,Yd,Ye);
    V = vertcat(Va,Vb,Vc,Vd,Ve);
    T = horzcat(X,Y,V);
    
    Constraint_SP = sparse(X,Y,V);
    
    %% Combine all the Constraint Matrices Together:
    RHS = vertcat(Rhs_conserve,Rhs_pinj_max,Rhs_pinj_min,Rhs_MaxFlow,Rhs_MinFlow,Rhs_MaxFlow_post,Rhs_MinFlow_post,Rhs_cutset);
    SIGN = vertcat(Sign_conserve,Sign_pinj_max,Sign_pinj_min,Sign_Maxflow,Sign_Minflow,Sign_Maxflow_post,Sign_Minflow_post,Sign_cutset);
    
    %% Use the quadratic cost coefficients:
    f_quad_gen = zeros(noofgen+noofload);
    for i = 1:noofgen
        c_quad = c(i,1);    
        f_quad_gen(i,i) = c_quad; % Additional soft constraint on delta_Pgi
    end    
    
    %% Set the model parameters:
    model.obj = f;
    model.Q = sparse(f_quad_gen); % Include quadratic cost coefficients
    model.A = Constraint_SP; 
    model.sense = SIGN;
    model.rhs = RHS;
    model.lb = Rhs_pinj_min; 

    clear params;
    params.outputflag = 0;
    result = gurobi(model, params);
    
    if ((strcmp(result.status,'OPTIMAL')==1) || (strcmp(result.status,'SUBOPTIMAL')==1))
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
        % Get the new flows for the branch:
        Branch(:,6) = flow_new(:,3);
        % Get the new generation values:
        Generator(:,2) = Generator_New(:,2);
        % Get the new load values:
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
        fprintf('Total amount of load shed = %f \n',round(sum(LoadNegativeChange(:,2))));
        fprintf('Total increase in dispatch  = %f \n', round(sum(GeneratorPositiveChange(:,2)))); 
        fprintf('Total decrease in dispatch  = %f \n', round(sum(GeneratorNegativeChange(:,2))));
        fprintf('Total change in cost of generation = $ %f \n',round(tot_change_cost));
        fprintf('-------------------------------------------- \n');
   else
        fprintf('-------------------------------------------- \n');
        fprintf('No feasible solution obtained! \n');
        fprintf('-------------------------------------------- \n');
   end
   time = toc;  

end

