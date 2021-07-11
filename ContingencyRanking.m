function [ T_sort, timeBomp ] = ContingencyRanking( Bus, Branch, Load, Generator, PTDF )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program Description: This program finds the contingency ranking
% based upon the PTDFs and branch ratings
%
% Author: Reetam Sen Biswas 
% Arizona State University
% 
% Last Modified: 03/20/2020; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

tic;
NoOfGen = length(Generator(:,1));
NoOfLoad = length(Load(:,1));
NoOfBranch = length(Branch(:,1));
NoOfBus = length(Bus(:,1));
C_gd = zeros(NoOfGen,NoOfLoad); 
Pl_max = Branch(:,7);
zero_col = zeros(NoOfBranch,1);
PTDF = horzcat(PTDF,zero_col);

for g = 1:NoOfGen
    gbus = Generator(g,1);
    for d = 1:NoOfLoad
        dbus = Load(d,1);
        if (gbus~=dbus)
            
            if (gbus~=NoOfBus)
                ptdf_lines_gbus = PTDF(:,gbus);        
            else
                ptdf_lines_gbus = zeros(NoOfBranch,1); 
            end
            
            if (dbus~=NoOfBus)
                ptdf_lines_dbus = PTDF(:,dbus);        
            else
                ptdf_lines_dbus = zeros(NoOfBranch,1); 
            end        
                
            ptdf_lines_gbus_dbus = ptdf_lines_gbus - ptdf_lines_dbus;
        
            value_ar = Pl_max./abs(ptdf_lines_gbus_dbus);
        
            value = min(value_ar);
        
            C_gd(g,d) = value;                        
        
        end
    end
end       

for nline = 1:NoOfBranch   
    
     gen_buses = Generator(:,1);
     load_buses = Load(:,1);
     
     ptdf_gen = PTDF(nline,gen_buses)';
     ptdf_load = PTDF(nline,load_buses);
     
     ptdf_gen_mat = repmat(ptdf_gen,[1,NoOfLoad]);
                     
     ptdf_load_mat = repmat(ptdf_load,[NoOfGen,1]);
     
     ptdf_gen_load_mat = ptdf_gen_mat - ptdf_load_mat;
     
     common = intersect(gen_buses,load_buses');
     
     for j = 1:length(common)
        gen_loc = find(gen_buses==common(j));
        load_loc = find(load_buses==common(j));
        ptdf_gen_load_mat(gen_loc,load_loc) = 0;
     end
     
     [r,c] = find(ptdf_gen_load_mat<0);
     ptdf_gen_load_mat_p = ptdf_gen_load_mat;
     
     for k = 1:length(r)
         ptdf_gen_load_mat_p(r(k),c(k)) = 0;
     end            
         
     ptdf_weight_positive = ptdf_gen_load_mat_p.*C_gd;
     Tp = sum(sum(ptdf_weight_positive));
     
     [r,c] = find(ptdf_gen_load_mat>0);
     ptdf_gen_load_mat_n = ptdf_gen_load_mat;
     
     for k = 1:length(r)
         ptdf_gen_load_mat_n(r(k),c(k)) = 0;
     end            
     
     ptdf_weight_negative = ptdf_gen_load_mat_n.*C_gd;
     Tn = sum(sum(ptdf_weight_negative));
     
    T(nline,1) = nline;
    T(nline,2) = Branch(nline,1);
    T(nline,3) = Branch(nline,2);
    T(nline,4) = max(Tp,abs(Tn));
     
end
T_max = max(T(:,4));
T(:,4) = T(:,4)./T_max;
T_sort = sortrows(T,4,'descend');
timeBomp = toc;

end

