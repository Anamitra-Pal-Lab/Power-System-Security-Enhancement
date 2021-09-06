function [ Vio, flag_vio, time ] = DC_RTCA_Ranking( mpc, T_sort, Rank_limit, Radial)
   
tic;
Vio = [];
count = 1;
count_non_rad_line = 1;
overload_factor = 1;
% overload_factor = 1.1875;
overload_margin = 5;
LimitOnContingency = 10;

% while (count_non_rad_line<Bompard_rank_limit) 
    
    for i = 1:length(mpc.branch(:,1))
        nbranch = T_sort(i,1);
        flag = IsPresent(Radial(:,1),nbranch);    
        if ((flag==0) && (mpc.branch(nbranch,11)==1))
            mpc.branch(nbranch,11) = 0;
            Res = rundcpf(mpc);
            diff = [];
            diff = overload_factor*Res.branch(:,6)-abs(Res.branch(:,14));
            loc = find(diff<-(0.0001+overload_margin));
            overload = diff(loc);
            if (isempty(loc)==0)
                Vio(count,1) = nbranch;            
                Vio(count,2) = Res.branch(nbranch,1);
                Vio(count,3) = Res.branch(nbranch,2);
                diff_neg = diff(loc);
                Vio(count,4) = sum(abs(diff_neg));
                Vio(count,5) = max(abs(diff_neg));
%             NoOfOverload = length(loc);
%             overloaded_branch(count,[1:NoOfOverload]) = loc';
%             overload_magnitude(count,[1:NoOfOverload]) = overload';
                count = count+1;
            end
            mpc.branch(nbranch,11) = 1;
            count_non_rad_line = count_non_rad_line + 1;
            if (count_non_rad_line>Rank_limit)
                break;
            end
        end        
    end
% end

if (isempty(Vio)==0)
    flag_vio = 1;
    Vio = sortrows(Vio,4,'descend');
    if length(Vio(:,1))>LimitOnContingency
        Vio_new = Vio(1:LimitOnContingency,:);
        Vio = []; Vio = Vio_new;
    end    
else
    flag_vio = 0;
end
time = toc;
 
end

