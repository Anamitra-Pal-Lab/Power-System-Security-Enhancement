function [ K, Tm, Cutset_FT, time ] = CreateInput_ODC( CL_Sp_vio, CutsetStack_vio,Branch) 
    
    tic;
    [noofcutset,col] = size(CL_Sp_vio);
    [xdim,ydim,zdim] = size(CutsetStack_vio);
    ncutset = 1;
    K = []; Tm = 0; Cutset_FT = [];    
    for r = 1:noofcutset        
%         Tm(ncutset,1) = floor(CL_Sp_vio(r,4));
        Tm(ncutset,1) = CL_Sp_vio(r,4);
        for i = 1:xdim
            F = CutsetStack_vio(i,1,r);
            T = CutsetStack_vio(i,2,r);
            if F~=0
               Cutset_FT(i,1,ncutset) = F;Cutset_FT(i,2,ncutset) = T;
               for j = 1:length(Branch(:,1))
                   if (F==Branch(j,1) && T==Branch(j,2)) || (F==Branch(j,2) && T==Branch(j,1))
                       K(ncutset,i) = j;
                       break;
                   end
               end               
            end
        end  
        ncutset = ncutset + 1;
    end 
    time = toc;   
    
end


