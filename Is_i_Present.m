% This function is used by the 'subgraph' script
function [ flag ] = Is_i_Present(i,S_k)

flag = 0;
for j=1:length(S_k)
    if S_k(j)==i
        flag=1;
        break;
    end
end

end

