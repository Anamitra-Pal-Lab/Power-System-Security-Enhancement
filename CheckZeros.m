function [ FlagVec ] = CheckZeros( Vec )
    FlagVec = [];
    for i = 1:length(Vec)
        if Vec(i)<10^-3 && Vec(i)>(-1)*10^3
            FlagVec(i,1) = 1;
        else
            FlagVec(i,1) = 0;
        end
    end

end

