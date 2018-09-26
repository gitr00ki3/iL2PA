function [ euc_dist ] = FindEuclidean( x1, x2, y1, y2 )
%FINDEUCLIDEAN Finds Euclidean distance between two points
    euc_dist= sqrt((x1-x2)^2+(y1-y2)^2);
end

