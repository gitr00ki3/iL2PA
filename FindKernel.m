function [ kern ] = FindKernel( dist, EPS_SQ )
%FINDKERNEL Find Kernel
    n= size(dist, 1);
    kern= zeros(n);
    for i = 1:n
        for j = i:n
            kern(i,j) = exp(- dist(i, j)^2/2*EPS_SQ);
            kern(j,i) = kern(i,j);
        end;
    end;
end