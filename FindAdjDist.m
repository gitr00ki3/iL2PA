function [ dist ] = FindAdjDist( nodes, adj )
%FINDADJDIST Finds distance based on adjacency of the graph
    n= size(nodes, 1);
    dist= zeros(n, n);
    [row, col]= find(adj);
    rowN= size(row, 1);
    for i= 1:rowN
        dist(row(i), col(i))= FindEuclidean(nodes(row(i), 1), nodes(col(i), 1), ...
            nodes(row(i), 2), nodes(col(i), 2));
    end;
end

