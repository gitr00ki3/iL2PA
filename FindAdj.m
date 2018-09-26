function [ adj, dist ] = FindAdj( nodes, range )
%FINDADJ Finds adjacency of a graph[nodes(nx2)] using the given range on calculated
%euclidean distance
    n= size(nodes, 1);      % Get number of nodes from the graph
    adj= zeros(n, n);
    dist= zeros(n, n);
    for i= 1:n
        for j= i+1:n
            d= FindEuclidean(nodes(i, 1), nodes(j, 1), nodes(i, 2), nodes(j, 2));
            if d<=range
                adj(i, j)= 1;
                dist(i, j)= d;
            end;
        end;
    end;
    adj= triu(adj)+triu(adj,1)';    % Copy upper triangular to lower
    dist= triu(dist)+triu(dist,1)';
end

