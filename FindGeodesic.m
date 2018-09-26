function [ g ] = FindGeodesic( dist )
%FINDGEODESIC Finds Geodesic on input distance matrix
    n= size(dist, 1);
    g = zeros(n,n);

    for i = 1:n
        for j = i+1:n
            if(dist(i,j) > 0)
                g(i,j) = dist(i,j);
                g(j,i) = dist(i,j);
            else
                g(i,j) = inf(1);
                g(j,i) = inf(1);
            end;
        end;
    end;

    for k = 1:n
        for i = 1:n
            for j = i:n
                g(i,j) = min(g(i,j), g(i,k)+g(k,j));
                g(j,i) = g(i,j);
            end;
        end;
    end;
end