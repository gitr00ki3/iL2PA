function [ lap ] = FindLaplacian( nodes, RANGE, SIG_SQ, C )
%FINDLAPLACIAN Finds Laplacian of a graph L=D-W
    n= size(nodes, 1);
    b= zeros(n, 1);
    w= zeros(n);

    for i = 1:n
        count = 0;
        sumD = 0;
        for j = 1:n
            if (i ~= j)
                d = FindEuclidean(nodes(i, 1), nodes(j, 1), ...
                    nodes(i, 2), nodes(j, 2));
                if (d <= RANGE)
                    sumD = sumD + d;
                    count = count + 1;
                end;
            end;
        end;
        if count ~= 0
            b(i,1) = sumD/count;
        end;
    end;
    
    for i = 1:n
        for j = i+1:n
            if(i ~= j)
                w(i, j) = (1/C)*(exp(- ( FindEuclidean(nodes(i, 1), nodes(j, 1), ...
                    nodes(i, 2), nodes(j, 2))^2) / SIG_SQ ) ...
                    * exp(- (b(i)-b(j))^2 / SIG_SQ ));
            end;
        end;
    end;
    
    w= triu(w)+triu(w, 1)';
    d= diag(sum(w, 2));
    lap= d-w;
    lap= triu(lap)+triu(lap, 1)';
end


