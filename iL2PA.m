%% Clear the previous run
close all;clear;clc;
load('0.3979MainReRun.mat');
rms_index= 0;
[adj, distO]= FindAdj(nodes, const.RANGE);        % Find the adjacency and distance
%% RSSI Trilateration based localization
noisy_nodes(1:const.AN, :)= nodes(1:const.AN, :);

figure('Name', 'RSSI Position calculation');
plot(nodes(1:const.AN, 1), nodes(1:const.AN, 2), '^', ...
    nodes(const.AN+1:const.N, 1), nodes(const.AN+1:const.N, 2), 'o', ...
    noisy_nodes(const.AN+1:const.N, 1), noisy_nodes(const.AN+1:const.N, 2), '+');%, ...
%     nodes(nodes_wo_an(:, 1), 1), nodes(nodes_wo_an(:, 1), 2), '*r');
hold on;
lineX= [nodes(const.AN+1:const.N, 1) noisy_nodes(const.AN+1:const.N, 1)];
lineY= [nodes(const.AN+1:const.N, 2) noisy_nodes(const.AN+1:const.N, 2)];
for i= 1:const.N-const.AN
    line(lineX(i, :), lineY(i, :), 'LineStyle', '-', 'Color', 'k');
    hold on;
end;
rectangle('Position', [const.MIN_D, const.MIN_D, ...
    const.MAX_D-const.MIN_D, const.MAX_D-const.MIN_D], 'EdgeColor', 'g');
legend('Anchor', 'Non anchor', 'Noisy');
rms_index= rms_index+1;
rms_error(rms_index)= myrms(nodes, noisy_nodes, const.N);
str = ['RMS Error= ', num2str(rms_error(rms_index))];
title(str);
%% Calculate Geodesic matrix on distN
% distN= distO;
distN= FindAdjDist(noisy_nodes, adj);
for i= 1:const.N
    for j= i+1:const.N
        if distN(i, j)~= 0
            distN(i, j)= exp(- distN(i, j)^2/2*const.EPS_SQ);
        end;
    end;
end;
distN= triu(distN)+triu(distN, 1)';
g= FindGeodesic(distN);

%% Calculate kernel
kern= FindKernel(g, const.EPS_SQ);
kernL= kern(1:const.N, 1:const.AN);
an_nodes= nodes(1:const.AN, :);
f= noisy_nodes;
%% iL2PA
pre_error= inf;
now_error= rms_error(rms_index);
while now_error<pre_error
    lap= FindLaplacian(f, const.RANGE, const.SIG_SQ, const.C);
    t1= kernL*transpose(kernL);
    t2= const.LAMBDA*kern;
    t3= (const.GAMMA*kern)*(lap*kern);
    t4= kernL*an_nodes;
    % Calculating alpha
    a = pinv(t1 + t2 + t3)*t4;
    % Calculating f and check if more regression needed
    f = kern * a;
    % Procrustes analysis
    f_an= f(1:const.AN, :);
    [d, Z, transform]= procrustes(an_nodes, f_an);
    newC= zeros(const.N, 2);
    newC(1:const.AN, :)= transform.c;
    for j= const.AN+1:const.N
        count= 0;
        for i= 1:const.AN
            if(adj(i,j) > 0)
                newC(j,:)= newC(j,:) + newC(i,:);
                count= count + 1;
            end;
        end;
        if(count~= 0)
            newC(j,:)= newC(j,:)/count;
        end;
    end;
    f= (transform.b * f * transform.T) + newC;
    f(1:const.AN, :)= an_nodes;
    rms_index= rms_index+1;
    rms_error(rms_index)= myrms(nodes, f, const.N);
    pre_error= now_error;
    now_error= rms_error(rms_index);
end;
%% Final result
figure('Name', 'Procrustes');
plot(nodes(1:const.AN, 1), nodes(1:const.AN, 2), '^', ...
    nodes(const.AN+1:const.N, 1), nodes(const.AN+1:const.N, 2), 'o', ...
    f(const.AN+1:const.N, 1), f(const.AN+1:const.N, 2), '+');
hold on;
lineX= [nodes(const.AN+1:const.N, 1) f(const.AN+1:const.N, 1)];
lineY= [nodes(const.AN+1:const.N, 2) f(const.AN+1:const.N, 2)];
for i= 1:const.N-const.AN
    line(lineX(i, :), lineY(i, :), 'LineStyle', '-', 'Color', 'k');
    hold on;
end;
rectangle('Position', [const.MIN_D, const.MIN_D, ...
    const.MAX_D-const.MIN_D, const.MAX_D-const.MIN_D], 'EdgeColor', 'g');
legend('Anchor', 'Original', 'Procrustes');
str = ['RMS Error= ', num2str(rms_error(rms_index))];
title(str);