function [n,der,center,radius,idx] = number_cluster(time,x,y,num_max,thresh,plot_bool)

%% Find probable number of clusters

dev = zeros(1,num_max);

time = 1000*time/max(time);

for i = 1:num_max
    [idx,c] = kmedoids([time;x;y]',i,'algorithm','clara','replicates',50,'Distance','seuclidean');
%     idx = clusterdata([time;x;y]',i);
%     c = 
    dev_int = zeros(1,i);
    for j = 1:i
        dev_int(j) = std(sqrt((x(idx==j)-c(j,1)).^2 + (y(idx==j)-c(j,2)).^2));
    end
    dev(i) = mean(dev_int);
end

dev = 1-dev/dev(1);
der = diff(diff(dev));

% n = min(find(der < thresh));
n = find(der == min(der))+1;
[idx,c] = kmedoids([time;x;y]',n,'algorithm','clara','replicates',100,'Distance','seuclidean');

%% Extract informations about clusters

radius = zeros(1,n);
center = zeros(2,n);
% qual = zeros(1,n);
for i = 1:n
    [center_int,radius_int] = minboundcircle(x(idx == i),y(idx == i),0);
    center(:,i) = center_int;
    radius(i) = radius_int;
end

%% Check if clusters are not overlapping

non_cluster = zeros(n,2);

for i = 1:n
    for j = i+1:n
        bool = is_cluster(center(:,i),radius(i),center(:,j),radius(j),thresh);
        if bool == 0
            non_cluster(i,:) = [i j];
        end
    end
end

idx_non_cluster = find(non_cluster(:,2));
non_cluster_pair = non_cluster(idx_non_cluster,:);
n_unique = numel(find(unique(non_cluster_pair)));
n_pairs = size(non_cluster_pair,1);

if numel(find(unique(non_cluster_pair))) > 0

    n_real = n  - numel(find(unique(non_cluster_pair))) +1;
    radius_real = zeros(1,n_real);
    center_real = zeros(2,n_real);
    idx_real = idx;
    
    for i = 1:n_pairs
%         non_cluster_pair(i,1)
        T1 = abs(min(time(find(idx_real==non_cluster_pair(i,1)))) - max(time(find(idx_real==non_cluster_pair(i,2)))));
        T2 = abs(max(time(find(idx_real==non_cluster_pair(i,1)))) - min(time(find(idx_real==non_cluster_pair(i,2)))));
        delta_time = min([T1;T2]);
        if delta_time < 2*mean(diff(time))
            idx_real( find(idx_real==non_cluster_pair(i,1)) ) = non_cluster_pair(i,2);
        end
    end

    idx_unique = unique(idx_real);

    for i = 1:n_real
        [center_int,radius_int] = minboundcircle(x(idx_real == idx_unique(i)),y(idx_real == idx_unique(i)),0);
        center_real(:,i) = center_int;
        radius_real(i) = radius_int;
    end

    
end

%% Threshold on length of track

%% Plot data
if numel(find(unique(non_cluster_pair))) > 0
    n = n_real;
    idx = idx_real;
    center = center_real;
    radius = radius_real;
end

if plot_bool == 1

idx_unique = unique(idx);

figure;
hold all
t = linspace(0,2*pi);
for i = 1:n
    xr = x(idx==idx_unique(i));
    yr = y(idx==idx_unique(i));
%     for j = 1:length(xr)-1
%         l = sqrt((xr(j+1)-xr(j)).^2 + (yr(j+1)-yr(j)).^2);
%     end
%     qual(i) = l/radius(i);
    plot(xr,yr);
    plot(center(1,i)+radius(i)*cos(t),center(2,i)+radius(i)*sin(t),'r');
end
axis equal

end
   


end