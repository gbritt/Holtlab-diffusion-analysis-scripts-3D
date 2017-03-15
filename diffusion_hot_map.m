function diffusion_hot_map

close all

[filename,path] =uigetfile('multiselect','off','.tif','Select file');
cd(path)

I = imread(filename);

[filename,path] =uigetfile('multiselect','off','.mat','Select file');
cd(path)

result = importdata(filename);

isolate_idx = [];

min_track_length_lin = 11;

num_tracks = size(result,1);
if num_tracks == 1
   num_tracks = size(result,2);
end

for j = 1:num_tracks
    if length(result(j).tracking.x) >= min_track_length_lin
        isolate_idx = [isolate_idx;j];
    end
end
    
size_rest_data = length(isolate_idx);
    
Dlin = [];
track = [];
x_m = [];
y_m = [];
Rg = [];

for j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    MSD = result(isolate_idx(j)).tracking.MSD;
    
    x = result(isolate_idx(j)).tracking.x;
    y = result(isolate_idx(j)).tracking.y;
    
    x_m = [x_m;mean(x)];
    y_m = [y_m;mean(y)];
    
%     Rg = [Rg; sqrt(numel(x)^(-1)*(sum((x-x_m).^2 + (y-y_m).^2)))];
    
    [yy_lin] = fit(time(1:10),MSD(1:10),'a*x',...
                'display','off','Startpoint',0.2);
    
    track = [track;result(isolate_idx(j)).tracking.frame];        
    Dlin = [Dlin;yy_lin.a/4];

    
end

J = zeros(size(I));

for i = 1:length(x_m)
    if floor(x_m(i)) > 1 && floor(y_m(i))>1 && Dlin(i) < 0.8 && Dlin(i) > 0.2
        J(floor(y_m(i)),floor(x_m(i))) = 1000*Dlin(i);
        J(floor(y_m(i))-1,floor(x_m(i))) = 1000*Dlin(i);
        J(floor(y_m(i))+1,floor(x_m(i))) = 1000*Dlin(i);
        J(floor(y_m(i)),floor(x_m(i))+1) = 1000*Dlin(i);
        J(floor(y_m(i))-1,floor(x_m(i))+1) = 1000*Dlin(i);
        J(floor(y_m(i))+1,floor(x_m(i))+1) = 1000*Dlin(i);
        J(floor(y_m(i)),floor(x_m(i))-1) = 1000*Dlin(i);
        J(floor(y_m(i))-1,floor(x_m(i))-1) = 1000*Dlin(i);
        J(floor(y_m(i))+1,floor(x_m(i))-1) = 1000*Dlin(i);
    end
end


figure
imagesc(I);
axis equal
colormap gray
axis([150 500 0 300])

figure
imagesc(J)
axis equal
% colormap jet
axis([150 500 0 300])

cmap = jet(800);
% Make values 0-5 black:
cmap(1,:) = zeros(1,3);
colormap(cmap);
colorbar




end