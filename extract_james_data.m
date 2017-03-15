function extract_james_data

% close all

[filename,path] = uigetfile('multiselect','on','.mat');
cd(path)

track_length = [];
power_fit = [];
K_fit = [];
D_lin = [];


for k = 1:length(filename)

result = struct();
result = importdata(filename{k});

isolate_idx = [];

num_tracks = size(result,1);
if num_tracks == 1
   num_tracks = size(result,2);
end

for j = 1:num_tracks
    track_length = [track_length;length(result(j).tracking.x)];
    if length(result(j).tracking.x) >= 50
        isolate_idx = [isolate_idx;j];
    end
end
    
size_rest_data = length(isolate_idx);

for j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    MSD = result(isolate_idx(j)).tracking.MSD;
    MSD = MSD-MSD(1);
    
    
    
    [yy,gof] = fit(time(2:end-1),MSD(2:end-1),'power1', 'display','off','lower',[0 0],'upper',[1 3]);
    
    
    
    if gof.rsquare > 0.70
        
%         f = figure;
%         plot(time,MSD,'b')
%         hold all
%         plot(yy)
%         pause
%         close(f)
       
        power_fit = [power_fit;yy.b];
        K_fit = [K_fit;yy.a];
        
        [yy,gof] = fit(time(1:3),MSD(1:3),'a*x', 'display','off');
        D_lin = [D_lin;yy.a];
        
    end

    
end

end

figure
plot(K_fit,power_fit,'bo')

figure
hist(D_lin,40)

figure
hist(track_length,100)



















end