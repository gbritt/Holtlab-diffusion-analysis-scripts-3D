function fit_lin_corr_data_v2(dt,conv,D0,min_track_length_lin,min_track_length_expo,sliding_size,sliding_step)

[filename,path] = uigetfile('multiselect','on','.mat');
cd(path)

flin = fittype('a*x');
fexp = fittype('a^2*(1-exp(-b*x))');

h = waitbar(0,'Fitting and extracting data...');

for k = 1:length(filename)

result = struct();
result = importdata(filename{k});

num_tracks = size(result,1);
if num_tracks == 1
   num_tracks = size(result,2);
end

%% First step, linear fit of all the trajectories that have a min length > min_track_length_lin
 
isolate_idx = [];
    
if min_track_length_lin < 11
    min_track_length_lin = 11;
end

for j = 1:num_tracks
    if length(result(j).tracking.x) >= min_track_length_lin
        isolate_idx = [isolate_idx;j];
    end
end
    
size_rest_data = length(isolate_idx);
    
Dlin = [];

parfor j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    MSD = result(isolate_idx(j)).tracking.MSD;
    
    [yy_lin] = fit(time(1:10),MSD(1:10),flin,...
                'display','off','Startpoint',[D0]);
            
    Dlin = [Dlin;yy_lin.a/4];

    
end

lin_fit = cell(1,1);
lin_fit = {Dlin};


%% Second step, corraled analysis of all the trajectories that have a min length > min_track_length_expo

isolate_idx = [];
    
if min_track_length_expo < min_track_length_lin
    min_track_length_expo = min_trakc_length_lin;
end

for j = 1:num_tracks
    if length(result(j).tracking.x) >= min_track_length_expo
        isolate_idx = [isolate_idx;j];
    end
end

size_rest_data = length(isolate_idx);

corr_fit = cell(size_rest_data,3);

    
parfor j = 1:size_rest_data
    
    MSD_length = length(result(isolate_idx(j)).tracking.x);
    chunks_number = floor((MSD_length-sliding_size)/sliding_step);
    Dexpo = zeros(1,chunks_number);
    R_c = zeros(1,chunks_number);
    
    
    for i = 1:chunks_number
        t_res(i) = result(isolate_idx(j)).tracking.time(1+(i-1)*sliding_step);
    end
    
    
    MSD_sliding = zeros(sliding_size,chunks_number);
    t_res = zeros(1,chunks_number);
       
    for i = 1:chunks_number
%         if exist_z == 0
            MSD_sliding(:,i) = calculate_MSD(result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                0,dt,conv);
            MSD_sliding(:,i) = MSD_sliding(:,i) - MSD_sliding(1,i);
%         elseif exist_z == 1
%             MSD_sliding(:,i) = calculate_MSD(result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
%                 result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
%                 result(isolate_idx(j)).tracking.z(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
%                 dt,conv);
%             MSD_sliding(:,i) = MSD_sliding(:,i) - MSD_sliding(1,i);
%         end
    end
    
    
   

for i = 1:chunks_number
    
    t = result(isolate_idx(j)).tracking.time(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
    t_fit = t-t(1);
    [yy_exp,gof_exp] = fit(t_fit(1:end-1),...
        MSD_sliding(1:end-1,i),fexp,'Lower',[0 0],'Robust','BiSquare','StartPoint',[2*sqrt(D0) 1],'display','off');
    [yy_lin,gof_lin] = fit(t_fit(1:end-1),...
        MSD_sliding(1:end-1,i),flin,'Robust','BiSquare','StartPoint',[D0],'display','off');

    lin_or_expo = gof_lin.rsquare/gof_exp.rsquare;
    
    if lin_or_expo < 0.7
                
        Dexpo(i) = yy_exp.b*yy_exp.a^2/4;
        R_c(i) = sqrt(5/2)*yy_exp.a;

    
    elseif lin_or_expo > 0.7    
        
        Dexpo(i) = yy_lin.a/4;
        R_c(i) = 10;
        
        
    end
    
    
        
   
end

    corr_fit(j,:) = {{Dexpo},{R_c},{isolate_idx(j)}};
    
end


row_headings = {'D_corr','R_c','track_index'};
extract_corr = cell2struct(corr_fit,row_headings,2);
extract_lin = cell2struct(lin_fit,'D_lin',1);

extract = struct('lin',extract_lin,'corr',extract_corr);

saving_name = strcat('extracted_',filename{k});

save(saving_name,'extract')



waitbar(k/length(filename))

end


disp('Extracted data saved')
close(h)











end