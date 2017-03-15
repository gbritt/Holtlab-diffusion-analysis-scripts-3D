function fit_lin_gyration_V2(dt,conv,D0,min_track_length_lin,min_track_length_expo,sliding_size,sliding_step,Tergo)

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
exist_z = result(1).tracking.z > 0;

%% Linear fit of all the trajectories that have a min length > min_track_length_lin
 
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
    
Dlin = zeros(1,size_rest_data);
Dlin_err = zeros(1,size_rest_data);
Dlin_gof_rsquare = zeros(1,size_rest_data);
Dlin_gof_rmse = zeros(1,size_rest_data);

parfor j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    MSD = result(isolate_idx(j)).tracking.MSD;
    
    [yy_lin,gof] = fit(time(1:10),MSD(1:10),flin,...
                'display','off','Startpoint',[D0]);
  
    Dlin(j) = yy_lin.a/4;
    err = confint(yy_lin);
    Dlin_err(j) = err(2)-err(1);
    Dlin_gof_rsquare(j) = gof.rsquare;
    Dlin_gof_rmse(j) = gof.rmse;
            
    

    
end

lin_fit = {{Dlin},{Dlin_err},{Dlin_gof_rsquare},{Dlin_gof_rmse}};


%% Corraled analysis of all the trajectories that have a min length > min_track_length_expo

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
    exist_z = result(j).tracking.z > 0
    for i = 1:chunks_number
        if exist_z == 0
            MSD_sliding(:,i) = calculate_MSD(result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                0,dt,conv);
            MSD_sliding(:,i) = MSD_sliding(:,i) - MSD_sliding(1,i);
        elseif exist_z == 1
            MSD_sliding(:,i) = calculate_MSD(result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.z(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                dt,conv);
            MSD_sliding(:,i) = MSD_sliding(:,i) - MSD_sliding(1,i);
        end
    end
    
    
   

for i = 1:chunks_number
    
    t = result(isolate_idx(j)).tracking.time(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
    t_fit = t-t(1);
    
    x = conv*result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
    y = conv*result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
    z = conv*result(isolate_idx(j)).tracking.z(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
    
    mean_x = mean(x);
    mean_y = mean(y);
    mean_z = mean(z);
  
    % end
    mean_x = mean(x);
    mean_y = mean(y);
    
   % if exist_z == 1 %need to troubleshoot parfor loop here
        Rg = sqrt(numel(x)^(-1)*(sum((x-mean_x).^2 + (y-mean_y).^2 + (z-mean_z).^2)))
        
%     else
%         Rg = sqrt(numel(x)^(-1)*(sum((x-mean_x).^2 + (y-mean_y).^2 )))
%     end
     yy_lin = fit(t_fit(1:11),...
        MSD_sliding(1:11,i),flin,'StartPoint',[D0],'display','off');

        R_c(i) = Rg;
        
        Dexpo(i) = yy_lin.a/4;

end

    corr_fit(j,:) = {{Dexpo},{R_c},{isolate_idx(j)}};
    
end

%% Calculating the ensemble average MSD and the mean-time average MSD

isolate_idx = [];
size_track = [];

    
if min_track_length_lin < 11
    min_track_length_lin = 11;
end

for j = 1:num_tracks

   if length(result(j).tracking.x) >= min_track_length_lin
        isolate_idx = [isolate_idx;j];
        size_track = [size_track;length(result(j).tracking.x)]; 
    end
end
    
size_rest_data = length(isolate_idx);

max_track = max(size_track);
MSD_e = zeros(1,max_track);
MSD_t = zeros(1,max_track);
x_MSD_ens = zeros(size_rest_data,max_track);
y_MSD_ens = zeros(size_rest_data,max_track);
z_MSD_ens = zeros(size_rest_data,max_track);
MSD_MSD_t = zeros(size_rest_data,max_track);

for i = 1:size_rest_data
    
    x_MSD_ens(i,1:size_track(i)) = result(isolate_idx(i)).tracking.x(1:end);
    y_MSD_ens(i,1:size_track(i)) = result(isolate_idx(i)).tracking.y(1:size_track(i));
    MSD_MSD_t(i,1:size_track(i)-1) = result(isolate_idx(i)).tracking.MSD(1:size_track(i)-1);    
    if exist_z == 1
        z_MSD_ens(i,1:size_track(i)) = result(isolate_idx(i)).tracking.z(1:size_track(i));
    end
end
for i = 1:max_track
    idx = find(x_MSD_ens(:,i) > 0);
    MSD_e(i) = conv^2*mean( ( (x_MSD_ens(idx,i) - x_MSD_ens(idx,1)).^2 + (y_MSD_ens(idx,i) - y_MSD_ens(idx,1)).^2 ) );
    if exist_z == 1
        MSD_e(i) = conv^2*mean( ( (x_MSD_ens(idx,i) - x_MSD_ens(idx,1)).^2 + (y_MSD_ens(idx,i) - y_MSD_ens(idx,1)).^2  ...
            + (z_MSD_ens(idx,i) - z_MSD_ens(idx,1)).^2 ) );
    end
    idx = find(MSD_MSD_t(:,i) > 0);
    MSD_t(i) = mean(MSD_MSD_t(idx,i));
end

MSD = {{MSD_e},{MSD_t}};

%% Calculating the ergodicity

isolate_idx = [];

for j = 1:num_tracks
   if length(result(j).tracking.x) >= Tergo %not enough input arguments?
        isolate_idx = [isolate_idx;j];
   end
end

size_rest_data = length(isolate_idx);

eps = zeros(size_rest_data,Tergo);

for i = 1:size_rest_data
    eps(i,1:Tergo) = result(isolate_idx(i)).tracking.MSD(1:Tergo)./MSD_e(1:Tergo)';
end

ergo = {eps};

%% Calculating the displacement distribution and orientation distribution

isolate_idx = [];

for j = 1:num_tracks
   if length(result(j).tracking.x) >=  min_track_length_lin
        isolate_idx = [isolate_idx;j];
   end
end

size_rest_data = length(isolate_idx);
disp_res = cell(1,Tergo);
disp_ori = cell(1,Tergo);

displ = [];
ori = [];

for i = 1:Tergo
    for j = 1:size_rest_data  
        displ = [displ;conv*displacement(result(isolate_idx(j)).tracking.x,result(isolate_idx(j)).tracking.y,result(isolate_idx(j)).tracking.z,i)']; % use numerical and z?
        ori = [ori;disp_corr(result(isolate_idx(j)).tracking.x,result(isolate_idx(j)).tracking.y,result(isolate_idx(j)).tracking.z,i)'];
    end
    disp_res{i} = {displ};
    disp_ori{i} = {ori};
end



%% Calculating the next step correlation

next_step_res = [];
step_next_size = [];

for j = 1:size_rest_data 
    [next_step_calc,step_next_size_calc] = next_step(result(isolate_idx(j)).tracking.x,result(isolate_idx(j)).tracking.y,result(isolate_idx(j)).tracking.z);
    next_step_res = [next_step_res;conv*next_step_calc'];
    step_next_size = [step_next_size;conv*step_next_size_calc'];
end

next = {{next_step_res},{step_next_size}};


%% Last step, saving the data


row_headings = {'D_corr','R_c','track_index'};
extract_corr = cell2struct(corr_fit,row_headings,2);
row_headings = {'D_lin','D_lin_err','gof_rsquare','gof_rmse'};
extract_lin = cell2struct(lin_fit,row_headings,2);
row_headings = {'MSD_ens','MSD_time'};
extract_MSD = cell2struct(MSD,row_headings,2);
extract_ergo = cell2struct(ergo,'ergo',1);
extract_step = cell2struct(disp_res,'step',1);
extract_ori = cell2struct(disp_ori,'ori',1);
row_headings = {'corr','step_size'};
extract_next = cell2struct(next,row_headings,2);

extract = struct('lin',extract_lin,...
    'corr',extract_corr,...
    'MSD',extract_MSD,...
    'ergo',extract_ergo,...
    'step',extract_step,...
    'ori',extract_ori,...
    'next',extract_next);

saving_name = strcat('analyzed_',filename{k});

save(saving_name,'extract')



waitbar(k/length(filename))

end


disp('Extracted data saved')
close(h)% this broke?











end
