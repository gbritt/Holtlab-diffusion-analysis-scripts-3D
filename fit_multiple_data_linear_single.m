function fit_multiple_data_linear_single(thresh_length)

[filename,path] = uigetfile('.mat','multiselect','off');
cd(path)

f_lin = fittype('a*x');

% h = waitbar(0,'Fitting and extracting data...');

% length(filename)

% for m = 1:length(filename)
    
    result = struct();

    result = importdata(filename);

    num_tracks = size(result,1);
    if num_tracks == 1
       num_tracks = size(result,2);
    end
    
    isolate_idx = [];
    
    for j = 1:num_tracks
        if length(result(j).tracking.x) >= thresh_length
            isolate_idx = [isolate_idx;j];
        end
    end
    
    size_rest_data = length(isolate_idx);
    
    D = [];

    for j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    MSD = result(isolate_idx(j)).tracking.MSD;
    
    [yy_lin] = fit(time(1:10),MSD(1:10),f_lin,...
                'display','off');
            
    D = [D;yy_lin.a];
    
    saving_name = strcat('extracted_lin_',filename);
    extract = struct('D',D);
    save(saving_name,'extract')

%     waitbar(m/length(filename))
    
    end



disp('Data saved')
% close(h)

end