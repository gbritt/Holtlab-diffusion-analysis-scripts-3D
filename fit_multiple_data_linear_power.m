function fit_multiple_data_linear_power(thresh_length)

[filename,path] = uigetfile('.mat','multiselect','on');
cd(path)

f_lin = fittype('a*x');

h = waitbar(0,'Fitting and extracting data...');

length(filename)

for m = 1:length(filename)
    
    result = struct();

    result = importdata(filename{m});

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
    Dalpha = [];
    alpha = [];

    for j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    MSD = result(isolate_idx(j)).tracking.MSD;
    
    MSD = MSD - MSD(1);
    time = time - time(1) + 0.00001;;
    
    yy_lin = fit(time(1:10),MSD(1:10),f_lin,...
                'display','off');
    yy_pow = fit(time(1:end-1),MSD(1:end-1),'power1',...
                'Startpoint',[yy_lin.a 1]);
            
    Dalpha = [Dalpha;yy_pow.a];
    alpha = [alpha;yy_pow.b];        
    D = [D;yy_lin.a/4];
    
    saving_name = strcat('extracted_lin_',filename{m});
    extract = struct('D',D,'Dalpha',Dalpha,'alpha',alpha);
    save(saving_name,'extract')

    waitbar(m/length(filename))
    
    end

end

disp('Data saved')
close(h)

end