function convert_u_track_data(dt,conv)

[filename,path] = uigetfile('.mat','multiselect','off');
cd(path)

data = struct();
data = importdata(filename);

track_pars = [11 0 0.1];

[traceList,info] = tracksFinal_2_traceList(data,track_pars);

num_traces = info.Ntraces;

    S = [];
    for i = 1:num_traces
        if i < 10
            S = [S;strcat('part0000',num2str(i))];
        elseif i >= 10 && i < 100
            S = [S;strcat('part000',num2str(i))];
        elseif i >= 100 && i < 999
            S = [S;strcat('part00',num2str(i))];
        elseif i >=1000 && i < 9999
            S = [S;strcat('part0',num2str(i))];
        elseif i >= 10000 && i < 99999
            S = [S;strcat('part',num2str(i))];
        end
    end

    field = cellstr(S);

    result = cell2struct(field','tracking',1);
for i = 1:num_traces
    
    x = traceList{i}.x;
    y = traceList{i}.y;
    frame = traceList{i}.frame;
    
    [MSD,time] = calculate_MSD_V2(x,y,0,dt,conv);
    

    result(i).tracking = struct('time',time',...
        'x',x,...
        'y',y,...
        'MSD',MSD',...
        'frame',frame);
    
end


saving_name = strcat('tracked_',filename,'.mat');
save(saving_name,'result')






end