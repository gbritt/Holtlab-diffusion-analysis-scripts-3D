function observe_tracking(exist_z)
% This function is designed to observe the tracking of particles in cells.
% Two finles need to be loaded for this. A structured matlab file,
% containing the time, x and y coordinates, MSD and frame number of a
% particular track, as well as the piture file. Exist_z takes a value 0 if
% these are 2D data, and a value 1 if there is a third dimension. 

close all
% clear all
clc




fig = figure('position',[20 20 1330 650]);


global I
global I_zoom
global result
global num_tracks
global hresults_button
global isolate_idx
global size_tif
global mean_x
global mean_y
global size_zoom
global conv
global dt
global filename
% global exist_z

isolate_idx = [];

h_load_images = uicontrol('Style','Pushbutton','String','Load .tif','Position',[30 100 70 25],...
    'Callback',{@load_image_callback},'visible','off');
h_dt_text = uicontrol('Style','Text','String','dt =','Position',[30 70 35 25],'visible','off');
h_dt_value = uicontrol('Style','Edit','String','0.1','Position',[65 70 35 25],'visible','off');
h_conv_text = uicontrol('Style','Text','String','conv =','Position',[30 40 35 25],'visible','off');
h_conv_value = uicontrol('Style','Edit','String','0.13','Position',[65 40 35 25],'visible','off');
h_3D_value = uicontrol('Style','popupmenu','String','2D|3D','Position',[30,10,75,25],'Visible','off',...
    'callback',{@threeD_callback});
h_slider_image = uicontrol('Style','slider','Position',[745 10 470 30],'visible','off',...
    'min',1,'callback',{@slider_image_callback},'value',1);
uicontrol('Style','Pushbutton','String','Load .mat','Position',[30 130 70 25],'Callback',{@load_mat_callback});
h_thresh_length = uicontrol('Style','Edit','String','40','Position',[540 300 70 25],'visible','off');
h_thresh_length_text = uicontrol('Style','text','String','Thresh. length','Position',[540 320 70 25],'visible','off');
h_replot_button = uicontrol('Style','Pushbutton','String','Replot','Position',[540 270 70 25],...
    'visible','off','callback',{@replot_callback});
h_isolate_button = uicontrol('Style','Pushbutton','String','Isolate','Position',[540 240 70 25],...
    'visible','off','callback',{@isolate_callback});
h_clear_isolate_button = uicontrol('Style','Pushbutton','String','Clear Isolate','Position',[540 210 70 25],...
    'visible','off','callback',{@clear_isolate_callback});
h_zoom_trace = uicontrol('Style','Text','String','Zoom Trace','Position',[540 180 70 25],'Visible','off');
h_zoom_trace_box = uicontrol('Style','Checkbox','Position',[540 160 25 25],'Visible','off','callback',{@zoom_trace_callback});
h_zoom_trace_size = uicontrol('Style','Edit','String','20','Position',[580 160 25 25],'Visible','off');
h_extract_isolate_curve = uicontrol('Style','Pushbutton','String','Extract Isolate','Position',[540 130 70 25],...
    'Visible','off','Callback',{@extract_isolate_curve});
h_sliding_MSD = uicontrol('Style','Pushbutton','String','Sliding MSD','Position',[540 100 70 25],...
    'Visible','off','Callback',{@sliding_MSD_callback});
h_sliding_step = uicontrol('Style','Edit','String','1','Position',[580 70 25 25],'Visible','off');
h_sliding_size = uicontrol('Style','Edit','String','15','Position',[540 70 25 25],'Visible','off');
h_save_figure = uicontrol('Style','Pushbutton','String','Save Data','Position',[540 40 70 25],...
    'Visible','off','Callback',{@save_figure_callback});
h_extract = uicontrol('Style','Pushbutton','String','Extract Data',...
    'Position',[540 10 70 25],'Callback',{@extract_data_sliding_callback},'Visible','off');

columnname =   {'Trajectory','Length','Visible ?', 'Isolate ?'};
columnformat = {'numeric','numeric',[],[]}; 
columneditable =  [false edit true true]; 
hresults_button = uitable('Position',...
            [180,10,350,330], ... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columneditable,...
            'visible','off');

function threeD_callback (source,eventdata)
        
exist_z = get(h_3D_value,'Value')-1;  
end

function save_figure_callback (source,eventdata)

dt = str2double(get(h_dt_value,'String'));
conv = str2double(get(h_conv_value,'String'));    
    
dat = get(hresults_button,'Data');
isolate_idx = find([dat{:,4}] == true);

size_rest_data = length(isolate_idx);

S = [];
for i = 1:size_rest_data
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

rest = cell2struct(field','result',1);



for i = 1:size_rest_data
    new_size = cell2mat(dat(isolate_idx(i),2));
    rest(i).result = struct('time',result(isolate_idx(i)).tracking.time(1:new_size),...
        'x',result(isolate_idx(i)).tracking.x(1:new_size),...
        'y',result(isolate_idx(i)).tracking.y(1:new_size),...
        'MSD',calculate_MSD(result(isolate_idx(i)).tracking.x(1:new_size),...
            result(isolate_idx(i)).tracking.y(1:new_size),0,dt,conv));
end


    
uisave('rest')
        
end

function extract_data_sliding_callback (source,eventdata)
        
dat = get(hresults_button,'Data');
isolate_idx = find([dat{:,4}] == true);    
    
size_rest_data = length(isolate_idx);

S = [];
for i = 1:size_rest_data
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

extract = cell2struct(field','result',1);


f = fittype('a^2*(1-exp(-b*x))');
sliding_step = str2double(get(h_sliding_step,'String'));
sliding_size = str2double(get(h_sliding_size,'String'));

h = waitbar(0,'Please wait...');

for j = 1:size_rest_data
    
    MSD_length = length(result(isolate_idx(j)).tracking.x);
    chunks_number = floor((MSD_length-sliding_size)/sliding_step);
    
    D = zeros(1,chunks_number);
    R_c = zeros(1,chunks_number);
    t_res = zeros(1,chunks_number);
    
    for i = 1:chunks_number
        t_res(i) = result(isolate_idx(j)).tracking.time(1+(i-1)*sliding_step);
    end
    
    MSD_sliding = zeros(sliding_size,chunks_number);
       
    for i = 1:chunks_number
        if exist_z == 0
            MSD_sliding(:,i) = calculate_MSD(result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                0,dt,conv);
        elseif exist_z == 1
                MSD_sliding(:,i) = calculate_MSD(result(isolate_idx(j)).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx(j)).tracking.z(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                dt,conv);
        end
        
    end
    

    for i = 1:chunks_number
          t = result(isolate_idx(j)).tracking.time(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
          t_fit = t-t(1);
          yy = fit(t_fit(1:end-1),...
              MSD_sliding(1:end-1,i),f,'Lower',[0 0],'Robust','BiSquare','StartPoint',[1,3]);

          D(i) = yy.b*yy.a/4;
          R_c(i) = sqrt(2*yy.a/5);
          
    end
        
    extract(j).result = struct('time',t_res,'D',D,'Rc',R_c);
end        
 
saving_name = strcat('extracted_',filename);

save(saving_name,'extract')
% uisave('extract')

disp('Data saved')
close(h)
        
end

function sliding_MSD_callback (source,eventdata)

dt = str2double(get(h_dt_value,'String'));
conv = str2double(get(h_conv_value,'String'));
    
if length(isolate_idx) == 1
    
    sliding_step = str2double(get(h_sliding_step,'String'));
    sliding_size = str2double(get(h_sliding_size,'String'));
    MSD_length = length(result(isolate_idx).tracking.x);
    chunks_number = floor((MSD_length-sliding_size)/sliding_step);
    
    D = zeros(1,chunks_number);
    R_c = zeros(1,chunks_number);
    t_res = zeros(1,chunks_number);
    
    for i = 1:chunks_number
        t_res(i) = result(isolate_idx).tracking.time(1+(i-1)*sliding_step);
    end
    
    
    MSD_sliding = zeros(sliding_size,chunks_number);

       
    for i = 1:chunks_number
        if exist_z == 0
            MSD_sliding(:,i) = calculate_MSD(result(isolate_idx).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                0,dt,conv);
        elseif exist_z == 1
            MSD_sliding(:,i) = calculate_MSD(result(isolate_idx).tracking.x(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx).tracking.y(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                result(isolate_idx).tracking.z(1+(i-1)*sliding_step:(i-1)*sliding_step + sliding_size),...
                dt,conv);
        end
    end
    
f = fittype('a^2*(1-exp(-b*x))');
   
figure('Position',[100 100 1200 500]);

for i = 1:chunks_number
    
    t = result(isolate_idx).tracking.time(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size);
    t_fit = t-t(1);
    yy = fit(t_fit(1:end-1),...
        MSD_sliding(1:end-1,i),f,'Lower',[0 0],'Robust','BiSquare','StartPoint',[1,3]);

    D(i) = yy.b*yy.a/4;
    R_c(i) = sqrt(2*yy.a/5);

    subplot(2,2,1)
    hold on
%     plot(t_fit+t(1),yy.a^2*(1-exp(-yy.b*t_fit)),'r')
    plot(result(isolate_idx).tracking.time(1+(i-1)*sliding_step:(i-1)*sliding_step+sliding_size),MSD_sliding(:,i),'b')
    
    subplot(2,2,2)
    plot(t_res(1:i),R_c(1:i),'r')
    
    subplot(2,2,3)
    plot(t_res(1:i),D(1:i),'r')
    
    subplot(2,2,4)
    plot(R_c(1:i),D(1:i),'ro')
end


    
    
end
    
end

function extract_isolate_curve (source,eventdata)

dt = str2double(get(h_dt_value,'String'));
conv = str2double(get(h_conv_value,'String'));
    
    if length(isolate_idx) == 1
       if exist_z == 0
        figure('Position',[100 100 1200 500])
        time = result(isolate_idx).tracking.time;
        MSD = result(isolate_idx).tracking.MSD;
        x = result(isolate_idx).tracking.x;
        y = result(isolate_idx).tracking.y;
        MSD_length = length(result(isolate_idx).tracking.x);
        SD = zeros(MSD_length,1);
        
        SD = conv^2*(result(isolate_idx).tracking.x - result(isolate_idx).tracking.x(1)).^2 + ...
        conv^2*(result(isolate_idx).tracking.y - result(isolate_idx).tracking.y(1)).^2;
        
        z = zeros(size(time));
        
        subplot(1,3,1)
        surface([time';time'],[MSD';MSD'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        subplot(1,3,2)
        surface([x';x'],[y';y'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        subplot(1,3,3)
        surface([time';time'],[SD';SD'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        colormap jet
        
       elseif exist_z == 1
        figure('Position',[100 100 1200 500])
        time = result(isolate_idx).tracking.time;
        MSD = result(isolate_idx).tracking.MSD;
        x = result(isolate_idx).tracking.x;
        y = result(isolate_idx).tracking.y;
        z_data = result(isolate_idx).tracking.z;
        MSD_length = length(result(isolate_idx).tracking.x);
        SD = zeros(MSD_length,1);
        
        SD = conv^2*(result(isolate_idx).tracking.x - result(isolate_idx).tracking.x(1)).^2 + ...
        conv^2*(result(isolate_idx).tracking.y - result(isolate_idx).tracking.y(1)).^2 + ...
        conv^2*(result(isolate_idx).tracking.z - result(isolate_idx).tracking.z(1)).^2;
        
        z = zeros(size(time));
        
        subplot(1,3,1)
        surface([time';time'],[MSD';MSD'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        subplot(3,3,2)
        surface([x';x'],[z_data';z_data'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        subplot(3,3,5)
        surface([y';y'],[z_data';z_data'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        subplot(3,3,8)
        surface([x';x'],[y';y'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        subplot(1,3,3)
        surface([time';time'],[SD';SD'],[z';z'],[time';time'],'facecol','no','edgecol','interp','linew',2)
        colormap jet
           
       end
    end

    
end

function zoom_trace_callback (source,eventdata)
   
    if get(h_zoom_trace_box,'Value') == 1
        if length(isolate_idx) == 1
            
    size_zoom = floor(str2double(get(h_zoom_trace_size,'String')));
    I_zoom = zeros(size_zoom+1,size_zoom+1,size_tif);
       
    mean_x = mean(result(isolate_idx).tracking.x);
    mean_y = mean(result(isolate_idx).tracking.y);
    
    rect_zoom = [floor(mean_x - size_zoom/2) floor(mean_y - size_zoom/2) size_zoom size_zoom];
    
    
           
            for i = 1:size_tif
            
               I_zoom(:,:,i) = imcrop(I(:,:,i),rect_zoom);
            end
            
        end
    end
end

function clear_isolate_callback (source,eventdata)
    
    if length(isolate_idx) >= 1
        
        dat = get(hresults_button,'Data'); 
        for i = 1:length(isolate_idx)
            dat(isolate_idx(i),4) = {false};
        end
        isolate_idx = [];
        
        set(hresults_button,'Data',dat)
        
    end
    
    thresh_length = str2double(get(h_thresh_length,'String'));

delete(subplot(2,2,1))
subplot(2,2,1)
hold all

dat=get(hresults_button,'Data');

for i = 1:num_tracks
    if length(result(i).tracking.x) >= thresh_length
        plot(result(i).tracking.time,result(i).tracking.MSD)
        dat(i,3) = {true};
    else
        dat(i,3) = {false};
    end
end

set(hresults_button,'Data',dat)
    
end

function isolate_callback (source,eventdata)

dt = str2double(get(h_dt_value,'String'));
conv = str2double(get(h_conv_value,'String'));    
    
dat = get(hresults_button,'Data');

isolate_idx = find([dat{:,4}] == true);
   
delete(subplot(2,2,1))
subplot(2,2,1)
hold all

% for i = 1:length(isolate_idx)
%     new_size = cell2mat(dat(isolate_idx(i),2));
%     if length(result(isolate_idx(i)).tracking.frame) > new_size
%         result(isolate_idx(i)).tracking.MSD(1:new_size) = calculate_MSD(result(isolate_idx(i)).tracking.x(1:new_size),...
%             result(isolate_idx(i)).tracking.y(1:new_size),0,dt,conv);
%     end
% end

for i = 1:length(isolate_idx)
%     plot(result(isolate_idx(i)).tracking.time(1:new_size),result(isolate_idx(i)).tracking.MSD(1:new_size))
    plot(result(isolate_idx(i)).tracking.time(),result(isolate_idx(i)).tracking.MSD())

end
    

end

function replot_callback (source,eventdata)
        
thresh_length = str2double(get(h_thresh_length,'String'));

delete(subplot(2,2,1))
subplot(2,2,1)
hold all

dat=get(hresults_button,'Data');


for i = 1:num_tracks
    if length(result(i).tracking.x) >= thresh_length
        plot(result(i).tracking.time,result(i).tracking.MSD)
        dat(i,3) = {true};
    else
        dat(i,3) = {false};
    end
end

set(hresults_button,'Data',dat)
        
        
end

function load_mat_callback (source,eventdata)

[filename,path] = uigetfile('multiselect','off','.mat');

cd(path)

result = struct();

result = importdata(filename);

num_tracks = size(result,1);
if num_tracks == 1
   num_tracks = size(result,2);
end

% exist_z = exist('result(1).tracking.z')

subplot(2,2,1)
hold all
for i = 1:num_tracks
    plot(result(i).tracking.time,result(i).tracking.MSD)
end


dat = cell(num_tracks,4);

for i = 1:num_tracks
    dat(i,1) =  {i};
    dat(i,2) =  {length(result(i).tracking.x)};
    dat(i,3) =  {true};
    dat(i,4) =  {false};

end

columnname =   {'Trajectory','Length','Visible ?', 'Isolate ?'};
columnformat = {'numeric','numeric',[],[]}; 
columneditable =  [false edit true true]; 
hresults_button = uitable('Position',...
            [180,10,350,330], 'Data', dat,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', columneditable);
        
set(h_load_images,'visible','on')   
set(h_thresh_length,'visible','on')
set(h_thresh_length_text,'visible','on')
set(h_replot_button,'visible','on')
set(h_isolate_button,'visible','on')
set(h_clear_isolate_button,'visible','on')
set(h_extract_isolate_curve,'visible','on')
set(h_sliding_step,'visible','on')
set(h_sliding_MSD,'visible','on')
set(h_sliding_size,'visible','on')
set(h_save_figure,'visible','on')
set(h_dt_text,'visible','on')
set(h_dt_value,'visible','on')
set(h_conv_text,'visible','on')
set(h_conv_value,'visible','on')
set(h_extract,'visible','on')
set(h_3D_value,'visible','on')
set(h_3D_value,'value',exist_z+1)

dt = str2double(get(h_dt_value,'String'));
conv = str2double(get(h_conv_value,'String'));


end

function load_image_callback (source,eventdata)
        
[filename,path] = uigetfile('multiselect','off','.tif');

cd(path)
size_tif = size(imfinfo(filename),1);

I_1 = imread(filename,'index',1);
m = size(I_1,1);
n = size(I_1,2);

I = zeros(m,n,size_tif);

for i = 1:size_tif
    I(:,:,i) = imread(filename,'index',i);
end

set(h_slider_image,'visible','on');
set(h_slider_image,'max',size_tif);



subplot(1,2,2)

imagesc(I(:,:,1))
axis equal
colormap gray

hold all

frame_num = floor(get(h_slider_image,'value'));
for i = 1:num_tracks
    n = find(result(i).tracking.frame == frame_num);
    if isempty(n) == 0
        plot(result(i).tracking.x(n),result(i).tracking.y(n),'o')
        text(result(i).tracking.x(n),result(i).tracking.y(n),num2str(i),'color','w')
    end
end

set(h_zoom_trace,'visible','on')
set(h_zoom_trace_box,'visible','on')
set(h_zoom_trace_size,'visible','on')  

        
end

function slider_image_callback (source,eventdata)

frame = floor(get(h_slider_image,'value'));

delete(subplot(1,2,2))

subplot(1,2,2)
if get(h_zoom_trace_box,'Value') == 0
    imagesc(I(:,:,frame))
elseif get(h_zoom_trace_box,'Value') == 1
    imagesc(I_zoom(:,:,frame))
end
axis equal
colormap gray

hold all

dat=get(hresults_button,'Data');
visible_idx = find([dat{:,3}] == true);

frame_num = floor(get(h_slider_image,'value'));
for i = 1:length(visible_idx)
    n = find(result(visible_idx(i)).tracking.frame == frame_num);
    if isempty(n) == 0 && isempty(isolate_idx) == 1
        plot(result(visible_idx(i)).tracking.x(1:n),result(visible_idx(i)).tracking.y(1:n),'-')
        text(result(visible_idx(i)).tracking.x(n),result(visible_idx(i)).tracking.y(n),num2str(visible_idx(i)),'color','w')
    end
end
    

if isempty(isolate_idx) == 0
    for i = 1:length(isolate_idx)
        n = find(result(isolate_idx(i)).tracking.frame == frame_num);
        if isempty(n) == 0
            if get(h_zoom_trace_box,'Value') == 1
                plot(result(isolate_idx(i)).tracking.x(1:n)-mean_x+size_zoom/2,result(isolate_idx(i)).tracking.y(1:n)-mean_y+size_zoom/2,'-')
            else
                plot(result(isolate_idx(i)).tracking.x(1:n),result(isolate_idx(i)).tracking.y(1:n),'-')
                text(result(isolate_idx(i)).tracking.x(n),result(isolate_idx(i)).tracking.y(n),num2str(isolate_idx(i)),'color','w')
            end
        end
    end
end
% delete(subplot(2,2,1))
% subplot(2,2,1)
% hold all
% 
% for i = 1:length(visible_idx)
%     plot(result(visible_idx(i)).tracking.time,result(visible_idx(i)).tracking.MSD)
%     plot([result(visible_idx(i)).tracking.time(frame_num) result(visible_idx(i)).tracking.time(frame_num)],[0 10],'k')
% end
    
end














end
