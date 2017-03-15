function results_tracking_V3(dt,conv)

% close all
% clear all


[filename,path] =uigetfile('multiselect','on','.txt','Select the file to convert');
cd(path)

% for m = 1:length(filename)
% filename = '/Users/Morgan/Desktop/GFA1 muflu device/P=0 - prepacked/ResultsP0.txt';
delimiter = '\t';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

conv = 0.130;
dt = 0.1;

Trajectory = dataArray{:, 2};
Frame = dataArray{:, 3};
x = dataArray{:, 4};
y = dataArray{:, 5};
z = dataArray{:, 6};
m0 = dataArray{:, 7};
m1 = dataArray{:, 8};
m2 = dataArray{:, 9};
m3 = dataArray{:, 10};
m4 = dataArray{:, 11};
NPscore = dataArray{:, 12};

idx = find(diff(Trajectory) > 0);

S = [];
for i = 1:length(idx)
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

MSD = zeros(1,length(idx));
time = zeros(1,length(idx));
x_res = zeros(1,length(idx));
y_res = zeros(1,length(idx));
frame_res = zeros(1,length(idx));

for i = 1:length(idx)
    
    if i == 1
        x_res(1:idx(i),i) = x(1:idx(i));
        y_res(1:idx(i),i) = y(1:idx(i));
        frame_res(1:idx(i),i) = Frame(1:idx(i));
        time(1:idx(i),i) = 0:dt:dt*(idx(i)-1);
        MSD(1:idx(i),i) = calculate_MSD(x(1:idx(i)),y(1:idx(i)),0,dt,conv);
      
    elseif i > 1
        x_res(1:(idx(i)-idx(i-1)),i) = x((idx(i-1)+1):idx(i));
        y_res(1:(idx(i)-idx(i-1)),i) = y((idx(i-1)+1):idx(i));
        frame_res(1:(idx(i)-idx(i-1)),i) = Frame((idx(i-1)+1):idx(i));
        time(1:(idx(i)-idx(i-1)),i) = 0:dt:dt*(idx(i)-idx(i-1)-1);
        MSD(1:(idx(i)-idx(i-1)),i) = calculate_MSD(x((idx(i-1)+1):idx(i)),y((idx(i-1)+1):idx(i)),0,dt,conv);
        
    end 
    
    result(i).tracking = struct('time',time(find(MSD(:,i))>0,i),...
        'x',x_res(find(MSD(:,i))>0,i),...
        'y',y_res(find(MSD(:,i))>0,i),...
        'MSD',MSD(find(MSD(:,i))>0,i),...
        'frame',frame_res(find(MSD(:,i))>0,i));
    
    
end

filename_int = filename;
name_file = filename_int(1:strfind(filename_int,'.txt')-1);

saving_name = strcat('tracked_',name_file,'.mat');
save(saving_name,'result')

% end

% uisave('result')
% save P=8_bis.mat result 

end

