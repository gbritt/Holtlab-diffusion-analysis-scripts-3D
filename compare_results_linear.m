function compare_results_linear
% Function that compare the results obtained for the hop diffusion in
% violin plots


[filename,path] = uigetfile('multiselect','on','.mat');
cd(path)

data_D = cell(1,length(filename));
name = cell(1,length(filename));

for i = 1:length(filename)
    
    res = importdata(filename{i});
    data_D{i} = res.D;
    name{i} = filename{i};
     
end
    
figure('position',[50 300 2000 600])
% subplot(1,5,1)
title('Diffusion calculated from linear fitting')
distributionPlot(data_D,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})')










end