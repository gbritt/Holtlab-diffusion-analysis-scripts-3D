function compare_results_linear_power
% Function that compare the results obtained for the hop diffusion in
% violin plots


[filename,path] = uigetfile('multiselect','on','.mat');
cd(path)

data_D = cell(1,length(filename));
% data_Dalpha = cell(1,length(filename));
% data_alpha = cell(1,length(filename));

name = cell(1,length(filename));
figure('position',[50 300 2000 600])
for i = 1:length(filename)
    
    res = importdata(filename{i});
    data_D{i} = res.D;
%     data_Dalpha{i} = res.Dalpha;
%     data_alpha{i} = res.alpha;
    name{i} = filename{i};
    subplot(1,2,2)
    title('D_\alpha VS \alpha')
    hold all
    plot(res.Dalpha,res.alpha,'o') 
    axis([0 10 0 2])

end
    

subplot(1,2,1)
title('Diffusion calculated from lienar fitting')
distributionPlot(data_D,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})')











end