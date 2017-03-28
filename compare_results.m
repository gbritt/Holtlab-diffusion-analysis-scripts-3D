function compare_results
% Function that compare the results obtained for the hop diffusion in
% violin plots


[filename,path] = uigetfile('multiselect','on','.mat');
cd(path)

data_Rc = cell(1,length(filename));
data_Din = cell(1,length(filename));
data_Dout = cell(1,length(filename));
data_Tin = cell(1,length(filename));
data_alpha = cell(1,length(filename));
name = cell(1,length(filename));

for i = 1:length(filename)
    
    res = importdata(filename{i});
    data_Rc{i} = res.Rc;
    data_Din{i} = res.Din;
    data_Dout{i} = res.Dout;
    data_Tin{i} = res.Tres;
    data_alpha{i} = res.alpha;
    name{i} = filename{i};
     
end
    
figure('position',[50 300 2500 600])
subplot(1,5,1)
title('Diffusion in corral')
distributionPlot(data_Din,'xName',name,'yLabel','Diffusion in corral (um^2/s)')
subplot(1,5,2)
title('Diffusion out of corral')
distributionPlot(data_Dout,'xName',name,'yLabel','Diffusion out of corral (um^2/s)')
subplot(1,5,3)
title('Size of the corral')
distributionPlot(data_Rc,'xName',name,'yLabel','Corral radius (um)')
subplot(1,5,4)
title('Residency time in the corral')
distributionPlot(data_Tin,'xName',name,'yLabel','Residency time (s)')
subplot(1,5,5)
title('Abnormality of diffusion out of the corral')
distributionPlot(data_alpha,'xName',name,'yLabel','Alpha')









end