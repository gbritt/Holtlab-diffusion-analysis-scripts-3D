function compare_fit_lin_gyration(x_names,pooled,window_size,dt)

[filename,path] = uigetfile('multiselect','on','.mat');
cd(path)

figure('position',[50 300 2000 600]);

data_D_lin = cell(1,length(filename));
data_D_lin_rest = cell(1,length(filename));
data_D_corr = cell(1,length(filename));
data_D_corr_nopx = cell(1,length(filename));
data_D_corr_rest = cell(1,length(filename));
data_D_corr_out = cell(1,length(filename));
data_R_c = cell(1,length(filename));
data_R_c_nopx = cell(1,length(filename));
name = cell(1,length(filename));
% name = x_names;

for n  = 1:length(filename)
    
    res = importdata(filename{n});
    
    data_D_lin{n} = res.lin.D_lin;
    mean_D_lin = mean(res.lin.D_lin);
    
    if isempty(x_names) == 1
        name{n} = filename{n};
    else
        name = x_names;
    end


    size_1 = size(res.corr,1); % number of trajectories that have been analyzed with the corraled fitting
    for i = 1:size_1
       data_D_corr_nopx{n} = [data_D_corr_nopx{n} res.corr(i).D_corr{1}(res.corr(i).R_c{:} > 0.13)];
       data_D_corr{n} = [data_D_corr{n} res.corr(i).D_corr{:}];
%        data_D_corr_out{n} = [data_D_corr_out{n} res.corr(i).D_corr{1}(res.corr(i).R_c{:} > 1)];
       data_R_c{n} = [data_R_c{n} res.corr(i).R_c{:}];
       data_R_c_nopx{n} = [data_R_c_nopx{n} res.corr(i).R_c{1}(res.corr(i).R_c{:} > 0.13)];
    end
    
    mean_D_lin = mean(res.lin.D_lin);
    subplot(1,4,4)
    title('Diffusion VS R_g')
    hold all
%     length(data_R_c{n})
%     length(data_D_corr{n}+(n-1)/2)
    if pooled == 1
        dscatter(data_R_c{n},data_D_corr{n}+(n-1)*0.05,'MARKER','o','smoothing',2);
        plot(data_R_c{n},1/(window_size*dt*4*0.17)*data_R_c{n}.^2 + (n-1)*0.05,'r.')
    elseif pooled == 0
        dscatter(data_R_c{n}',data_D_corr{n}'+(n-1)*0.05,'MARKER','o','smoothing',2);
        plot(data_R_c{n},1/(window_size*dt*4*0.17)*data_R_c{n}.^2 + (n-1)*0.05,'r.')

    end
%     set(hs,'XLabel','Radius of Gyration (\mum)','YLabel','Diffusion coefficient (\mum^2.s^{-1})')

    
    
%     plot(data_R_c{n},data_D_corr{n}+(n-1)/2,'o')
    
end

% name
% % data_D_lin_rest = data_D_lin(data_D_lin{:}<0.2);
% data_D_corr_nopx_rest = data_D_corr_nopx(data_D_corr_nopx{:}<0.2);

subplot(1,4,1)
title('Diffusion calculated from linear fitting')
h = distributionPlot(data_D_lin,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
    'histOpt',2,'histOri','right','xyOri','flipped','showMM',4);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end
hold all
h = distributionPlot(data_D_corr_nopx,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
    'histOpt',2,'histOri','right','xyOri','flipped','color','r','showMM',0);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end

subplot(1,4,2)
title('Diffusion from linear sliding fitting')
h = distributionPlot(data_D_corr_nopx,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
    'histOpt',0,'histOri','right','xyOri','flipped','color','r','showMM',4);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end
hold all
h = distributionPlot(data_D_corr,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
    'histOpt',0,'histOri','right','xyOri','flipped','color','k','showMM',0);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end

% subplot(1,4,3)
% title('Diffusion outside of the corral calculated from corraled fitting')
% h = distributionPlot(data_D_corr_out,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
%     'histOpt',0,'histOri','right','xyOri','flipped','showMM',0);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end
% hold all
% h = distributionPlot(data_D_corr_nopx,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
%     'histOpt',0,'histOri','right','xyOri','flipped','color','r','showMM',0);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end

subplot(1,4,3)
title('Radius of gyration')
h = distributionPlot(data_R_c_nopx,'xName',name,'yLabel','Radius of confinement (\mum)',...
    'histOpt',0,'histOri','right','xyOri','flipped','color','r','showMM',0);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end
hold all
h = distributionPlot(data_R_c,'xName',name,'yLabel','Radius of confinement (\mum)',...
    'histOpt',0,'histOri','right','xyOri','flipped','showMM',4);
% for i = 1:length(filename)
%     set(h{1}(i),'FaceAlpha',0.5)
% end




end