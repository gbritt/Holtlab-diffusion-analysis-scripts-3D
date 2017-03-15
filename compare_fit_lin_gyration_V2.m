function compare_fit_lin_gyration_V2(multiple,x_names,dt,steps_hist,size_step,n_max)

if multiple == 0
    [filename,path] = uigetfile('.mat');
    num = size(filename,1);
elseif multiple == 1
    [filename,path] = uigetfile('multiselect','on','.mat');
    num = size(filename,2);
end



cd(path)

figure('position',[50 300 2500 600]);

if multiple == 1
    data_R_c = cell(1,length(filename));
    data_D_lin = cell(1,length(filename));
    name = cell(1,length(filename));
    mean_D = zeros(1,length(filename));
    med_D = zeros(1,length(filename));
    err_mean_D = zeros(1,length(filename));
elseif multiple == 0
    data_R_c = [];
    data_D_lin = [];
    name = [];
    mean_D = 0;
    med_D = 0;
    err_mean_D = 0;
end

lambda = zeros(size_step,num);
Dalpha = zeros(1,num);
alpha = zeros(1,num);

for n  = 1:num
    
    if multiple == 1
        res = importdata(filename{n});
    elseif multiple == 0
        res = importdata(filename);
    end
    
    if multiple == 1
        data_D_lin{n} = [res.lin.D_lin{1}(:)];
    elseif multiple == 0
        data_D_lin = [res.lin.D_lin];
        [pD,nbinD] = hist(data_D_lin(data_D_lin > 0),20);
%         yyD = fit(nbinD(1:10)',pD(1:10)'/trapz(nbinD(1:10),pD(1:10)),'a*exp(-a*x)');
%         D0 = yyD.a;
    end
    
    subplot(1,7,2)
    title('Ensemble averaged MSD')
    mean_D(n) = mean(res.lin.D_lin{1}(:));
    med_D(n) = median(res.lin.D_lin{1}(:));
    err_mean_D(n) = std(res.lin.D_lin{1}(:))/sqrt(numel(res.lin.D_lin{1}(:))-1);
    t = 0:dt:(n_max-1)*dt;
    hold all
    plot(t,res.MSD.MSD_ens{1}(1:n_max),'LineWidth',2)
%     yyMSD = fit(t(2:20)',res.MSD.MSD_ens{1}(2:20)','power1','Startpoint',[mean_D 1]);
%     Dalpha(n) = yyMSD.a;
%     alpha(n) = yyMSD.b;
%     plot(t,4*mean_D*t,'k')
%     plot(yyMSD,'r')
    xlabel('Time (s)')
    ylabel('MSD (\mum^2)')
    set(gca,'XScale','log','YScale','log')
    
    subplot(1,7,4)
    title('Distribution of displacement')
    for i = 1:length(steps_hist)
        step = res.step(steps_hist(i)).step{1}(res.step(steps_hist(i)).step{1}(:) > 0);
        [p,nbin] = hist(step,20);
        hold all
        plot(nbin,(10^(n-1))*p/trapz(nbin,p),'o-','LineWidth',2)
%         plot(nbin,p/trapz(nbin,p),'o-','LineWidth',2)

    end
    xlabel('\deltax (\mum)')
    ylabel('\phi(\deltax)')
    set(gca,'YScale','log')
    
    subplot(1,7,5)
    title('Distribution of angle')
    for i = 1:length(steps_hist)
        step = res.ori(steps_hist(i)).ori{1}(res.ori(steps_hist(i)).ori{1}(:) ~= 0);
        [p,nbin] = hist(acos(step),50);
% rose(acos(step),50)
        hold all
        plot(nbin(4:end-4),0.5*(n-1)+p(4:end-4)/trapz(nbin,p),'o-','LineWidth',2)
% plot(nbin(4:end-4),p(4:end-4)/trapz(nbin,p),'o-','LineWidth',2)
%         plot([-1 1],[0.5+(n-1) 0.5+(n-1)],'k')
    end
    xlabel('\theta')
    ylabel('\phi(\theta)')
%     set(gca,'YScale','log')

    subplot(1,7,6)
    title('Next-step correlation')
    hold all
    long = min(size(res.next.step_size{1},1),size(res.next.corr{1},1));
    [x,y,err] = binning_2(res.next.step_size{1}(1:long),res.next.corr{1}(1:long),20);
    errorbar(x,y,err,'o-','LineWidth',2)
%     yy = fit(x(1:5)',y(1:5),'a*x','weight',err(1:5));
% %     plot(yy,'r')
%     lambda = yy.a
    xlabel('Step size (\mum)')
    ylabel('Correlation next step (\mum)')
    
    subplot(1,7,7)
    title('Distribution of ergodicity parameter')
    for i = 1:length(steps_hist)
        ergo = res.ergo.ergo(:,steps_hist(i)+1);
        [p,nbin] = hist(ergo,15);
        hold all
        plot(nbin,(10^(n-1))*p/trapz(nbin,p),'o-','LineWidth',2)
    end
%     plot(nbinD/mean_D,pD/trapz(nbinD,pD),'k','LineWidth',2)
    xlabel('\xi')
    ylabel('\phi(\xi)')
    set(gca,'YScale','log')
%     
%     if multiple == 1
%         data_D_lin{n} = res.lin.D_lin;
%     elseif multiple == 0
%         data_D_lin = res.lin.D_lin;
%     end
   
    
    if isempty(x_names) == 1
        if multiple == 1
            name{n} = filename{n};
        elseif multiple == 0
            name = filename;
        end
    else
        name = x_names;
    end


    size_1 = size(res.corr,1); % number of trajectories that have been analyzed with the corraled fitting
    for i = 1:size_1
        if multiple == 1
            data_R_c{n} = [data_R_c{n} res.corr(i).R_c{:}];
        elseif multiple == 0
            data_R_c = [data_R_c;res.corr(i).R_c{1}(:)];
        end
    end
    
    
%     Fitting the results of the step size distribution
    
%     for i = 1:size_step
%         step = res.step(i).step{1}(res.step(i).step{1}(:) > 0);
%         [p,nbin] = hist(step,20);
%         pnorm = p/trapz(nbin,p);
%         yy = fit(nbin(3:end)',pnorm(3:end)','exp1');
%         lambda(i,n) = -1/yy.b;
%     end

end




vec_color = zeros(3,7);
vec_color(:,1) = [0 0.45 0.74];
vec_color(:,2) = [0.85 0.33 0.1];
vec_color(:,3) = [0.93 0.69 0.13];
vec_color(:,4) = [0.49 0.18 0.56];
vec_color(:,5) = [0.47 0.67 0.19];
vec_color(:,6) = [0.3 0.75 0.93];
vec_color(:,7) = [0.64 0.08 0.18];

subplot(1,7,1)
title('Diffusion calculated from linear fitting')
if multiple == 1
%     h = distributionPlot(data_D_lin,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
%         'histOpt',2,'histOri','right','xyOri','flipped','showMM',4);
    h = distributionPlot(data_D_lin,'xName',name,'yLabel','Diffusion (\mum^2.s^{-1})',...
        'histOpt',2,'histOri','right','xyOri','flipped','showMM',4);
    for i = 1:length(h{1}(:))
        if i > 7
            set(h{1}(i),'FaceColor',vec_color(:,i-7));
        else
            set(h{1}(i),'FaceColor',vec_color(:,i));
        end
    end
%     
%     x1 = data_D_lin{1}(:);
%     x2 = data_D_lin{2}(:);
    
%     [a,b,c] = kstest2(x1,x2)
elseif multiple == 0
    plot(nbinD,pD/trapz(nbinD,pD),'bo-','LineWidth',2)
    hold all
    plot(yyD,'r')
end

subplot(1,7,3)
title('Radius of gyration')
h = distributionPlot(data_R_c,'xName',name,'yLabel','Radius of gyration (\mum)',...
    'histOpt',0,'histOri','right','xyOri','flipped','color','k','showMM',0);
    for i = 1:length(h{1}(:))
        if i > 7
            set(h{1}(i),'FaceColor',vec_color(:,i-7));
        else
            set(h{1}(i),'FaceColor',vec_color(:,i));
        end
    end

    
    figure
errorbar(1:1:length(filename),mean_D,err_mean_D,'bo')
hold all
plot(1:1:length(filename),med_D,'ro')


% figure
% hold all
% time_step = dt:dt:(size_step)*dt;
% exponent = zeros(1,num);
% diff = zeros(1,num);
% for i = 1:num
%     yy = fit(time_step()',lambda(:,i),'power1');
%     exponent(i) = yy.b;
%     diff(i) = yy.a^2;
%     plot(time_step,lambda(:,i),'o-','LineWidth',2)
%     plot(yy,'r')
%     
%     plot(time_step,sqrt(Dalpha(i)*time_step.^(alpha(i)/2)),'k')
% end
% set(gca,'YScale','log','XScale','log')
% 
% figure
% plot(exponent,'bo-','LineWidth',2)
% figure
% plot(diff,'bo-','LineWidth',2)

% figure
% plot(nbinD,pD/trapz(nbinD,pD),'bo-','LineWidth',2)
% yy = fit(nbinD(1:10)',pD(1:10)'/trapz(nbinD(1:10),pD(1:10)),'a*exp(-a*x)')
% hold all
% plot(yy,'r')
% set(gca,'YScale','log')








end