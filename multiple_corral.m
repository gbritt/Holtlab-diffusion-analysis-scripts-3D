function multiple_corral

% close all

D = 0.2;
p = -1;
a0 = 1;
b0 = 1000;

num = 1000;

n_part = 1;

S = [];
for i = 1:length(n_part)
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
% h = waitbar(0,'Simulating...');
%% Loop on multiple particles

for n = 1:n_part


%% Initiation of variables
x0 = 0;
y0 = 0;

c = [0,0];

x = [];
y = [];
idx = [];

a = abs(a0+randn);
b = abs(b0+randn);
[xv,yv] = corral_diffusion(a,b,D,p,x0,y0,[x0,y0],num);

x = [x,xv];
y = [y,yv];
% idx = [idx,idx_int];

%% Run for multiple hop 
% for i = 1:10   
%     
% %     c = [0,0];
%     if exit == 1
%         c(1) = c(1) + 0.5*a;
%         a = abs(a0+randn);
%         c(1) = c(1) + 0.5*a;
%         c(2) = y(end);
%     elseif exit == 2
%         c(1) = c(1) - 0.5*a;
%         a = abs(a0+randn);
%         c(1) = c(1) - 0.5*a;
%         c(2) = y(end);
%     elseif exit == 3
%         c(2) = c(2) + 0.5*a;
%         a = abs(a0+randn);
%         c(2) = c(2) + 0.5*a;
%         c(1) = x(end);
%     elseif exit == 4
%         c(2) = c(2) - 0.5*a;
%         a = abs(a0+randn);
%         c(2) = c(2) - 0.5*a;
%         c(1) = x(end);
%     end
%     [xv,yv,idx_int,exit] = corral_diffusion(a,D,p,x(end),y(end),c);
%     x = [x,xv];
%     y = [y,yv];
% %     idx = [idx,idx(i)+idx_int];
%     
%     
% end

% figure
% hold all
% plot(x(1:idx(1)),y(1:idx(1)))
% for i = 1:length(idx)-1
%     plot(x(idx(i):idx(i+1)),y(idx(i):idx(i+1)))
%     
% end
% axis equal
MSD = calculate_MSD(x,y,0,1,1);
time = 0:1:length(MSD);



% figure
% plot(MSD)
% hold all
% v = 0:1:length(MSD);
% MSD_th = 4*D*v;
% plot(MSD_th,'k')
% % for i = 1:length(idx)
% %     plot(idx(i),MSD(idx(i)),'ro')
% %     
% % end

result(n).tracking = struct('time',time',...
        'x',x',...
        'y',y',...
        'MSD',MSD');
    
% waitbar(n/n_part)

end

disp_res = cell(1,num);
disp_ori = cell(1,num);

displ = [];
ori = [];

for i = 1:num
    for j = 1:n_part  
        displ = [displ;displacement(result(j).tracking.x,result(j).tracking.y,i)'];
        ori = [ori;disp_corr(result(j).tracking.x,result(j).tracking.y,i)'];
    end
    disp_res{i} = {displ};
    disp_ori{i} = {ori};
end

stat = cell2struct(disp_ori,'ori',1);
extract = struct('stat',stat);

% result.ori = struct('ori',disp_ori);
% result.displ = struct('disp',disp_res);

uisave('stat')

end