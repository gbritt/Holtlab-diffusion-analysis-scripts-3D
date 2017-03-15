function brownian_diffusion

slope = zeros(1,50);

for k = 1:50

n = 500;

% t = 0:0.1:199*0.1;
x = zeros(1,n);
y = zeros(1,n);


D = 10;
% D = gamrnd(2.5,0.01,1,100);
Rg = zeros(1,n);

% S = [];
% for i = 1:100
%     if i < 10
%         S = [S;strcat('part0000',num2str(i))];
%     elseif i >= 10 && i < 100
%         S = [S;strcat('part000',num2str(i))];
%     elseif i >= 100 && i < 999
%         S = [S;strcat('part00',num2str(i))];
%     elseif i >=1000 && i < 9999
%         S = [S;strcat('part0',num2str(i))];
%     elseif i >= 10000 && i < 99999
%         S = [S;strcat('part',num2str(i))];
%     end
% end
% 
% field = cellstr(S);
% 
% result = cell2struct(field','tracking',1);

for m = 1:length(D)

for i = 1:n-1
    
   theta = randi(360)*2*pi/360;
   step = abs(randn);
    
    x(i+1) = x(i) + D(m)*cos(theta);
    y(i+1) = y(i) + D(m)*sin(theta);
    
    

end
% 
MSD  = calculate_MSD(x,y,0,0.1,1);
% 
% result(m).tracking = struct('time',t',...
%         'x',x(m,:)',...
%         'y',y(m,:)',...
%         'MSD',MSD');

x = x;
y = y;

for i = 1:n-1
mean_x = mean(x(1:i));
mean_y = mean(y(1:i));
    
Rg(i) = sqrt(numel(x(1:i))^(-1)*(sum((x(1:i)-mean_x).^2 + (y(1:i)-mean_y).^2)));
end

end



% saving_name = strcat('tracked_sim_3.mat');
% save(saving_name,'result')

% figure
% plot(x,y)
% axis equal
%
% hold all
% plot(MSD,Rg.^2,'ro')
ft = fittype('a*x');
% hold all
yy = fit(MSD(1:floor(n/5))',Rg(1:floor(n/5)).^2',ft);
slope(k) = yy.a;
% plot(yy)
end
% 
mean(slope)




end