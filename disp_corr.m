function [ori] = disp_corr(x,y,z,step)
%may crash if no z
ori = zeros(1,floor((length(x)/step-2)));
% step_next_size = zeros(1,floor((length(x)-2)));
% next_step = zeros(1,floor((length(x)-2)));

for i = 1+step:step:floor(length(x)-step)
   
    ori(i) = ( (x(i+step)-x(i))*(x(i)-x(i-step)) + (y(i+step)-y(i))*(y(i)-y(i-step)) + (z(i+step)-z(i))*(z(i)-z(i-step)) ) /...
        sqrt(  ( (x(i+step)-x(i))^2 + (y(i+step)-y(i))^2 + (z(i+step)-z(i))^2 ) * ...
        ( (x(i-step)-x(i))^2 + (y(i-step)-y(i))^2 + (z(i-step)-z(i))^2 )  );
    
%     step_size(i) = sqrt( ( (x(i+step)-x(i))^2 + (y(i+step)-y(i))^2 ) );
    
end

% for i = 2:lenght(x)-1
%     next_step(i) = ( (x(i+1)-x(i))*(x(i)-x(i-1)) + (y(i+1)-y(i))*(y(i)-y(i-1)) ) / sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 );
%     step_next_size(i) = sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 );
% end


end