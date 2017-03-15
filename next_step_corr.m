function [step_size,disp,C] = next_step_corr(x,y,step)

C = zeros(1,floor((length(x)-2)/step));
step_size = zeros(1,floor((length(x)-2)/step));
disp = zeros(1,floor((length(x)-2)/step));

for i = 1+step:step:length(x)-step
    
    C(i) = ( (x(i)-x(i-step))*(x(i+step)-x(i)) + (y(i)-y(i-step))*(y(i+step)-y(i)) ) /...
        sqrt( ( (x(i)-x(i-step))^2 + (y(i)-y(i-step))^2 ) );
    step_size(i) = sqrt( ( (x(i)-x(i-step))^2 + (y(i)-y(i-step))^2 ) );
    disp(i) = sqrt( ( (x(i)-x(1))^2 + (y(i)-y(1))^2 ) );
    
end


end