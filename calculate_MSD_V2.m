function [MSD,time] = calculate_MSD_V2(x,y,z,dt,conv)

tau_max = length(x);

if z == 0
    z = zeros(size(x));
end

MSD = zeros(1,tau_max);


for tau = 1:tau_max; 
   sqdisp_x = conv^2*((x(1+tau:end) - x(1:end-tau)).^2); 
   sqdisp_y = conv^2*((y(1+tau:end) - y(1:end-tau)).^2);
   sqdisp_z = conv^2*((z(1+tau:end) - z(1:end-tau)).^2);
   sqdisp = sqdisp_x + sqdisp_y + sqdisp_z;
   MSD(tau) = mean(sqdisp);
end

time = 0:dt:dt*(length(MSD)-1);

end