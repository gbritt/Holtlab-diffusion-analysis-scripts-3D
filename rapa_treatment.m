function rapa_treatment

close all

dt = 0.01;

t = 0:dt:10;
kd = 1;
kd2 = 0.5*kd;
k = 0.1;
eps = 1;
kp = 2;

v0 = 20;
vv0 = (kd2-k)/kd2*v0;
D0 = (kd + kd2)*eps/kp;

D = zeros(1,length(t));
v = zeros(1,length(t));
vv = zeros(1,length(t));

D(1) = D0;
v(1) = v0;
vv(1) = vv0;

for i = 1:length(t)-1
    v(i+1) = v(i) + k*v(i)*dt;
    vv(i+1) = vv(i) + ( (kd2-k)*v(i) - kd2*vv(i) )*dt;
    D(i+1) = D(i) + ( (kd - kd2 + 2*k*v(i)/(v(i)-vv(i)))*D(i) - kp/eps*D(i).^2 )*dt;
    
end

figure
plot(t,D)

figure
plot(t,vv)
hold all
plot(t,v)

end