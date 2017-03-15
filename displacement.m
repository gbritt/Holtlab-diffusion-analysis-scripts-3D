function displ = displacement(x,y,z,step)
%may break if no z

displ = zeros(1,length(x));


for i = 1:step:length(x)-step
   displ(i) = sqrt( ( (x(i+step)-x(i))^2 + (y(i+step)-y(i))^2 + ...
       (z(i+step)-z(i))^2) );
end





end