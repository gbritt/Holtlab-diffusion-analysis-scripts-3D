function find_radius(x,y)

[n,der,center,radius,idx] = number_cluster(x,y,10,0.4);

figure;
hold all
t = linspace(0,2*pi);
for i = 1:n
    plot(x(idx==i),y(idx==i))
    plot(center(1,i)+radius(i)*cos(t),center(2,i)+radius(i)*sin(t),'r')
end
axis equal



c = 0;
        
for j = 1:n
            
    x_red = x(idx == j);
    y_red = y(idx == j);
    [n_red,der,center,radius,idx] = number_cluster(x_red,y_red,10,0.4); 
    for i_red = 1:n_red
        plot(x_red(idx==i_red),y_red(idx==i_red))
        plot(center(1,i_red)+radius(i_red)*cos(t),center(2,i_red)+radius(i_red)*sin(t),'r')
        c = c+1
    end
    
    
end






end