function [x,y] = corral_diffusion(a,b,D,p,x0,y0,c,num)

% close all
n = num; % number of steps

% a = 4; % corral size
% D = 0.1; % diffusion coefficient
% dt = 1; % delta time

% p = 0.01;

x = zeros(1,n);
y = zeros(1,n);

x(1) = x0;
y(1) = y0;

for i = 1:n-1
   %% Cases
    if x(i) > (a/2+c(1)) || x(i) < (-a/2+c(1)) % exit from the right or left
        proba = rand;
        if proba > 1-p
            xv = x(1:i);
            yv = y(1:i);
            idx = i;
            if x(i) > (a/2+c(1))
                exit = 1;
            elseif x(i) < (-a/2+c(1))
                exit = 2;
            end
%             plot(xv,yv)
%             hold all
%             plot(x(i),y(i),'ro')
%             axis equal
            return
        end
        x(i+1) = x(i-1);
        if x(i) > c(1)
            dy = (a/2 + c(1) - x(i-1))*(y(i)-y(i-1))/(x(i)-x(i-1));
            x(i) = a/2 + c(1);
            y(i) = y(i-1) + dy;
        elseif x(i) < c(1)
            dy = (-a/2 + c(1) - x(i-1))*(y(i)-y(i-1))/(x(i)-x(i-1));
            x(i) = -a/2 + c(1);
            y(i) = y(i-1) + dy;
        end
        y(i+1) = y(i-1) + 2*dy;
    elseif y(i) > (b/2+c(2)) || y(i) < (-b/2+c(2)) % exit from the right or left
        proba = rand;
        if proba > 1-p
            xv = x(1:i);
            yv = y(1:i);
            idx = i;
            if y(i) > (b/2+c(2))
                exit = 3;
            elseif y(i) < (-b/2+c(2))
                exit = 4;
            end
%             plot(xv,yv)
%             hold all
%             plot(x(i),y(i),'ro')
%             axis equal
            return
        end
        y(i+1) = y(i-1);
        if y(i) > c(2)
            dx = (b/2 +c(2) - y(i-1))*(x(i)-x(i-1))/(y(i)-y(i-1));
            y(i) = b/2 + c(2);
            x(i) = x(i-1) + dx;
        elseif y(i) < c(2)
            dx = (-b/2 +c(2) - y(i-1))*(x(i)-x(i-1))/(y(i)-y(i-1));
            y(i) = -b/2 + c(2);
            x(i) = x(i-1) + dx;
        end
        x(i+1) = x(i-1) + 2*dx;
    else
        dr = randn;
        theta = randi([0 36000])*pi/18000;
   
        x(i+1) = x(i) + 2*D*dr*cos(theta);
        y(i+1) = y(i) + 2*D*dr*sin(theta);
    end
 
    
end
% 
plot(x,y)
axis equal
axis([-a/2 a/2 -b/2 b/2])


















end