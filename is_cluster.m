function [boolean] = is_cluster(center1,radius1,center2,radius2,thresh)

dist = sqrt( (center1(1)-center2(1))^2 + (center1(2)-center2(2))^2 );
mean_radius = mean([radius1;radius2]);

if dist >= thresh*mean_radius
    boolean = 1;
elseif dist < thresh*mean_radius
    boolean = 0;
end

end