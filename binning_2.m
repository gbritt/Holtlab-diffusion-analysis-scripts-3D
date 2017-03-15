function [x,y,err] = binning_2(xe,ye,n_bin)

d = max(xe)-min(xe);
l = length(xe);

bin_size = d/n_bin;

x = min(xe):bin_size:max(xe);
y = zeros(length(x),1);
err = zeros(length(x),1);



for w = 1:n_bin
    calc = zeros(l,1);
    for v = 1:l
        if (xe(v) < (min(xe) + w*bin_size)) && (xe(v) > (min(xe) + (w-1) * bin_size))
            calc(v) = ye(v);
        end
    end
%     calc(1)
%     pause
        ind = find(calc ~= 0);
        calc_int = calc(ind);
        y(w) = mean(calc_int);
        err(w) = std(calc_int)/sqrt(length(ind)-1);

end



