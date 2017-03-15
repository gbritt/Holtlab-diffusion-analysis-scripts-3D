function fit_multiple_data(D0,dt,conv,thresh_length)

[filename,path] = uigetfile('.mat','multiselect','on');
cd(path)


f_exp = fittype('a^2*(1-exp(-b*x))');
f_lin = fittype('a*x');

h = waitbar(0,'Fitting and extracting data...');

for m = 1:length(filename)
    
    result = struct();

    result = importdata(filename{m});

    num_tracks = size(result,1);
    if num_tracks == 1
       num_tracks = size(result,2);
    end
    
    isolate_idx = [];
    
    for j = 1:num_tracks
        if length(result(j).tracking.x) >= thresh_length
            isolate_idx = [isolate_idx;j];
        end
    end
    
    size_rest_data = length(isolate_idx);
    
    Din = [];
    Dout = [];
    Rc = [];
    tres = [];
    alpha = [];
    

    for j = 1:size_rest_data

    time = result(isolate_idx(j)).tracking.time;
    x = result(isolate_idx(j)).tracking.x;
    y = result(isolate_idx(j)).tracking.y;
    j
    [n,der,center,radius,idx] = number_cluster(time',x',y',5,0.8,0);
    idx_unique = unique(idx);
    for i = 1:n
            xr = x(idx==idx_unique(i));
            yr = y(idx==idx_unique(i));
            
            if length(xr) > 10
            
            [MSDr,tr] = calculate_MSD_V2(xr,yr,0,dt,conv);
            
            size_r = 4*round(length(xr)/5);
            if round(size_r/2) == size_r/2
                weight = [20*ones(1,size_r/2),ones(1,size_r/2)];
            else
                weight = [20*ones(1,round(size_r/2)),ones(1,floor(size_r/2))];
            end
            
            [yy_lin,gof_lin] = fit(tr(1:size_r)',MSDr(1:size_r)',f_lin,...
                'Startpoint',[D0],'weight',weight,'display','off');
            [yy_exp,gof_exp] = fit(tr(1:size_r)',MSDr(1:size_r)',f_exp,...
                'Robust','Bisquare',...
                'Weight',weight,'display','off');
 %                 'Startpoint',[(conv*radius(i)),D0/((conv*radius(i))^2)],... 
            lin_or_expo = gof_lin.rsquare/gof_exp.rsquare;
            
            if lin_or_expo < 0.7
                Din = [Din;yy_exp.a^2*yy_exp.b/4];
                Rc = [Rc;sqrt(5/2)*yy_exp.a];
                tres = [tres;length(xr)*dt];
            elseif lin_or_expo > 0.7
                yy_pow = fit(tr(2:size_r)',MSDr(2:size_r)','power1',...
                'Startpoint',[D0 1],'weight',weight(2:end),'display','off',...
                'Lower',[0 0]);
                Dout = [Dout;yy_lin.a];
                alpha = [alpha;yy_pow.b];
            end
            
            end
    end

    end
    
    saving_name = strcat('extracted_',filename{m});
    extract = struct('Din',Din,'Dout',Dout,'Rc',Rc,'Tres',tres,'alpha',alpha);
    save(saving_name,'extract')

    waitbar(m/length(filename))

end

disp('Data saved')
close(h)

end