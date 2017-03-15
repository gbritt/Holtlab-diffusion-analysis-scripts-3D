function pull_result_data

[filename,path] = uigetfile('.mat','multiselect','on');
cd(path)

Din = [];
Dout = [];
tres = [];
Rc = [];
alpha = [];

for i = 1:length(filename)
    
    result = struct();
    result = importdata(filename{i});

    Din = [Din,result.Din'];
    Dout = [Din,result.Dout'];
    tres = [tres,result.Tres'];
    Rc = [Rc,result.Rc'];
    alpha = [alpha,result.alpha'];
    
end

Din = Din';
Dout = Dout';
tres = tres';
Rc = Rc';
alpha = alpha';

pull = struct('Din',Din,'Dout',Dout,'Rc',Rc,'Tres',tres,'alpha',alpha);

uisave('pull')



end