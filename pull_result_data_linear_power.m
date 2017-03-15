function pull_result_data_linear_power

[filename,path] = uigetfile('.mat','multiselect','on');
cd(path)

D = [];
Dalpha = [];
alpha = [];


for i = 1:length(filename)
    
    result = struct();
    result = importdata(filename{i});

    D = [D,result.D'];
    Dalpha = [Dalpha,result.Dalpha'];
    alpha = [alpha,result.alpha'];

    
end

D = D';

pull = struct('D',D,'Dalpha',Dalpha,'alpha',alpha);

uisave('pull')



end