function pull_result_data_linear

[filename,path] = uigetfile('.mat','multiselect','on');
cd(path)

D = [];


for i = 1:length(filename)
    
    result = struct();
    result = importdata(filename{i});

    D = [D,result.D'];

    
end

Din = D';

pull = struct('D',D);

uisave('pull')



end