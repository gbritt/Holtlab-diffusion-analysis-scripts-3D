function gyration

[filename,path] = uigetfile('multiselect','off','.mat');
cd(path)

res = importdata(filename);

size_1 = size(res.corr,1);

figure
hold all

for i = 1:size_1
       
    
%        data_D_corr{n} = [data_D_corr{n} res.corr(i).D_corr{:}];
       plot(res.corr(i).R_c{:},res.corr(i).D_corr{:},'o')
       
end








end