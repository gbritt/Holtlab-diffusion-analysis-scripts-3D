function pool_data_fit_lin_corr

[filename,path] = uigetfile('.mat','multiselect','on');
cd(path)

lin_fit = cell(1,1);
Dlin = [];

corr_fit = cell(1,2);
Dcorr = [];
Rc = [];


for i = 1:length(filename)
    
    res = importdata(filename{i});

    Dlin = [Dlin,res.lin.D_lin'];
    
    size_1 = size(res.corr,1);
    for j = 1:size_1
        Dcorr = [Dcorr;res.corr(j).D_corr{1}(:)];
        Rc = [Rc;res.corr(j).R_c{1}(:)];
    end
  
end

lin_fit = {Dlin};
corr_fit = {{Dcorr},{Rc}};

row_headings = {'D_corr','R_c'};
extract_corr = cell2struct(corr_fit,row_headings,2);
extract_lin = cell2struct(lin_fit,'D_lin',1);

pool = struct('lin',extract_lin,'corr',extract_corr);

uisave('pool')








end