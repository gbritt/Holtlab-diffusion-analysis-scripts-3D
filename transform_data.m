function transform_data(saving_name)

[filename,path] = uigetfile('multiselect','off','.tif');
cd(path)
size_tif = size(imfinfo(filename),1);

I_BF = (imread(filename,'index',1));
imwrite((I_BF),strcat('BF_',saving_name,'.tif'));


m = size(I_BF,1);
n = size(I_BF,2);

I_TIRF = zeros(m,n);

for i = 2:size_tif
    I_TIRF = zeros(m,n);
    if i == 2
        I_TIRF = (imread(filename,'index',i));
        imwrite(I_TIRF,strcat('TIRF_',saving_name,'.tif'));
    else
        I_TIRF = (imread(filename,'index',i));
        imwrite(I_TIRF,strcat('TIRF_',saving_name,'.tif'),'WriteMode','Append');
    end
end




end