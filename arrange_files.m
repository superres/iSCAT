clc;clear;close all; warning off;

fs0='raw_'; %fs0='E:\Papers\my own\iOTF\rebuttal\simulated\raw_';
%%stack files of 100 Mb were separated into 5 stacks 20 Mb to avoid 25Mb
%%upload limitation of GitHub
for i=1:9
    fs1=sprintf('%d.tif',i);
    fs1=[fs0,fs1];
    a=readTIF(fs1);
    for j=1:5
        a1=a(:,:,j*20-19:j*20);
        fs2=sprintf('%d.%d.tif',i,j);
        fs2=[fs0,fs2];
        writeTIF(fs2,a1);
    end
end

%%5 separated stacks of 20 Mb were added up as one stach (100 Mb)
for i=1:9
    for j=1:5
        fs2=sprintf('%d.%d.tif',i,j);
        fs2=[fs0,fs2];
        a(:,:,j*20-19:j*20)=readTIF(fs2);
    end
    fs1=sprintf('%d.tif',i);
    fs1=[fs0,fs1];
    writeTIF(fs1,a);
end

function file=readTIF(filename)   
    sizef=size(imfinfo(filename),1);
    file = imread(filename, 1); 
    for idx = 2 : sizef
        tiffTemp = imread(filename, idx);
        file = cat(3 , file, tiffTemp);
    end 
    file=double(file);
end
