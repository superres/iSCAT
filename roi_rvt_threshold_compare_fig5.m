clc;clear;close all; warning off;

t=[0.06,0.07,0.1];
pages=1; r=1:25;
a0=readTIF('compare_median10-1.tif');
titles1={'Ground truth','Raw','Deblurred'};
titles2={'Ground truth RVT','Raw RVT','Deblurred RVT'};
for i=1:3
    getParticles(a0(:,:,i),r,t(i),15,5,i,titles1{i},titles2{i});
end


function getParticles(a,r,t,idx,title1,title2)
    b=RVT(a,r);
    b=normalMaxMin(b); b(b<t)=0;
    subplot(2,3,idx); 
    imshow(a,[]); title(title1); 
    pl=subplot(2,3,idx+3); imshow(b,[]); title(title2); 
end
