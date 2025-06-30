clc; close all;
clear;

[psf,otf]=getPsf(520,0.9,0.058,0.9,0.9,128);
a0=readTIF('ipsf1.tif');
r=64; l=size(a0,1)*0.5;
a0=a0(l-r:l+r+1,l-r:l+r+1);
a=0.5*(a0(r,:)+a0(r+1,:));
a1=getFFT(a0);
r=size(a1,1)*0.5;
a2=0.5*(a1(r,:)+a1(r+1,:));
r=size(otf,1)*0.5;
a3=0.5*(otf(r,:)+otf(r+1,:));
xx1=1:max(size(a3));
xx=1:size(a0,2);
m=2; n=2;
subplot(m,n,1); imshow(a0,[]);
subplot(m,n,2); plot(xx,a);
subplot(m,n,3); imshow(a1,[]);
subplot(m,n,4); hold on;
plot(xx,a2);
plot(xx1,a3); hold off;

function b=getFFT(a) 
    b=log(1+abs(fftshift(fft2(a))));
end
function [psf,otf]=getPsf(lambda,NA,pixelsize,attStrength,attFWHM,imgSize)
    % example
    % lambda = 525; % 单位：纳米 % 发射波长
    % NA = 1.49; pixelsize = 65e-3; % 单位：微米
    % RL = 5;
    % attStrength = 0.9; % 默认值
    % attFWHM=1.0;
    % imgSize=512;
    % 
    % [psf,otf]=getPsf(lambda,NA,pixelsize,attStrength,attFWHM,imgSize);

    cutoff=1000/(0.5*lambda/NA);
    % imgSize=size(Iraw,1);
    cyclesPerMicron=1/(imgSize*pixelsize);
    sampleLateral=ceil(cutoff/cyclesPerMicron)+1;
    
    for I=1:sampleLateral
        v=abs(I-1)/sampleLateral;
        vals(I)=(1/pi)*(2*acos(v)-sin(2*acos(v)))*power(1,v);
        dist=abs(I-1)*cyclesPerMicron;
        valsOnlyAtt(I)=1-attStrength*exp(-power(dist,2)/(power(0.5*attFWHM,2)));
        valsAtt(I)=vals(I)*valsOnlyAtt(I);
    end
    
    y=1:imgSize;
    x=1:imgSize;
    [x,y]=meshgrid(x,y);
    rad=hypot(y-imgSize/2-1,x-imgSize/2-1);
    cycl=rad.*cyclesPerMicron;
    onlyatt=1-attStrength*(exp(-power(cycl,2)/(power(0.5*attFWHM,2))));
    
    mask=cycl>cutoff;
    cycl(mask)=0;
    
    pos=cycl./cyclesPerMicron;
    cpos=pos+1;
    
    lpos=floor(cpos);
    hpos=ceil(cpos);
    f=cpos-lpos;
    
    retl=vals(lpos).*(1-f);
    reth=vals(hpos).*f;
    otf=retl+reth;
    
    mask1=ceil(cpos)>sampleLateral;
    otf(mask1)=0;
    otf(mask)=0;
    otfatt=otf.*onlyatt;
    
    psf = ifftn(otf);
    psf = circshift(psf,floor([imgSize,imgSize]/2));
       
    psf = abs(psf);
end


