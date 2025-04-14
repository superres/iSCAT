clc;clear;close all; warning off;

%%% this script was to illustrate the input raw (blurred) image, noise image
%%% and deblurred image in the spatial and frequency domain.
%%% the filter form was shown both for the original notch filter introduced
%%% from HIFI-SIM and the iOTF limitted notch filter.

%% subplot dimension
m=2; n=3; 

%%%raw image in the spatial and the frequency domain
a0=readTIF('E:\Papers\my own\iOTF\rebuttal\simulated\raw6.tif');
b=fftshift(fft2(a0));
ss=size(a0,1);
figure;
cm1=subplot(m,n,1); imshow(a0,[]); colormap(cm1,'gray'); 
title('Blurred');
c1=colorbar; 
maxa=max(max(a0(:,:,1)));mina=min(min(a0(:,:,1))); ll=maxa-mina;
mina=ceil(100*mina+ll)/100; maxa=floor(100*maxa-ll)/100; ll=maxa-mina;
set(c1,'YTick',mina:ll:maxa);
set(c1,'FontSize',16);
cm2=subplot(m,n,n+1); imshow(log(1+abs(b)),[]); colormap(cm2,'parula'); 
title('FT(blurred)');
c4=colorbar;
maxa=max(max(log(1+abs(b))));mina=min(min(log(1+abs(b)))); ll=maxa-mina;
mina=ceil(100*mina+ll)/100; maxa=floor(100*maxa-ll)/100; ll=maxa-mina;
set(c4,'YTick',mina:ll:maxa);
set(c4,'FontSize',16);
%%%raw image in the spatial and the frequency domain

%%%notched iotf
mask01=NotchFilters( ss,0.08,520,0.29,0,0 );
e=b.*mask01;
f=real(ifft2(fftshift(e)));
f=1+4*(f-mean(mean(f)));

%% deblurred image in the spatial and the frequency domain
cm3=subplot(m,n,3); imshow(f,[]); colormap(cm3,'gray'); 
title('Deblurred');
c3=colorbar; 
maxa=max(max(f));mina=min(min(f)); ll=maxa-mina;
mina=ceil(100*mina+ll)/100; maxa=floor(100*maxa-ll)/100; ll=maxa-mina;
set(c3,'YTick',mina:ll:maxa);
set(c3,'FontSize',16);
cm4=subplot(m,n,n+3); imshow(log(1+abs(e)),[]); colormap(cm4,'parula'); 
title('FT(deblurred)');
c6=colorbar;
maxa=max(max(log(1+abs(e))));mina=min(min(log(1+abs(e)))); ll=maxa-mina;
mina=ceil(100*mina+ll)/100; maxa=floor(100*maxa-ll)/100; ll=maxa-mina;
set(c6,'YTick',mina:ll:maxa);
set(c6,'FontSize',16);
%% deblurred image in the spatial and the frequency domain

mask1=zeros(ss,ss);
mask1(mask01==0)=1;
c=b.*mask1;
d=abs(ifft2(fftshift(c)));

%% noise in the spatial and the frequency domain
cm5=subplot(m,n,2); imshow(d,[]); colormap(cm5,'gray'); 
title('Background noise');
c2=colorbar; 
maxa=max(max(d));mina=min(min(d));ll=maxa-mina;
mina=ceil(100*mina+ll)/100; maxa=floor(100*maxa-ll)/100; ll=maxa-mina;
set(c2,'YTick',mina:ll:maxa);
set(c2,'FontSize',16);
cm6=subplot(m,n,n+2); imshow(log(1+abs(c)),[]); colormap(cm6,'parula'); 
title('FT(noise)');
c5=colorbar;
maxa=max(max(log(1+abs(c))));mina=min(min(log(1+abs(c)))); ll=maxa-mina;
mina=ceil(100*mina+ll)/100; maxa=floor(100*maxa-ll)/100; ll=maxa-mina;
set(c5,'YTick',mina:ll:maxa);
set(c5,'FontSize',16);
%% noise in the spatial and the frequency domain

%% notch filter and iOTF limitted notch filter
figure;
mask00=NotchFilters( 128,0.08,520,0.9,0,0 );
cm7=subplot(1,2,1); imshow(mask00,[]); colormap(cm7,'gray'); 
title('Notch filter');
cm8=subplot(1,2,2); imshow(mask01,[]); colormap(cm8,'gray'); 
title('iOTF boundary limitted notch filter');
mmax(1)=max(max(mask00)); mmin(1)=min(min(mask00));
mmax(2)=max(max(mask01)); mmin(2)=min(min(mask01));
%% notch filter and iOTF limitted notch filter

function otfAtt= NotchFilters( imgSize,micronsPerPixel,lambda,NA,kx,ky )
    % demo parameters micronsPerPixel=0.07;NPixel=512;lambda=520;NA=0.4;
    cutoff=1000/(0.5*lambda/NA);
    cyclesPerMicron=1/(imgSize*micronsPerPixel);
    w=imgSize;
    h=imgSize;
    siz=[h w];
    cnt=siz/2+1;
    kx=kx+cnt(2);
    ky=ky+cnt(1);

    y=1:h;
    x=1:w;
    [x,y]=meshgrid(x,y);
    rad=hypot(y-ky,x-kx);
    cycl=rad.*cyclesPerMicron;
    otfAtt=(1-exp(-power(cycl,4)/(2*power(cutoff,4))));

    cnt=[imgSize/2+1,imgSize/2+1];
    [x,y]=meshgrid(1:imgSize,1:imgSize);
    rad=sqrt((y-cnt(1)).^2+(x-cnt(2)).^2);

    Mask=double(rad<1.0*(cutoff/cyclesPerMicron+1));
    otfAtt=otfAtt.*Mask;
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
