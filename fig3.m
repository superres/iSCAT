clc; close all;
clear;

tic;
%% 1. prerecognised particle from RVT and DoM plugin in imageJ
%% 2. the output data format of RVT is 32 bit, while the available format for DoM is 8 or 16 bit
%%    so there was a data convertion from 32 bit to 16 bit in imageJ (Image->Type->16-bit)
p00=csvread('points_from_RVT_DOM_PIXEL2_SNR7.csv',1,1);
p0=round(p00(:,1:2));
%% 3. apply 2D Gaussian fitting to get accurate localization and peak intensity of particles,
%%    an iteration method was introduced and the initial parameters were determined from the result of RVT and DoM 
fs0='MED_raw_10.tif';
p=getRefsPosi(fs0,15,4,p0);
%% 4. generate isolated particle image in the object space,
%%    the localizations and intensity of isolated pixels were fitted above
a0=zeros(512,512,100);
for i=1:size(p,1)
    a0(p(i,2),p(i,1),:)=p(i,3);
end
%% 5. load iPSF model as convolution kernel (iPSF generated with the python code https://github.com/SandoghdarLab/PiSCAT) 
ipsf00=readTIF('ipsf.tif');
%% 6. the simulated signal image was obtained by applying convolution of the iPSF and particle image in the object space
iotf0=fftshift(fft2(ipsf00));
a0 = real(fftshift(ifft2( fft2(a0).*fftshift(iotf0))))+1;
a12=a0(:,:,1);
% imshow(a12,[]);
V = 0.0001;
%% 7. the background intensity of simulated signal image was normalized 1
%%    this simulated signal image stack were taken as the ground truth images in simulation
a0=getBackgroundNormalized1(a12,V);
a10=a0(:,:,1);
a11=getNotched1(a10,p,512,0.08,520,0.36,27);

mm=2; nn=3;
% subplot(mm,nn,1); imshow(a10,[]);
% subplot(mm,nn,2); imshow(a11,[]);
% subplot(mm,nn,3); imshow(a12,[]);

l=18;
nf=23;
xx=1:2*l+1; xx=xx*0.058;
[a20,a30]=getCutLine(a10,p(nf,:),l);
[a21,a31]=getCutLine(a11,p(nf,:),l);
[a22,a32]=getCutLine(a12,p(nf,:),l);

subplot(mm,nn,1); imshow(a30,[]); colorbar;
subplot(mm,nn,2); imshow(a31,[]); colorbar;
subplot(mm,nn,3); imshow(a32,[]); colorbar;

subplot(mm,nn,1+nn); hold on;
plot(xx,a20,'b'); plot(xx,a22,'k'); hold off;
subplot(mm,nn,2+nn); hold on;
plot(xx,a21,'r'); plot(xx,a22,'k'); hold off;

toc;

tic;
vs=1:9; vs=vs*0.00006;
for i=1:9
    a10 = imnoise(a12, 'gaussian', 0, vs(i));
    a11=getNotched1(a10,p,512,0.08,520,0.36,15);
    ssims(1,i)=ssim(a10,a12);
    ssims(2,i)=ssim(a11,a12);
end

xxx=1:9; xxx=xxx*80;subplot(mm,nn,3+nn); hold on;
scatter(vs,ssims(1,:),'filled'); scatter(vs,ssims(2,:),'filled'); hold off;
legend('Raw','Deblurred');
toc;

function [line_intensity,a1]=getCutLine(a,p,r)
    x=p(1)-r:p(1)+r;
    y=p(2)-r:p(2)+r;
    a1=a(y,x);
    line_intensity=a(p(2),x);
end



function a=getBackgroundNormalized1(a,V)
    % V = 0.0003; % 噪声方差
    a = imnoise(a, 'gaussian', 0, V);
end

function p=getRefsPosi(fs0,l,l0,p)
    b=readTIF(fs0);
    for ii=1:size(p,1)
        x1=max(1,p(ii,1)-l0):min(size(b,1),p(ii,1)+l0);
        y1=max(1,p(ii,2)-l0):min(size(b,2),p(ii,2)+l0);
        rs=eps+1-b(y1,x1); rs(rs<0)=0;
        x=max(1,p(ii,1)-l):min(size(b,1),p(ii,1)+l);
        y=max(1,p(ii,2)-l):min(size(b,2),p(ii,2)+l);
        rois{ii}=zeros(length(y),length(x));
        rois{ii}(min(y1)-min(y)+1:max(y1)-min(y)+1,min(x1)-min(x)+1:max(x1)-min(x)+1)=rs;
        ps=fit2DGaussian(rois{ii},[p(ii,1)-x(1),p(ii,2)-y(1)]);
        p(ii,3)=ps(1);
    end
    p=sortrows(p,3);
end

function p=fit2DGaussian(a,ps)
    x=1:size(a,2); y=1:size(a,1);
    a=a+eps;
    [X Y]=meshgrid(x,y);
    p0=[0.5,ps(1),ps(2),3,3];
    options = optimoptions('fmincon','MaxIter',500,'Display','off');
    f=@(p) sum(sum(abs(a-p(1)*exp(-0.5*(X-p(2)).*(X-p(2))./(p(4)*p(4))-0.5*(Y-p(3)).*(Y-p(3))./(p(5)*p(5))))));
    p=fmincon(f,p0,[],[],[],[],[0.0000001,ps(1)-1,ps(2)-1,1,1],[1,ps(2)+1,ps(2)+1,6,6],[],options);
end

function e=getNotched1(a,p,npixel,pixelsize,wavelength,na,r)   
    %% a parameter in line 115 was introduced to compensate for intensity loss of the
    %% Fourier domain filter, since any components in Fourier domain
    %% contribute to the particle intensity and back ground noise intensity
    %% an optimal compensate parameter was determined by reaching the highest
    %% structure similarity index measure (SSIM) score with the ground truth
    m_compensate=14;
    for i=1:size(a,3)
        %% apply the filter in the frequency domain
        b=fftshift(fft2(a(:,:,i)));
        mask=NotchFilter( npixel,pixelsize,wavelength,na,0,0 );
        c=b.*mask;
        d=real(ifft2(fftshift(c)));
        %% compensate the intensity with parameter in line 115
        [amean,~]=getBackgrounds(d,p,r);
        e(:,:,i)=1+m_compensate*(d-amean);
    end
end
function otfAtt= NotchFilter( imgSize,micronsPerPixel,lambda,NA,kx,ky )
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

    Mask=double(rad<cutoff/cyclesPerMicron+1);
    otfAtt=otfAtt.*Mask;
end