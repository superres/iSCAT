clc;close all; clear;

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

%%% experimental curve comparison
fs0='raw\';
for i=1:9
    fs1=sprintf('raw_%d.tif',i);
    fs1=[fs0,fs1];
    exp_raw=readTIF(fs1);
    fs1=sprintf('raw_%d.tif',2);
    fs1=[fs0,fs1];
    exp_raw7=readTIF(fs1);
    exp_med=getMedian(exp_raw7);
    exp_deblurred=getNotched1(exp_raw,p,512,0.08,520,0.36,15);

    for j=1:9
        ssims2(:,i)=getCurve(exp_med,exp_raw(:,:,j),exp_deblurred(:,:,j),p(24,:),20,i,j);
    end
end
%%%%% fig.4(e)
xx=1:9; xx=80*xx;
point_size=200;
% figure; 
subplot(1,4,3:4); hold on;
for i=1:2
    scatter(xx,ssims2(i,:),point_size,'filled');
end
xlabel('Exposure time \mus'); ylabel('SSIM'); title('e','Position',[-75,1.01]);
ylim([min(min(ssims2))*0.99,1.01]);
legend('Blurred','Deblurred');
set(gca,'FontName','Helvetica','FontSize',16);
%%%%% fig.4(e)
toc;
%%% experimental curve comparison

function e=getNotched1(a,p,npixel,pixelsize,wavelength,na,r)   
    %% a parameter in line 115 was introduced to compensate for intensity loss of the
    %% Fourier domain filter, since any components in Fourier domain
    %% contribute to the particle intensity and back ground noise intensity
    %% an optimal compensate parameter was determined by reaching the highest
    %% structure similarity index measure (SSIM) score with the ground truth
    m_compensate=4;
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

function file=readTIF(filename)   
    sizef=size(imfinfo(filename),1);
    file = imread(filename, 1); 
    for idx = 2 : sizef
        tiffTemp = imread(filename, idx);
        file = cat(3 , file, tiffTemp);
    end 
    file=double(file);
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
function ssims=getCurve(a0,a1,a2,p,r0,idx,idy)
    line_width=3;
    x=p(1)-r0:p(1)+r0;
    y=p(2)-r0:p(2)+r0;
    a01=a0(y,x);
    
    [r,c]=find(a01==min(min(a01)));
    x=x(1)-1+c-r0:x(1)-1+c+r0;
    y=y(1)-1+r-r0:y(1)-1+r+r0;
    xx=1:2*r0+1; xx=xx*0.08;
    a(1,:)=a0(y(1)-1+r,x);
    a(2,:)=a1(y(1)-1+r,x);
    a(3,:)=a2(y(1)-1+r,x);

    if idx==1&&idy==1
        %%%%% fig. 3/4(a-d)
        figure; 
        subplot(2,4,1); imshow(a1(y,x,1),[]); c1=colorbar; colormap('gray');
        maxa=max(max(a1(y,x,1)));mina=min(min(a1(y,x,1)));
        maxa=floor(maxa*100)/100; mina=ceil(mina*100)/100; l=maxa-mina;
        xlabel('Blurred'); title('a','Position',[-0.7,0.95]);
        set(c1,'YTick',mina:l:maxa);
        set(c1,'FontName','Helvetica','FontSize',16);
        set(gca,'YTick',mina:l:maxa);
        set(gca,'FontName','Helvetica','FontSize',16);
        
        c2=subplot(2,4,2); hold on;
        maxa=max(max(a));mina=min(min(a));
        maxa=floor(maxa*100)/100; mina=ceil(mina*100)/100; l=maxa-mina;
        
        xlim([0,max(xx)*1.05]); ylim([mina*0.99,1.01*maxa]);
        xlabel('Position \mum'); ylabel('Normalized intensity'); title('b','Position',[-0.8,1.038]);
        set(c2,'XTick',0:1:floor(max(xx)),'YTick',mina:l:maxa);
        set(c2,'FontName','Helvetica','FontSize',16);
        plot(xx,a(1,:),'k--','LineWidth',line_width);
        plot(xx,a(2,:),'b','LineWidth',line_width); hold off;
        
        subplot(2,4,5); imshow(a2(y,x,1),[]); c3=colorbar; colormap('gray'); 
        maxa=max(max(a2(y,x,1)));mina=min(min(a2(y,x,1)));
        maxa=floor(maxa*100)/100; mina=ceil(mina*100)/100; l=maxa-mina;
        xlabel('Deblurred'); title('c','Position',[-0.7,0.95]);
        set(c3,'YTick',mina:l:maxa);
        set(c3,'FontName','Helvetica','FontSize',16);
        set(gca,'YTick',mina:l:maxa);
        set(gca,'FontName','Helvetica','FontSize',16);
        
        c4=subplot(2,4,6); hold on;
        maxa=max(max(a));mina=min(min(a));
        maxa=floor(maxa*100)/100; mina=ceil(mina*100)/100; l=maxa-mina;
        xlim([0,max(xx)*1.05]); ylim([mina*0.99,1.01*maxa]);
        xlabel('Position \mum'); ylabel('Normalized intensity'); title('d','Position',[-0.8,1.038]);
        set(c4,'XTick',0:1:floor(max(xx)));
        set(c4,'YTick',mina:l:maxa);
        set(c4,'FontName','Helvetica','FontSize',16);
        plot(xx,a(1,:),'k--','LineWidth',line_width);
        plot(xx,a(3,:),'r','LineWidth',line_width); hold off;
        %%%%% fig. 3/4(a-d)
    end

    ssims(1)=ssim(a1(y,x,1),a0(y,x));
    ssims(2)=ssim(a2(y,x,1),a0(y,x));
end

