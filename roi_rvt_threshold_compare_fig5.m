clc;clear;close all; warning off;

t=[0.06,0.07,0.1];
pages=1; r=1:25;
a0=readTIF('E:\data\2024.9.4\4-2\compares.tif');
for i=1:3
    pp0{i}=getParticles(a0(:,:,i),r,t(i),15,5,i);
end

[rois,p0]=getROIs(a0,pp0{1},15);
ref=getAveraged(rois);

for i=2:3
    fs1=[names0{i},'1.tif'];
    fs1=[fs0,fs1];
    a0=readsTIF(fs1,1);
    [rois,p1]=getROIs(a0,pp0{2},15);
    for j=1:size(rois,3)
        pcc(j,i-1)=getPearson(ref,rois(:,:,j));
        % subplot(1,2,1); imshow(ref);
        % subplot(1,2,2); imshow(rois(:,:,j));
        psnrs(j,i-1)=ssim(ref,rois(:,:,j));
    end
end

function [rois,ps]=getROIs(a,p,r)
    j=0;
    for i=1:size(p,1)
        x=max(1,p(i,3)-r):min(size(a,2),p(i,3)+r);
        y=max(1,p(i,2)-r):min(size(a,1),p(i,2)+r);
        b=a(y,x);
        if size(b,1)+size(b,2)==4*r+2
            j=j+1;
            rois(:,:,j)=normalMaxMin(b);
            ps(:,j)=[p(i,3),p(i,2)];
        end
    end
end


function p=fit2DGaussian(a,ps)
    x=1:size(a,2); y=1:size(a,1);
    a=a+eps;
    [X Y]=meshgrid(x,y);
    p0=[1,ps(1),ps(2),3,3];
    options = optimoptions('fmincon','MaxIter',500,'Display','off');
    f=@(p) sum(sum(abs(a-p(1)*exp(-0.5*(X-p(2)).*(X-p(2))./(p(4)*p(4))-0.5*(Y-p(3)).*(Y-p(3))./(p(5)*p(5))))));
    p=fmincon(f,p0,[],[],[],[],[0.00001,ps(1)-1,ps(2)-1,1,1],[1,ps(2)+1,ps(2)+1,6,6],[],options);
end
function pp=getParticles(a,r,t,l,l0,idx)
    b=RVT(a,r);
    b=normalMaxMin(b); b(b<t)=0;
    subplot(2,3,idx+3); imshow(b,[]); colormap("turbo"); hold on;
    c=b;
    maxc=max(max(c));
    i=0;
    r0=1;
    while(maxc>0)   
        i=1+i;  
        [p(i,1),p(i,2)]=find(c==maxc);
        % rectangle('Position',[p(i,2)-r0,p(i,1)-r0,2*r0,2*r0],'Curvature',[1 1],'FaceColor','r');
        x1=max(1,p(i,2)-l0):min(size(c,2),p(i,2)+l0);
        y1=max(1,p(i,1)-l0):min(size(c,1),p(i,1)+l0);
        rs=eps+1-a(y1,x1); rs(rs<0)=0;
        % 
        % % drawPoints(b,p(i,:),1);
        % 
        x=max(1,p(i,2)-l):min(size(c,2),p(i,2)+l);
        y=max(1,p(i,1)-l):min(size(c,1),p(i,1)+l);
        aa=1-a(y,x);
        rois{i}=zeros(length(y),length(x));
        rois{i}(min(y1)-min(y)+1:max(y1)-min(y)+1,min(x1)-min(x)+1:max(x1)-min(x)+1)=rs;
        % % subplot(4,4,i); imshow(aa,[]); %imshow(rois{i},[]);
        % 
        ps(:,i)=fit2DGaussian(rois{i},[p(i,2)-x(1),p(i,1)-y(1)]);
        pp(i,2)=ps(3,i)+y(1)-1; pp(i,3)=ps(2,i)+x(1)-1; pp(i,1)=ps(1,i);

        % rectangle('Position',[pp(i,2)-r0,pp(i,1)-r0,2*r0,2*r0],'Curvature',[1 1],'FaceColor','r');

        c(y,x)=0;
        maxc=max(max(c));
    end
    pp=round(pp);
    drawPoints(b,a,pp,1,idx);
end

function drawPoints(a,b,p,r,idx)
    p=round(p);
    % figure;  subplot(1,2,1); imshow(a,[]);
    subplot(2,3,idx); 
    % fig=figure;
    imshow(b,[]); colormap("gray");
    
    % hold on; 
    % for i=1:size(p,1)
    %     x=max(1,p(i,3)-r); y=max(1,p(i,2)-r);
    %     rectangle('Position',[x,y,2*r,2*r],'Curvature',[1 1],'FaceColor','r');
    % end
    % hold off;

    % frames=frame2im(getframe(fig));
    % fs1=sprintf('E:/data/2024.9.4/4-2/procesed%d.bmp',idx);
    % imwrite(frames,fs1);
end
function f1=getF1Score(p0,p,r0)
    m=size(p0,1); n=size(p,1);
    tp=0;
    for i=1:m
        for j=1:n
            r=hypot(p0(i,1)-p(j,1),p0(i,2)-p(j,2));
            if r<r0
                tp=tp+1;
            end
        end
    end
    fp=n-tp; fn=m-tp;
    f1=2*tp/(2*tp+fp+fn);
end
