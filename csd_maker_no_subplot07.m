function [cax] = csd_maker_no_subplot07(data,time,muavg,timerange,artrange,caxset,axpos,subpl)

factor=1;

% opengl ('neverselect')

try
ti(1)=max(find(time<=timerange(1)));
ti(2)=max(find(time<=timerange(2)));
catch
    ti(1) = 1;
    ti(2) = length(time);
end
ar(1)=max(find(time<=artrange(1)));
ar(2)=max(find(time<=artrange(2)));

if sum(abs(artrange))>0
    artremove=mean(data(:,max(find(time<=-25)):max(find(time<=-3))),2);
    for i=ar(1):ar(2)
        data(:,i)=artremove;
    end
end

num_of_points=size(data,2);
xi=1:1:num_of_points*factor;
x=1:num_of_points*factor/num_of_points:num_of_points*factor;

for i=1:size(data,1)
    yi(i,:) = interp1(x,squeeze(data(i,:)),xi);
end
clear data

csd=yi(:,ti(1):ti(2));

GRID_SCALE2         = 1000;
SHADING             = 'flat';
% opengl 'neverselect';
x                   = zeros(1,size(csd,1));
y                   = linspace(0.5,-0.5,size(csd,1));
Xmapsize            = 1; 
Ymapsize            = 1;
rmax                = 0.5;
xmin        = min(-.5*Xmapsize,min(x)); xmax = max(0.5*Xmapsize,max(x));
ymin        = min(-.5*Ymapsize,min(y)); ymax = max(0.5*Ymapsize,max(y));

[x,y]       = meshgrid(linspace(xmin,xmax,size(csd,2)),y);
[Xi,Yi]     = meshgrid(linspace(xmin,xmax,size(csd,2)),linspace(ymax,ymin,GRID_SCALE2));

Zi          = interp2(x,y,csd,Xi,Yi,'cubic');
delta = Xi(2)-Xi(1);

set(gcf,'CurrentAxes',subpl)
set(gca,'XColor','r','YColor','r','Visible','off');
set(gca,'Xlim',[-rmax*1*Xmapsize rmax*1*Xmapsize],'Ylim',[-rmax*1*Ymapsize  rmax*1*Ymapsize]);
hold on
surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);

ax = gca;

if muavg==0 || muavg==4
    colormap(subpl,hot)
    if max(max(csd))>0
        caxis([-max(max(csd))*0.2 max(max(csd))*0.8])
    else
        caxis([max(max(csd)) -max(max(csd))])
    end
        
elseif muavg==3
    colormap (subpl,hot)
    if abs(min(min(csd)))>max(max(csd))
        caxis([min(min(csd))*0.6 abs(min(min(csd)))*0.6])
    else
        caxis([-max(max(csd))*0.6 max(max(csd))*0.6])    
    end
elseif muavg==1

    colormap (subpl,jet)
    if abs(min(min(csd)))>max(max(csd))
        caxis([min(min(csd))*0.8 abs(min(min(csd)))*0.8])
    else
        caxis([-max(max(csd))*0.8 max(max(csd))*0.8])    
    end
    
 
else
    colormap (subpl,jet)
    try
    caxis([min(min(csd))*0.8 -min(min(csd))*0.2])
    catch
    end
end


% if isempty(caxset)
%     
% else
%     caxis(caxset)
% end

cax=caxis;

cmap=colormap;
if muavg==1
    cmap=flipud(cmap);
end
colormap(ax, cmap);


set(gca,'Position',axpos)
handles.mapax = axes( 'Position',axpos,...
        'Color','none',...
        'XColor','k','YColor','k',...
        'XLim', [time(ti(1)) time(ti(2))], 'YLim',[1 size(csd,1)],...
        'Ytick',[1:size(csd,1)],'YDir', 'reverse',...
        'Fontsize',16);
        H=gca;
        H.LineWidth=1.5; %axis thickness
hold on        
