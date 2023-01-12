% windows01

function [seg, segvec] = windows02(y,sizelimit,edge1, siglength)

%%%% y          - data points (vector) of something found significant or e.g. above or below some limit
%%%% sizelimit  - minimum number of data points in a segment
%%%% edge1      - edges to reject around non-segment periods in number of data points
%%%% siglength  - length of original signal in number of data points



seg     = [];
segvec  = [];

if length(y)>1

    
    y1=diff(y);
   
    
    a=1; b=a; c=1; seg = [];
    while a<length(y1)+1
        while y1(a)<2 && a<length(y1)
            a=a+1;
        end
        seg(c,1)=y(b); seg(c,2)=y(a);
        a=a+1; b = a; c = c+1;
    end
    

    
    if size(seg,1)>1
        seg(c-1,2)=seg(c-1,2)+1;
    end
    
    
    x =diff(seg')+1;
    
    z=find(x<sizelimit);
    seg(z,:)=[];
    

    
    for i1=1:size(seg,1)
        seg(i1,1)=seg(i1,1)-edge1;
        seg(i1,2)=seg(i1,2)+edge1;
    end
    
%     z=find(seg(:,1)<1);
%     seg(z,:)=[];
%     z=find(seg(:,2)>siglength);
%     seg(z,:)=[];
    z=find(seg(:,1)<1);
    seg(z,1)=1;
    z=find(seg(:,2)>siglength);
    seg(z,2)=siglength;


    for ii=1:size(seg,1)
        segvec=[segvec seg(ii,1):seg(ii,2)];
%         try
%             segvec=[segvec seg(ii,1):seg(ii,2)];
%         catch
%             segvec=[segvec; seg(ii,1):seg(ii,2)];
%             % segvec=[segvec];
%         end
    end
    
end

y = unique(segvec);
    
seg     = [];
segvec  = [];

if length(y)>1
    
    y1=diff(y);
    a=1; b=a; c=1; seg = [];
    while a<length(y1)+1
        while y1(a)<2 && a<length(y1)
            a=a+1;
        end
        seg(c,1)=y(b); seg(c,2)=y(a);
        a=a+1; b = a; c = c+1;
    end
    if size(seg,1)>=1
        seg(c-1,2)=seg(c-1,2)+1;
    end
    
    
    x =diff(seg')+1;
    z=find(x<sizelimit);
    seg(z,:)=[];
    
    for i1=1:size(seg,1)
        seg(i1,1)=seg(i1,1)-edge1;
        seg(i1,2)=seg(i1,2)+edge1;
    end
    
    %     z=find(seg(:,1)<1);
%     seg(z,:)=[];
%     z=find(seg(:,2)>siglength);
%     seg(z,:)=[];
    z=find(seg(:,1)<1);
    seg(z,1)=1;
    z=find(seg(:,2)>siglength);
    seg(z,2)=siglength;
    
    for ii=1:size(seg,1)
        
        segvec=[segvec seg(ii,1):seg(ii,2)];
        
%         try
%             segvec=[segvec seg(ii,1):seg(ii,2)];
%         catch
%             segvec=[segvec; seg(ii,1):seg(ii,2)];
%             % segvec=[segvec];
%         end
    end
    
end
