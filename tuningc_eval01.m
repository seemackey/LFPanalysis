directory = 'C:\Users\peter\Dropbox\temp\tono\';
filenames = {};

if isempty(filenames)
    filenames=[];
    f=dir( [ directory 'tun01*.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
end 


%%%%% for a1, belt eval
tonogaps    = [0 1 2 3];
tunindex01  = [];
tunindex01c = [];
FFs         = [];
arej_ratio  = [];
locm        = [];
locc        = [];
locg        = [];
maxamp      = [];
bfavg       = [];
bfsm         = [];
bfsc        = [];
gran_auto   = [];
trialnos    = [];
sigresponses = [];
sigresponse_channels = [];
sigresponse_none = [];
a0 = 0;

bool_manual     = 0;  %%% 0 for auto, 1 for manual
gran_manual     = [];

for filecik=1:length(filenames)
    load([directory filenames{filecik}])
    
    %     %%% correction
    %     params_tc.seleegs.time = [4 24];
    
    
    arej_ratio(filecik)=length(params_tc.deleted_trials)/length(trig_new.trig1s{1});
    locm(filecik,:)=tunamp.max_loc.m;
    locc(filecik,:)=tunamp.max_loc.c_h;
    locg(filecik,:)=tunamp.max_loc.c_hgh;
    trialnos(filecik)=length(trig_new.trig1);
    
    if bool_manual == 1
        tunamp.max_loc.m(1) = gran_manual(filecik);
        tunamp_mean = tunamp.mean.m(tunamp.max_loc.m(1),1:end-1);
        [c,tunamp.max_loc.m(2)]=max(tunamp_mean);
    end
    
    tunamp_mean = tunamp.mean.m(tunamp.max_loc.m(1),:);
    tonoarray = 1:length(tunamp_mean)-1;
    bf = tunamp.max_loc.m(2);
    bfsm(filecik)=bf;
    for i1=1:length(tonogaps)
        z = find(tonoarray<bf-tonogaps(i1) | tonoarray>bf+tonogaps(i1));
        
        if mean(tunamp_mean(z))>0
            tunindex01(filecik,i1)=tunamp_mean(bf)/mean(tunamp_mean(z));
        else
            zz = find(tunamp_mean<0);
            tunamp_mean(zz)=0;
            tunindex01(filecik,i1)=tunamp_mean(bf)/mean(tunamp_mean(z));
        end
    end
    
    tunamp_mean = tunamp.mean.c_h(tunamp.max_loc.c_h(1),:);
    tonoarray = 1:length(tunamp_mean)-1;
    bf = tunamp.max_loc.c_h(2);
    bfsc(filecik)=bf;
    
    for i1=1:length(tonogaps)
        z = find(tonoarray<bf-tonogaps(i1) | tonoarray>bf+tonogaps(i1));
        
        if mean(tunamp_mean(z))>0
            tunindex01c(filecik,i1)=tunamp_mean(bf)/mean(tunamp_mean(z));
        else
            zz = find(tunamp_mean<0);
            tunamp_mean(zz)=0;
            tunindex01c(filecik,i1)=tunamp_mean(bf)/mean(tunamp_mean(z));
        end
    end
    
    FFs(filecik)=FanoFactor.m(tunamp.max_loc.m(1),tunamp.max_loc.m(2));
    
    avg = squeeze(avgs.m(tunamp.max_loc.m(1),tunamp.max_loc.m(2),:));
    maxamp(filecik) = max(avg(max(find(avgs.time<=params_tc.seleegs.time(1))):max(find(avgs.time<=params_tc.seleegs.time(2)))));
    bfavg(filecik,:)=avg;
    
    
    ttype1 = trig_new.ttype1;
    a=0; trigtype=[];
    for i=0:3000
        z=find(trig_new.ttype1==i);
        if ~isempty(z)
            a=a+1;
            trigtype(a)=i;
        end
    end

    t_baseline=max(find(avgs.time<=params_tc.baseline(1))):max(find(avgs.time<=params_tc.baseline(2)));
    t_sel = max(find(avgs.time<=params_tc.seleegs.time(1))):max(find(avgs.time<=params_tc.seleegs.time(2)));
%     amps_baseline = squeeze(mean(eegs.eegm(:,:,t_baseline),3));
%     amps_sel = squeeze(mean(eegs.eegm(:,:,t_sel),3));
%     sigresponse_m = [];
%     for trigcik=1:length(trigtype)
%         z = find(trig_new.ttype1==trigtype(trigcik));
%         trialnumbers(trigcik)=length(z);
%         for chancik = 1:size(eegs.eegm,1)
%             sigresponse_m(chancik,trigcik) = signrank(amps_baseline(chancik,z),amps_sel(chancik,z));
%         end
%     end
%     sigresponses(filecik) = sigresponse_m(tunamp.max_loc.m(1),tunamp.max_loc.m(2));
%     a = 0;
%     chs = [];
%     for chancik = 1:size(eegs.eegm,1)
%         if min(sigresponse_m(chancik,:))<0.05/size(eegs.eegm,1)/length(trigtype)
%             a=a+1;
%             chs(a)=chancik;
%         end
%     end
%     sigresponse_channels{filecik} = chs;
%     if isempty(chs)
%         a0=a0+1;
%     sigresponse_none(a0)=filecik;
%     end
end


gran_auto = locm(:,1);


%%%% belt
z_belt = find(tunindex01(:,3)<10 | maxamp'<2);
filenames_belt=filenames(z_belt);
% zz = z;
% for i1=1:length(z)
%     if bfs(z(i1))<8
%         if  tunindex01(z(i1),3)<5 | maxamp(i1)'<2
%         else
%             z(i1)=0;
%         end
%     end
% end
% z(z==0)=[];
 fprintf(['ratio of belt ' num2str(length(z_belt)/length(tunindex01)) '\n'])

%%%% a1
z_a1 = find(tunindex01(:,3)>10 & maxamp'>2);
filenames_a1=filenames(z_a1);
% zz = z;
% for i1=1:length(z)
%     if bfs(z(i1))<8
%         if  tunindex01(z(i1),3)>5 & maxamp(z)'>2
%         else
%             z(i1)=0;
%         end
%     end
% end
% z(z==0)=[];
 fprintf(['ratio of a1 ' num2str(length(z_a1)/length(tunindex01)) '\n'])

% z = find(tunindex01(:,2)<7);
% length(z)/length(tunindex01)

 clear a a0 amps amps_baseline amps_sel avg avgs bf bool bool_manual chancik chs eegs f FanoFactor ff filecik i i1 
 clear locc locg locm params_tc sigresponse_m sigresponses t_baseline t_sel tonoarray tonogaps trialnumbers trig_new trigcik trigtype ttype1 tunamp tunamp_mean z zz
 