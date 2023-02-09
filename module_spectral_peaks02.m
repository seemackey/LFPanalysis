
function [peaks,surfcohere,ccoef,params_speaks,posall] = module_spectral_peaks02(wraw, wlfp)

bool.coherence  = 0;

waitb           = waitbar(0);


winsize         = 1000;
winstep         = 1000;
% corrtype        = 'Pearson';
corrtype        = 'Spearman';

bool_ampcollect = 1;
bool_arej       = 0;
arej_limit      = 10;

peaks.c.frq = [];
peaks.c.amp = [];
peaks.c.diff = [];
peaks.c.dist = [];
peaks.c.ampmean = [];
peaks.c.chans = [];
peaks.c.spectra = [];

peaks.avgc.itc.frq      = [];
peaks.avgc.itc.time     = [];
peaks.avgc.itc.r        = [];
peaks.avgc.itc.p        = [];
peaks.avgc.itc.chans    = [];
peaks.avgc.itc.trigcik1 = [];
peaks.avgc.itc.strigtype   = [];

peaks.avgc.po.frq       = [];
peaks.avgc.po.time      = [];
peaks.avgc.po.amp       = [];
peaks.avgc.po.chans     = [];
peaks.avgc.po.trigcik1  = [];
peaks.avgc.po.strigtype  = [];

peaks.avgm.itc.frq      = [];
peaks.avgm.itc.time     = [];
peaks.avgm.itc.r        = [];
peaks.avgm.itc.p        = [];
peaks.avgm.itc.chans    = [];
peaks.avgm.itc.trigcik1 = [];
peaks.avgm.itc.strigtype   = [];

peaks.avgm.po.frq       = [];
peaks.avgm.po.time      = [];
peaks.avgm.po.amp       = [];
peaks.avgm.po.chans     = [];
peaks.avgm.po.trigcik1  = [];
peaks.avgm.po.strigtype  = [];

peaks.avge.itc.frq      = [];
peaks.avge.itc.time     = [];
peaks.avge.itc.r        = [];
peaks.avge.itc.p        = [];
peaks.avge.itc.chans    = [];
peaks.avge.itc.trigcik1 = [];
peaks.avge.itc.strigtype   = [];

peaks.avge.po.frq       = [];
peaks.avge.po.time      = [];
peaks.avge.po.amp       = [];
peaks.avge.po.chans     = [];
peaks.avge.po.trigcik1  = [];
peaks.avge.po.strigtype  = [];

peaks.avgb.itc.frq      = [];
peaks.avgb.itc.time     = [];
peaks.avgb.itc.r        = [];
peaks.avgb.itc.p        = [];
peaks.avgb.itc.chans    = [];
peaks.avgb.itc.trigcik1 = [];
peaks.avgb.itc.strigtype   = [];

peaks.avgb.po.frq       = [];
peaks.avgb.po.time      = [];
peaks.avgb.po.amp       = [];
peaks.avgb.po.chans     = [];
peaks.avgb.po.trigcik1  = [];
peaks.avgb.po.strigtype  = [];

peaks.avgh.itc.frq      = [];
peaks.avgh.itc.time     = [];
peaks.avgh.itc.r        = [];
peaks.avgh.itc.p        = [];
peaks.avgh.itc.chans    = [];
peaks.avgh.itc.trigcik1 = [];
peaks.avgh.itc.strigtype   = [];

peaks.avgh.po.frq       = [];
peaks.avgh.po.time      = [];
peaks.avgh.po.amp       = [];
peaks.avgh.po.chans     = [];
peaks.avgh.po.trigcik1  = [];
peaks.avgh.po.strigtype  = [];

peaks.avgl.itc.frq      = [];
peaks.avgl.itc.time     = [];
peaks.avgl.itc.r        = [];
peaks.avgl.itc.p        = [];
peaks.avgl.itc.chans    = [];
peaks.avgl.itc.trigcik1 = [];
peaks.avgl.itc.strigtype   = [];

peaks.avgl.po.frq       = [];
peaks.avgl.po.time      = [];
peaks.avgl.po.amp       = [];
peaks.avgl.po.chans     = [];
peaks.avgl.po.trigcik1  = [];
peaks.avgl.po.strigtype  = [];

peaks.m.frq = [];
peaks.m.amp = [];
peaks.m.diff = [];
peaks.m.dist = [];
peaks.m.ampmean = [];
peaks.m.chans = [];
peaks.m.spectra = [];

peaks.b.frq = [];
peaks.b.amp = [];
peaks.b.diff = [];
peaks.b.dist = [];
peaks.b.ampmean = [];
peaks.b.chans = [];
peaks.b.spectra = [];

peaks.e.frq = [];
peaks.e.amp = [];
peaks.e.diff = [];
peaks.e.dist = [];
peaks.e.ampmean = [];
peaks.e.chans = [];
peaks.e.spectra = [];

peaks.h.frq = [];
peaks.h.amp = [];
peaks.h.diff = [];
peaks.h.dist = [];
peaks.h.ampmean = [];
peaks.h.chans = [];
peaks.h.spectra = [];

peaks.l.frq = [];
peaks.l.amp = [];
peaks.l.diff = [];
peaks.l.dist = [];
peaks.l.ampmean = [];

peaks.l.chans = [];
peaks.l.spectra = [];

a_cntc=0;
a_cntm=0;
a_cntb=0;
a_cnth=0;
a_cnte=0;
a_cntl=0;

a_avgc=0;
a_avgm=0;
a_avgb=0;
a_avge=0;
a_avgl=0;
a_avgh=0;

posall.m        = [];
posall.c        = [];


%%%%%%%% sliding variables
stimbegin = 1;
stimend = size(wraw.cntc_po,3);

wsize=winsize/1000*wraw.adrate;
wstep=winstep/1000*wraw.adrate;

time_slide      =((stimbegin:wstep:(stimend-wsize))-stimbegin)/wraw.adrate;

%%%%%% compensate for new, 22 channel wavelet
% if size(wraw.cntc_po,1)>21
%     wraw.cntc_po=wraw.cntc_po(1:21,:,:);
%     wraw.cntc_ph=wraw.cntc_ph(1:21,:,:);
% end
% if size(wraw.cntm_po,1)>21
%     wraw.cntm_po=wraw.cntc_po(1:21,:,:);
%     wraw.cntm_ph=wraw.cntc_ph(1:21,:,:);
% end
% if size(wraw.cntc,1)>21
%     wraw.cntc=wraw.cntc(1:21,:);
%     wraw.cntc=wraw.cntc(1:21,:);
% end
% if size(wraw.cntm,1)>21
%     wraw.cntm=wraw.cntm(1:21,:);
%     wraw.cntm=wraw.cntm(1:21,:);
% end

posall.lfp      = [];
posall.m        = [];
posall.c        = [];
posall.b        = [];
posall.e        = [];
posall.slmua    = [];
posall.time_slide = time_slide;

poe = [];
poc = [];
pob = [];

pom = [];
poh = [];
slmua = [];


% sliding power for csd
a2=['moving csd'];
waitbar(0, waitb,a2);
if ~isempty(wraw.cntc_po)
    poc = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        poc(:,:,a)=squeeze(mean(wraw.cntc_po(:,:,i1:i1+wsize-1),3));
    end
    posall.c        = poc;
else
    poc = [];
end

% sliding power for bipolar
a2=['moving bipolar'];
waitbar(0, waitb,a2);
if ~isempty(wraw.cntb_po)
    pob = zeros(size(wraw.cntb_po,1),size(wraw.cntb_po,2),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        pob(:,:,a)=squeeze(mean(wraw.cntb_po(:,:,i1:i1+wsize-1),3));
    end
    
else pob =[];
    posall.b        = pob;
end

% sliding power for field
a2=['moving field'];
waitbar(0, waitb,a2);
if ~isempty(wraw.cnte_po)
    poe = zeros(size(wraw.cnte_po,1),size(wraw.cnte_po,2),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        poe(:,:,a)=squeeze(mean(wraw.cnte_po(:,:,i1:i1+wsize-1),3));
    end
    posall.e        = poe;
    
else 
    poe =[];
end

% sliding power for mua
a2=['moving mua'];
waitbar(0, waitb,a2);
if ~isempty(wraw.cntm_po)
    pom = zeros(size(wraw.cntm_po,1),size(wraw.cntm_po,2),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        pom(:,:,a)=squeeze(mean(wraw.cntm_po(:,:,i1:i1+wsize-1),3));
    end
    posall.m        = pom;
    % sliding mua
    slmua = zeros(size(wraw.cntm,1),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        slmua(:,a)=squeeze(mean(wraw.cntm(:,i1:i1+wsize-1),2));
    end
    posall.slmua        = slmua;
else
    slmua =[];
end


% sliding power for hgamma
a2=['moving hgamma'];
waitbar(0, waitb,a2);
if ~isempty(wraw.cnth_po)
    poh = zeros(size(wraw.cnth_po,1),size(wraw.cnth_po,2),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        poh(:,:,a)=squeeze(mean(wraw.cnth_po(:,:,i1:i1+wsize-1),3));
    end
    posall.h        = poh;
    % sliding hgamma
    slhgamma = zeros(size(wraw.cnth,1),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        slhgamma(:,a)=squeeze(mean(wraw.cnth(:,i1:i1+wsize-1),2));
    end
    posall.slhgamma        = slhgamma;
else
    slhgamma =[];
end

a2=['moving surface'];
waitbar(0, waitb,a2);
if ~isempty(wlfp.cnt_po)
    % sliding power for lfp
    pol = zeros(size(wlfp.cnt_po,1),size(wlfp.cnt_po,2),length(time_slide));
    a = 0;
    for i1=stimbegin:wstep:(stimend-wsize)
        a = a+1;
        pol(:,:,a)=squeeze(mean(wlfp.cnt_po(:,:,i1:i1+wsize-1),3));
    end
    posall.lfp        = pol;
else 
    pol =[];
end


%%%% finding peaks in the whole averaged power

a2=['peakfind'];
waitbar(0, waitb,a2);
%%%%% for csd
if ~isempty(wraw.cntc_po)
    for chancik = 1:size(wraw.cntc_po,1)
        sp1 = squeeze(mean(wraw.cntc_po(chancik,:,:),3));
        v=diff(sp1);
        zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zci_pos = zci(v(zci)>0)+1;
        zci_neg = zci(v(zci)<0)+1;
        peaks.c.frq = [peaks.c.frq; zci_pos];
        peaks.c.amp = [peaks.c.amp; sp1(zci_pos)'];
        peaks.c.chans = [peaks.c.chans; zeros(length(zci_pos),1)+chancik];
        
        peaks.c.ampmean = [peaks.c.ampmean; zeros(length(zci_pos),1)+mean(sp1)];
        
        z = find(v(zci)>0);
        for i1 = 1:length(zci_pos)
            a_cntc=a_cntc+1;
            
            peaks.c.spectra(a_cntc,:)=sp1;
            
            try
                d(1) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)-1)+1);
                dist(1) = zci(z(i1))-zci(z(i1)-1);
            catch
                d(1) = -1;
                dist(1) = -1;
            end
            try
                d(2) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)+1)+1);
                dist(2) = zci(z(i1)+1)-zci(z(i1));
            catch
                d(2) = -1;
                dist(2) = -1;
            end
            peaks.c.diff(a_cntc,:) = d;
            peaks.c.dist(a_cntc,:) = dist;
        end
    end
end

%%%%% for bipolar
if ~isempty(wraw.cntb_po)
    for chancik = 1:size(wraw.cntb_po,1)
        sp1 = squeeze(mean(wraw.cntb_po(chancik,:,:),3));
        v=diff(sp1);
        zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zci_pos = zci(v(zci)>0)+1;
        zci_neg = zci(v(zci)<0)+1;
        peaks.b.frq = [peaks.b.frq; zci_pos];
        peaks.b.amp = [peaks.b.amp; sp1(zci_pos)'];
        peaks.b.chans = [peaks.b.chans; zeros(length(zci_pos),1)+chancik];
        
        peaks.b.ampmean = [peaks.b.ampmean; zeros(length(zci_pos),1)+mean(sp1)];
        
        z = find(v(zci)>0);
        for i1 = 1:length(zci_pos)
            a_cntb=a_cntb+1;
            
            peaks.b.spectra(a_cntb,:)=sp1;
            
            try
                d(1) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)-1)+1);
                dist(1) = zci(z(i1))-zci(z(i1)-1);
            catch
                d(1) = -1;
                dist(1) = -1;
            end
            try
                d(2) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)+1)+1);
                dist(2) = zci(z(i1)+1)-zci(z(i1));
            catch
                d(2) = -1;
                dist(2) = -1;
            end
            peaks.b.diff(a_cntb,:) = d;
            peaks.b.dist(a_cntb,:) = dist;
        end
    end
end

%%%%% for field
if ~isempty(wraw.cnte_po)
    for chancik = 1:size(wraw.cnte_po,1)
        sp1 = squeeze(mean(wraw.cnte_po(chancik,:,:),3));
        v=diff(sp1);
        zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zci_pos = zci(v(zci)>0)+1;
        zci_neg = zci(v(zci)<0)+1;
        peaks.e.frq = [peaks.e.frq; zci_pos];
        peaks.e.amp = [peaks.e.amp; sp1(zci_pos)'];
        peaks.e.chans = [peaks.e.chans; zeros(length(zci_pos),1)+chancik];
        
        peaks.e.ampmean = [peaks.e.ampmean; zeros(length(zci_pos),1)+mean(sp1)];
        
        z = find(v(zci)>0);
        for i1 = 1:length(zci_pos)
            a_cnte=a_cnte+1;
            
            peaks.e.spectra(a_cnte,:)=sp1;
            
            try
                d(1) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)-1)+1);
                dist(1) = zci(z(i1))-zci(z(i1)-1);
            catch
                d(1) = -1;
                dist(1) = -1;
            end
            try
                d(2) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)+1)+1);
                dist(2) = zci(z(i1)+1)-zci(z(i1));
            catch
                d(2) = -1;
                dist(2) = -1;
            end
            peaks.e.diff(a_cnte,:) = d;
            peaks.e.dist(a_cnte,:) = dist;
        end
    end
end

%%%%% for mua
if ~isempty(wraw.cntm_po)
    for chancik = 1:size(wraw.cntm_po,1)
        sp1 = squeeze(mean(wraw.cntm_po(chancik,:,:),3));
        v=diff(sp1);
        zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zci_pos = zci(v(zci)>0)+1;
        zci_neg = zci(v(zci)<0)+1;
        peaks.m.frq = [peaks.m.frq; zci_pos];
        peaks.m.amp = [peaks.m.amp; sp1(zci_pos)'];
        peaks.m.chans = [peaks.m.chans; zeros(length(zci_pos),1)+chancik];
        
        peaks.m.ampmean = [peaks.m.ampmean; zeros(length(zci_pos),1)+mean(sp1)];
        
        z = find(v(zci)>0);
        for i1 = 1:length(zci_pos)
            a_cntm=a_cntm+1;
            
            peaks.m.spectra(a_cntm,:)=sp1;
            
            try
                d(1) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)-1)+1);
                dist(1) = zci(z(i1))-zci(z(i1)-1);
            catch
                d(1) = -1;
                dist(1) = -1;
            end
            try
                d(2) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)+1)+1);
                dist(2) = zci(z(i1)+1)-zci(z(i1));
            catch
                d(2) = -1;
                dist(2) = -1;
            end
            peaks.m.diff(a_cntm,:) = d;
            peaks.m.dist(a_cntm,:) = dist;
        end
    end
end


%%%%% for hgamma
if ~isempty(wraw.cnth_po)
    for chancik = 1:size(wraw.cnth_po,1)
        sp1 = squeeze(mean(wraw.cnth_po(chancik,:,:),3));
        v=diff(sp1);
        zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zci_pos = zci(v(zci)>0)+1;
        zci_neg = zci(v(zci)<0)+1;
        peaks.h.frq = [peaks.h.frq; zci_pos];
        peaks.h.amp = [peaks.h.amp; sp1(zci_pos)'];
        peaks.h.chans = [peaks.h.chans; zeros(length(zci_pos),1)+chancik];
        
        peaks.h.ampmean = [peaks.h.ampmean; zeros(length(zci_pos),1)+mean(sp1)];
        
        z = find(v(zci)>0);
        for i1 = 1:length(zci_pos)
            a_cnth=a_cnth+1;
            
            peaks.h.spectra(a_cnth,:)=sp1;
            
            try
                d(1) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)-1)+1);
                dist(1) = zci(z(i1))-zci(z(i1)-1);
            catch
                d(1) = -1;
                dist(1) = -1;
            end
            try
                d(2) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)+1)+1);
                dist(2) = zci(z(i1)+1)-zci(z(i1));
            catch
                d(2) = -1;
                dist(2) = -1;
            end
            peaks.h.diff(a_cnth,:) = d;
            peaks.h.dist(a_cnth,:) = dist;
        end
    end
end

%%%%% for lfp
if ~isempty(wlfp.cnt_po)
    for chancik = 1:size(wlfp.cnt_po,1)
        sp1 = squeeze(mean(wlfp.cnt_po(chancik,:,:),3));
        v=diff(sp1);
        zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zci_pos = zci(v(zci)>0)+1;
        zci_neg = zci(v(zci)<0)+1;
        peaks.l.frq = [peaks.l.frq; zci_pos];
        peaks.l.amp = [peaks.l.amp; sp1(zci_pos)'];
        peaks.l.chans = [peaks.l.chans; zeros(length(zci_pos),1)+chancik];
        
        peaks.l.ampmean = [peaks.l.ampmean; zeros(length(zci_pos),1)+mean(sp1)];
        
        z = find(v(zci)>0);
        for i1 = 1:length(zci_pos)
            a_cntl=a_cntl+1;
            
            peaks.l.spectra(a_cntl,:)=sp1;
            
            try
                d(1) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)-1)+1);
                dist(1) = zci(z(i1))-zci(z(i1)-1);
            catch
                d(1) = -1;
                dist(1) = -1;
            end
            try
                d(2) = sp1(zci(z(i1))+1)-sp1(zci(z(i1)+1)+1);
                dist(2) = zci(z(i1)+1)-zci(z(i1));
            catch
                d(2) = -1;
                dist(2) = -1;
            end
            peaks.l.diff(a_cntl,:) = d;
            peaks.l.dist(a_cntl,:) = dist;
        end
    end
end

%%%%%%%%% peaksearch on averages

%%% csd
if ~isempty(wraw.avg.poc)
    for trigcik1 = 1:length(wraw.avg.poc)
        for trigcik2 = 1:size(wraw.avg.poc{trigcik1},1)
            for chancik = 1:size(wraw.avg.poc{trigcik1},2)
                xmap = squeeze(wraw.avg.itcc{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap,1);
                    v=diff(xmap(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap,2);
                    v=diff(xmap(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc = [];
                peakcounter = 0;
                
                for frqcik = 1:size(xmap,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                r = zeros(1,size(peakloc,1));
                p = zeros(1,size(peakloc,1));
                for i1=1:size(peakloc,1)
                    r(i1) = xmap(peakloc(i1,1),peakloc(i1,2));
                    p(i1)  = squeeze(wraw.avg.itcpc{1}(trigcik2,chancik,peakloc(i1,1),peakloc(i1,2)));
                end
                
                xmap_po = squeeze(wraw.avg.poc{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap_po,1);
                    v=diff(xmap_po(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap_po,2);
                    v=diff(xmap_po(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc_po = [];
                peakcounter = 0;
                for frqcik = 1:size(xmap_po,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc_po(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                amp = zeros(1,size(peakloc_po,1));
                for i1=1:size(peakloc_po,1)
                    amp(i1) = xmap_po(peakloc_po(i1,1),peakloc_po(i1,2));
                end
                
                if ~isempty(peakloc)
                    peaks.avgc.itc.frq      = [peaks.avgc.itc.frq;  peakloc(:,1)];
                    peaks.avgc.itc.time     = [peaks.avgc.itc.time;  peakloc(:,2)];
                    peaks.avgc.itc.r        = [peaks.avgc.itc.r;  r'];
                    peaks.avgc.itc.p        = [peaks.avgc.itc.p;  p'];
                    peaks.avgc.itc.chans    = [peaks.avgc.itc.chans; zeros(size(peakloc,1),1)+chancik];
                    peaks.avgc.itc.trigcik1 = [peaks.avgc.itc.trigcik1; zeros(size(peakloc,1),1)+trigcik1];
                    peaks.avgc.itc.strigtype = [peaks.avgc.itc.strigtype; zeros(size(peakloc,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
                
                if ~isempty(peakloc_po)
                    peaks.avgc.po.frq      = [peaks.avgc.po.frq;  peakloc_po(:,1)];
                    peaks.avgc.po.time     = [peaks.avgc.po.time;  peakloc_po(:,2)];
                    peaks.avgc.po.amp       = [peaks.avgc.po.amp; amp'];
                    peaks.avgc.po.chans    = [peaks.avgc.po.chans; zeros(size(peakloc_po,1),1)+chancik];
                    peaks.avgc.po.trigcik1 = [peaks.avgc.po.trigcik1; zeros(size(peakloc_po,1),1)+trigcik1];
                    peaks.avgc.po.strigtype = [peaks.avgc.po.strigtype; zeros(size(peakloc_po,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
            end
        end
    end
end

%%% mua
if ~isempty(wraw.avg.pom)
    for trigcik1 = 1:length(wraw.avg.pom)
        for trigcik2 = 1:size(wraw.avg.pom{trigcik1},1)
            for chancik = 1:size(wraw.avg.pom{trigcik1},2)
                xmap = squeeze(wraw.avg.itcm{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap,1);
                    v=diff(xmap(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap,2);
                    v=diff(xmap(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc = [];
                peakcounter = 0;
                
                for frqcik = 1:size(xmap,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                r = zeros(1,size(peakloc,1));
                p = zeros(1,size(peakloc,1));
                for i1=1:size(peakloc,1)
                    r(i1) = xmap(peakloc(i1,1),peakloc(i1,2));
                    p(i1)  = squeeze(wraw.avg.itcpm{1}(trigcik2,chancik,peakloc(i1,1),peakloc(i1,2)));
                end
                
                xmap_po = squeeze(wraw.avg.pom{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap_po,1);
                    v=diff(xmap_po(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap_po,2);
                    v=diff(xmap_po(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc_po = [];
                peakcounter = 0;
                for frqcik = 1:size(xmap_po,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc_po(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                amp = zeros(1,size(peakloc_po,1));
                for i1=1:size(peakloc_po,1)
                    amp(i1) = xmap_po(peakloc_po(i1,1),peakloc_po(i1,2));
                end
                
                if ~isempty(peakloc)
                    peaks.avgm.itc.frq      = [peaks.avgm.itc.frq;  peakloc(:,1)];
                    peaks.avgm.itc.time     = [peaks.avgm.itc.time;  peakloc(:,2)];
                    peaks.avgm.itc.r        = [peaks.avgm.itc.r;  r'];
                    peaks.avgm.itc.p        = [peaks.avgm.itc.p;  p'];
                    peaks.avgm.itc.chans    = [peaks.avgm.itc.chans; zeros(size(peakloc,1),1)+chancik];
                    peaks.avgm.itc.trigcik1 = [peaks.avgm.itc.trigcik1; zeros(size(peakloc,1),1)+trigcik1];
                    peaks.avgm.itc.strigtype = [peaks.avgm.itc.strigtype; zeros(size(peakloc,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
                
                if ~isempty(peakloc_po)
                    peaks.avgm.po.frq      = [peaks.avgm.po.frq;  peakloc_po(:,1)];
                    peaks.avgm.po.time     = [peaks.avgm.po.time;  peakloc_po(:,2)];
                    peaks.avgm.po.amp       = [peaks.avgm.po.amp; amp'];
                    peaks.avgm.po.chans    = [peaks.avgm.po.chans; zeros(size(peakloc_po,1),1)+chancik];
                    peaks.avgm.po.trigcik1 = [peaks.avgm.po.trigcik1; zeros(size(peakloc_po,1),1)+trigcik1];
                    peaks.avgm.po.strigtype = [peaks.avgm.po.strigtype; zeros(size(peakloc_po,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
            end
        end
    end
end

%%% lfp
if ~isempty(wraw.avg.poe)
    for trigcik1 = 1:length(wraw.avg.poe)
        for trigcik2 = 1:size(wraw.avg.poe{trigcik1},1)
            for chancik = 1:size(wraw.avg.poe{trigcik1},2)
                xmap = squeeze(wraw.avg.itce{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap,1);
                    v=diff(xmap(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap,2);
                    v=diff(xmap(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc = [];
                peakcounter = 0;
                
                for frqcik = 1:size(xmap,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                r = zeros(1,size(peakloc,1));
                p = zeros(1,size(peakloc,1));
                for i1=1:size(peakloc,1)
                    r(i1) = xmap(peakloc(i1,1),peakloc(i1,2));
                    p(i1)  = squeeze(wraw.avg.itcpe{1}(trigcik2,chancik,peakloc(i1,1),peakloc(i1,2)));
                end
                
                xmap_po = squeeze(wraw.avg.poe{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap_po,1);
                    v=diff(xmap_po(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap_po,2);
                    v=diff(xmap_po(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc_po = [];
                peakcounter = 0;
                for frqcik = 1:size(xmap_po,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc_po(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                amp = zeros(1,size(peakloc_po,1));
                for i1=1:size(peakloc_po,1)
                    amp(i1) = xmap_po(peakloc_po(i1,1),peakloc_po(i1,2));
                end
                
                if ~isempty(peakloc)
                    peaks.avge.itc.frq      = [peaks.avge.itc.frq;  peakloc(:,1)];
                    peaks.avge.itc.time     = [peaks.avge.itc.time;  peakloc(:,2)];
                    peaks.avge.itc.r        = [peaks.avge.itc.r;  r'];
                    peaks.avge.itc.p        = [peaks.avge.itc.p;  p'];
                    peaks.avge.itc.chans    = [peaks.avge.itc.chans; zeros(size(peakloc,1),1)+chancik];
                    peaks.avge.itc.trigcik1 = [peaks.avge.itc.trigcik1; zeros(size(peakloc,1),1)+trigcik1];
                    peaks.avge.itc.strigtype = [peaks.avge.itc.strigtype; zeros(size(peakloc,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
                
                if ~isempty(peakloc_po)
                    peaks.avge.po.frq      = [peaks.avge.po.frq;  peakloc_po(:,1)];
                    peaks.avge.po.time     = [peaks.avge.po.time;  peakloc_po(:,2)];
                    peaks.avge.po.amp       = [peaks.avge.po.amp; amp'];
                    peaks.avge.po.chans    = [peaks.avge.po.chans; zeros(size(peakloc_po,1),1)+chancik];
                    peaks.avge.po.trigcik1 = [peaks.avge.po.trigcik1; zeros(size(peakloc_po,1),1)+trigcik1];
                    peaks.avge.po.strigtype = [peaks.avge.po.strigtype; zeros(size(peakloc_po,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
            end
        end
    end
end

%%% bipolar
if ~isempty(wraw.avg.pob)
    for trigcik1 = 1:length(wraw.avg.pob)
        for trigcik2 = 1:size(wraw.avg.pob{trigcik1},1)
            for chancik = 1:size(wraw.avg.pob{trigcik1},2)
                xmap = squeeze(wraw.avg.itcb{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap,1);
                    v=diff(xmap(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap,2);
                    v=diff(xmap(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc = [];
                peakcounter = 0;
                
                for frqcik = 1:size(xmap,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                r = zeros(1,size(peakloc,1));
                p = zeros(1,size(peakloc,1));
                for i1=1:size(peakloc,1)
                    r(i1) = xmap(peakloc(i1,1),peakloc(i1,2));
                    p(i1)  = squeeze(wraw.avg.itcpb{1}(trigcik2,chancik,peakloc(i1,1),peakloc(i1,2)));
                end
                
                xmap_po = squeeze(wraw.avg.pob{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap_po,1);
                    v=diff(xmap_po(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap_po,2);
                    v=diff(xmap_po(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc_po = [];
                peakcounter = 0;
                for frqcik = 1:size(xmap_po,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc_po(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                amp = zeros(1,size(peakloc_po,1));
                for i1=1:size(peakloc_po,1)
                    amp(i1) = xmap_po(peakloc_po(i1,1),peakloc_po(i1,2));
                end
                
                if ~isempty(peakloc)
                    peaks.avgb.itc.frq      = [peaks.avgb.itc.frq;  peakloc(:,1)];
                    peaks.avgb.itc.time     = [peaks.avgb.itc.time;  peakloc(:,2)];
                    peaks.avgb.itc.r        = [peaks.avgb.itc.r;  r'];
                    peaks.avgb.itc.p        = [peaks.avgb.itc.p;  p'];
                    peaks.avgb.itc.chans    = [peaks.avgb.itc.chans; zeros(size(peakloc,1),1)+chancik];
                    peaks.avgb.itc.trigcik1 = [peaks.avgb.itc.trigcik1; zeros(size(peakloc,1),1)+trigcik1];
                    peaks.avgb.itc.strigtype = [peaks.avgb.itc.strigtype; zeros(size(peakloc,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
                
                if ~isempty(peakloc_po)
                    peaks.avgb.po.frq      = [peaks.avgb.po.frq;  peakloc_po(:,1)];
                    peaks.avgb.po.time     = [peaks.avgb.po.time;  peakloc_po(:,2)];
                    peaks.avgb.po.amp       = [peaks.avgb.po.amp; amp'];
                    peaks.avgb.po.chans    = [peaks.avgb.po.chans; zeros(size(peakloc_po,1),1)+chancik];
                    peaks.avgb.po.trigcik1 = [peaks.avgb.po.trigcik1; zeros(size(peakloc_po,1),1)+trigcik1];
                    peaks.avgb.po.strigtype = [peaks.avgb.po.strigtype; zeros(size(peakloc_po,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
            end
        end
    end
end


%%% hgamma
if ~isempty(wraw.avg.poh)
    for trigcik1 = 1:length(wraw.avg.poh)
        for trigcik2 = 1:size(wraw.avg.poh{trigcik1},1)
            for chancik = 1:size(wraw.avg.poh{trigcik1},2)
                xmap = squeeze(wraw.avg.itch{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap,1);
                    v=diff(xmap(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap,2);
                    v=diff(xmap(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc = [];
                peakcounter = 0;
                
                for frqcik = 1:size(xmap,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                r = zeros(1,size(peakloc,1));
                p = zeros(1,size(peakloc,1));
                for i1=1:size(peakloc,1)
                    r(i1) = xmap(peakloc(i1,1),peakloc(i1,2));
                    p(i1)  = squeeze(wraw.avg.itcph{1}(trigcik2,chancik,peakloc(i1,1),peakloc(i1,2)));
                end
                
                xmap_po = squeeze(wraw.avg.poh{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap_po,1);
                    v=diff(xmap_po(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap_po,2);
                    v=diff(xmap_po(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc_po = [];
                peakcounter = 0;
                for frqcik = 1:size(xmap_po,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc_po(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                amp = zeros(1,size(peakloc_po,1));
                for i1=1:size(peakloc_po,1)
                    amp(i1) = xmap_po(peakloc_po(i1,1),peakloc_po(i1,2));
                end
                
                if ~isempty(peakloc)
                    peaks.avgh.itc.frq      = [peaks.avgh.itc.frq;  peakloc(:,1)];
                    peaks.avgh.itc.time     = [peaks.avgh.itc.time;  peakloc(:,2)];
                    peaks.avgh.itc.r        = [peaks.avgh.itc.r;  r'];
                    peaks.avgh.itc.p        = [peaks.avgh.itc.p;  p'];
                    peaks.avgh.itc.chans    = [peaks.avgh.itc.chans; zeros(size(peakloc,1),1)+chancik];
                    peaks.avgh.itc.trigcik1 = [peaks.avgh.itc.trigcik1; zeros(size(peakloc,1),1)+trigcik1];
                    peaks.avgh.itc.strigtype = [peaks.avgh.itc.strigtype; zeros(size(peakloc,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
                
                if ~isempty(peakloc_po)
                    peaks.avgh.po.frq      = [peaks.avgh.po.frq;  peakloc_po(:,1)];
                    peaks.avgh.po.time     = [peaks.avgh.po.time;  peakloc_po(:,2)];
                    peaks.avgh.po.amp       = [peaks.avgh.po.amp; amp'];
                    peaks.avgh.po.chans    = [peaks.avgh.po.chans; zeros(size(peakloc_po,1),1)+chancik];
                    peaks.avgh.po.trigcik1 = [peaks.avgh.po.trigcik1; zeros(size(peakloc_po,1),1)+trigcik1];
                    peaks.avgh.po.strigtype = [peaks.avgh.po.strigtype; zeros(size(peakloc_po,1),1)+wraw.avg.strigtype{trigcik1}(trigcik2)];
                end
            end
        end
    end
end

%%% surface
if ~isempty(wlfp.avg.po)
    for trigcik1 = 1:length(wlfp.avg.po)
        for trigcik2 = 1:size(wlfp.avg.po{trigcik1},1)
            for chancik = 1:size(wlfp.avg.po{trigcik1},2)
                xmap = squeeze(wlfp.avg.itc{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap,1);
                    v=diff(xmap(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap,2);
                    v=diff(xmap(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc = [];
                peakcounter = 0;
                
                for frqcik = 1:size(xmap,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                r = zeros(1,size(peakloc,1));
                p = zeros(1,size(peakloc,1));
                for i1=1:size(peakloc,1)
                    r(i1) = xmap(peakloc(i1,1),peakloc(i1,2));
                    p(i1)  = squeeze(wlfp.avg.itcp{1}(trigcik2,chancik,peakloc(i1,1),peakloc(i1,2)));
                end
                
                xmap_po = squeeze(wlfp.avg.po{trigcik1}(trigcik2,chancik,:,:));
                peaks1 = {};
                peaks2 = {};
                for frqcik = 1:size(xmap_po,1);
                    v=diff(xmap_po(frqcik,:));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks1{frqcik}=zci_pos;
                end
                for tcik = 1:size(xmap_po,2);
                    v=diff(xmap_po(:,tcik));
                    zci = find(v(:).*circshift(v(:), [-1 0]) <= 0);
                    zci_pos = zci(v(zci)>0)+1;
                    peaks2{tcik}=zci_pos;
                end
                peakloc_po = [];
                peakcounter = 0;
                for frqcik = 1:size(xmap_po,1);
                    for i1=1:length(peaks1{frqcik})
                        fpeaks = peaks2{peaks1{frqcik}(i1)};
                        for i2=1:length(fpeaks)
                            if fpeaks(i2)==frqcik
                                peakcounter =  peakcounter+1;
                                peakloc_po(peakcounter,:) = [frqcik peaks1{frqcik}(i1)];
                            end
                        end
                    end
                end
                amp = zeros(1,size(peakloc_po,1));
                for i1=1:size(peakloc_po,1)
                    amp(i1) = xmap_po(peakloc_po(i1,1),peakloc_po(i1,2));
                end
                if ~isempty(peakloc)
                    peaks.avgl.itc.frq      = [peaks.avgl.itc.frq;  peakloc(:,1)];
                    peaks.avgl.itc.time     = [peaks.avgl.itc.time;  peakloc(:,2)];
                    peaks.avgl.itc.r        = [peaks.avgl.itc.r;  r'];
                    peaks.avgl.itc.p        = [peaks.avgl.itc.p;  p'];
                    peaks.avgl.itc.chans    = [peaks.avgl.itc.chans; zeros(size(peakloc,1),1)+chancik];
                    peaks.avgl.itc.trigcik1 = [peaks.avgl.itc.trigcik1; zeros(size(peakloc,1),1)+trigcik1];
                    peaks.avgl.itc.strigtype = [peaks.avgl.itc.strigtype; zeros(size(peakloc,1),1)+wlfp.avg.strigtype{trigcik1}(trigcik2)];
                end
                
                if ~isempty(peakloc_po)
                    peaks.avgl.po.frq      = [peaks.avgl.po.frq;  peakloc_po(:,1)];
                    peaks.avgl.po.time     = [peaks.avgl.po.time;  peakloc_po(:,2)];
                    peaks.avgl.po.amp       = [peaks.avgl.po.amp; amp'];
                    peaks.avgl.po.chans    = [peaks.avgl.po.chans; zeros(size(peakloc_po,1),1)+chancik];
                    peaks.avgl.po.trigcik1 = [peaks.avgl.po.trigcik1; zeros(size(peakloc_po,1),1)+trigcik1];
                    peaks.avgl.po.strigtype = [peaks.avgl.po.strigtype; zeros(size(peakloc_po,1),1)+wlfp.avg.strigtype{trigcik1}(trigcik2)];
                end
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%% amplitude correlations

a2=['amp correlations'];
waitbar(0, waitb,a2);

ccoef.r.cm = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),2)-100;
ccoef.p.cm = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),2)-100;
ccoef.r.ch = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),2)-100;
ccoef.p.ch = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),2)-100;
ccoef.r.cl = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
ccoef.p.cl = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
ccoef.r.ml = zeros(size(wraw.cntc_po,1)+1,size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
ccoef.p.ml = zeros(size(wraw.cntc_po,1)+1,size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
ccoef.r.mh = zeros(size(wraw.cntm_po,1),size(wraw.cntm_po,2))-100;
ccoef.p.mh = zeros(size(wraw.cntm_po,1),size(wraw.cntm_po,2))-100;

if bool.coherence==1
    surfcohere.p = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
    surfcohere.r = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
    surfcohere.meandiff = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),size(wlfp.cnt_ph,1))-100;
else
    surfcohere.p = [];
    surfcohere.r = [];
    surfcohere.meandiff = [];
end

poc_mean = squeeze(mean(poc,2));

zdel = [];
if bool_arej == 1
    for chancik1 = 1:size(poc,1)
        zdel = [zdel, find(poc_mean(chancik1,:)>arej_limit)];
    end
    zdel = unique(zdel);
    poc(:,:,zdel)=[];
    slmua(:,zdel)=[];
    if ~isempty(wlfp.cnt_po)
        pol(:,:,zdel)=[];
    end
end



if ~isempty(poc) & ~isempty(slmua)
    
    %%%%% csd channels with their own MUA channels and the MUA averaged across all channels
    for chancik1 = 1:size(poc,1)
        for frqcik1 = 1:size(poc,2)
            [R,P] = corr(squeeze(poc(chancik1,frqcik1,:)),slmua(chancik1,:)','type',corrtype);
            ccoef.r.cm(chancik1,frqcik1,1)= R;
            ccoef.p.cm(chancik1,frqcik1,1)= P;
            [R,P] = corr(squeeze(poc(chancik1,frqcik1,:)),slmua(end,:)','type',corrtype);
            ccoef.r.cm(chancik1,frqcik1,2)= R;
            ccoef.p.cm(chancik1,frqcik1,2)= P;
        end
    end
    
end

if ~isempty(poc) & ~isempty(slhgamma)
    
    %%%%% csd channels with their own HGamma channels and the HGamma averaged across all channels
    for chancik1 = 1:size(poc,1)
        for frqcik1 = 1:size(poc,2)
            [R,P] = corr(squeeze(poc(chancik1,frqcik1,:)),slhgamma(chancik1,:)','type',corrtype);
            ccoef.r.ch(chancik1,frqcik1,1)= R;
            ccoef.p.ch(chancik1,frqcik1,1)= P;
            [R,P] = corr(squeeze(poc(chancik1,frqcik1,:)),slhgamma(end,:)','type',corrtype);
            ccoef.r.ch(chancik1,frqcik1,2)= R;
            ccoef.p.ch(chancik1,frqcik1,2)= P;
        end
    end
end

if ~isempty(pol) & ~isempty(poc)
    %%%%% csd channels with all LFP channels
    for chancik1 = 1:size(poc,1)
        for frqcik1 = 1:size(poc,2)
            for chancik2 = 1:size(pol,1)
                try
                    [R,P] = corr(squeeze(poc(chancik1,frqcik1,:)),squeeze(pol(chancik2,frqcik1,:))','type',corrtype);
                catch
                    [R,P] = corr(squeeze(poc(chancik1,frqcik1,:)),squeeze(pol(chancik2,frqcik1,:)),'type',corrtype);
                end
                ccoef.r.cl(chancik1,frqcik1,chancik2)= R;
                ccoef.p.cl(chancik1,frqcik1,chancik2)= P;
            end
        end
    end
end

if ~isempty(pom) && ~isempty(poh)
    %%%%% mua channels with hgamma channels at same frequency and
    %%%%% channel
    for chancik1 = 1:size(pom,1)
        for frqcik1 = 1:size(pom,2)
            
            try
                [R,P] = corr(squeeze(pom(chancik1,frqcik1,:)),squeeze(poh(chancik1,frqcik1,:))','type',corrtype);
            catch
                [R,P] = corr(squeeze(pom(chancik1,frqcik1,:)),squeeze(poh(chancik1,frqcik1,:)),'type',corrtype);
            end
            ccoef.r.mh(chancik1,frqcik1)= R;
            ccoef.p.mh(chancik1,frqcik1)= P;
        end
    end
end

if ~isempty(pol) & ~isempty(slmua)
    
    %%%%% mua channels with all LFP channels, the frequency variable is
    %%%%% lfp frequency
    for chancik1 = 1:size(slmua,1)
        for frqcik1 = 1:size(pol,2)
            for chancik2 = 1:size(pol,1)
                try
                    [R,P] = corr(squeeze(slmua(chancik1,:))',squeeze(pol(chancik2,frqcik1,:))','type',corrtype);
                catch
                    [R,P] = corr(squeeze(slmua(chancik1,:))',squeeze(pol(chancik2,frqcik1,:)),'type',corrtype);
                end
                ccoef.r.ml(chancik1,frqcik1,chancik2)= R;
                ccoef.p.ml(chancik1,frqcik1,chancik2)= P;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% coherence with surface

% if bool.coherence==1
%     a2=['coherence with surface'];
%     waitbar(0, waitb,a2);
%
%     %%%%  something wrong here, everything is too significant....
%
%     if ~isempty(wlfp.cnt_po) & ~isempty(wraw.cntc_ph)
%         for chancik1 = 1:size(wraw.cntc_ph,1)
%             for frqcik1 = 1:size(wraw.cntc_ph,2)
%                 for chancik2 = 1:size(wlfp.cnt_ph,1)
%                     phase1  = squeeze(wraw.cntc_ph(chancik1,frqcik1,:));
%                     phase2  = squeeze(wlfp.cnt_ph(chancik2,frqcik1,:));
%                     if length(phase1)>length(phase2)
%                         phase1=phase1(1:length(phase2));
%                     end
%                     if length(phase2)>length(phase1)
%                         phase2=phase2(1:length(phase1));
%                     end
%                     phdiff  = exp(1i.*(phase1-phase2));
%                     si      = mean(phdiff);
%                     mean_r  = abs(si);
%                     mean_ph = angle(si);
%                     [surfcohere.p(chancik1,frqcik1,chancik2), surfcohere.r(chancik1,frqcik1,chancik2)] = rayleigh(angle(phdiff));
%                     [surfcohere.meandiff(chancik1,frqcik1,chancik2), ang_dev] = angstat(angle(phdiff));
%                 end
%             end
%         end
%     end
%
%
% end



ccoef.zdel = zdel;
ccoef.arej_ratio = length(zdel)/(length(zdel)+size(slmua,2));

params_speaks.winsize = winsize;
params_speaks.winstep = winstep;
params_speaks.adrate = wraw.adrate;
params_speaks.corrtype = corrtype;

params_speaks.bool_ampcollect = bool_ampcollect;
params_speaks.bool_arej       = bool_arej;
params_speaks.arej_limit      = arej_limit;


close (waitb)

% save('data_spectral_peaks02_spont_clean01.mat','peaks','frq','surfcohere','ccoef','params_speaks','posall','-mat','-v7.3')

