%ana_stats_EEG_IACS
clearvars
files = dir('./EEG/imported/ICA/Cleaned/*.mat');
nchan = 59;
nfrq = 200; % interpolate power spectrum to this number of frequency bins
srb1 = 80; % short-range channel LOWER limit in mm, channels closer than this are NEIGHBORS
srb2 = 130; % short-range channel UPPER limit in mm, channels farther than this are LONG-RANGE

% Lower (X1, inclusive) and upper (X2, exclusive) frequency band limits, in Hz
D1 = 1; D2 = 4; % delta
T1 = 4; T2 = 8; % theta
A1 = 8; A2 = 12; % alpha
B1 = 12; B2 = 30; % beta
G1 = 30; G2 = 50; % gamma

% REMEMBER: EEG power transformed in this order: 1) channel-mean (within
% subject), 2) log10-scaling, 3) subject mean. Different order will change
% results. 

% Get channel distances
load('easycapM11_59ch.mat')
POS = lay.cfg.elec.chanpos; % 3D Cartesian coordinates
ED = nan(size(POS,1),size(POS,1)); % preallocation for Euclidian distances
for ich = 1:size(POS,1)
    for jch = 1:size(POS,1)
        ED(ich,jch) = sqrt((POS(ich,1)-POS(jch,1))^2 + (POS(ich,2)-POS(jch,2))^2 + ...
            (POS(ich,3)-POS(jch,3))^2);
    end
end

source = ED == 0;
neighbors = ED <= srb1 & ED > 0;
shortrange = ED > srb1 & ED < srb2;
longrange = ED >= srb2;

%% load EEG measures
for ifile = 1:length(files)
    fprintf('Now loading %s ...\n',files(ifile).name)
    load(sprintf('./EEG/imported/ICA/Cleaned/%s',files(ifile).name))
    if ~exist('LZC_all','var')
        LZC_all = nan(length(files),3); % preallocation
        LZC_allOcc = nan(length(files),3); % preallocation
        LZC_all = nan(length(files),3); % preallocation
        Oidx = contains(data.label,'O'); % occipital channels
        LZC_chans = nan(length(files),length(data.label),3); % preallocation
        DdwPLI = nan(nchan,nchan,length(files),3); % delta dwPLI
        TdwPLI = nan(nchan,nchan,length(files),3); % theta dwPLI
        AdwPLI = nan(nchan,nchan,length(files),3); % alpha dwPLI
        BdwPLI = nan(nchan,nchan,length(files),3); % beta dwPLI
        GdwPLI = nan(nchan,nchan,length(files),3); % gamma dwPLI
    end
    if ~exist('PWR_all','var')
        PWR_all = nan(size(data.PWR,1),size(data.PWR,2),size(data.PWR,3),length(files));
    end
    if ~exist('foi','var')
        foi = data.foi;
    else
        assert(all(foi == data.foi),'Frequency vectors don''t match')
    end
    
    % Loop over trials and extract measures
    for itrl = 1:3
        LZC_chans(ifile,:,itrl) = nanmean(data.LZC{itrl},2);
        LZC_all(ifile,itrl) = nanmean(data.LZC{itrl}(:));
        LZC_allOcc(ifile,itrl) = nanmean(data.LZC{itrl}(Oidx,:),[1,2]);
        if all(~isnan( data.PWR(:,:,itrl) )) % check that there are no nans in output, otherwise conclude there was not enough data for this subject/condition
            PWR_all(:,:,itrl,ifile) = data.PWR(:,:,itrl);
        end
        tmp = data.dwPLI{itrl}; 
        % For eah frequency band, do a simple average, or integral? 
        DdwPLI(:,:,ifile,itrl) = squeeze(trapz(foi(foi >= D1 & foi < D2),...
            data.dwPLI{itrl}(:,:,foi >= D1 & foi < D2),3));
        TdwPLI(:,:,ifile,itrl) = squeeze(trapz(foi(foi >= T1 & foi < T2),...
            data.dwPLI{itrl}(:,:,foi >= T1 & foi < T2),3));
        AdwPLI(:,:,ifile,itrl) = squeeze(trapz(foi(foi >= A1 & foi < A2),...
            data.dwPLI{itrl}(:,:,foi >= A1 & foi < A2),3));
        BdwPLI(:,:,ifile,itrl) = squeeze(trapz(foi(foi >= B1 & foi < B2),...
            data.dwPLI{itrl}(:,:,foi >= B1 & foi < B2),3));
        GdwPLI(:,:,ifile,itrl) = squeeze(trapz(foi(foi >= G1 & foi < G2),...
            data.dwPLI{itrl}(:,:,foi >= G1 & foi < G2),3));
        
        if ~exist('dwPLI_SR_Pre','var')
            dwPLI_SR_Pre    = nan(length(files),length(foi));
            dwPLI_SR_During = nan(length(files),length(foi));
            dwPLI_SR_Post   = nan(length(files),length(foi));
            dwPLI_LR_Pre    = nan(length(files),length(foi));
            dwPLI_LR_During = nan(length(files),length(foi));
            dwPLI_LR_Post   = nan(length(files),length(foi));
        end
        

        for ifrq = 1:length(foi)
            tmp = data.dwPLI{itrl}(:,:,ifrq);
            switch itrl
                case 1
                    dwPLI_SR_Pre(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
                    dwPLI_LR_Pre(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
                case 2
                    dwPLI_SR_During(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
                    dwPLI_LR_During(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
                case 3
                    dwPLI_SR_Post(ifile,ifrq) = squeeze(nanmean(tmp(shortrange),[1 2]));
                    dwPLI_LR_Post(ifile,ifrq) = squeeze(nanmean(tmp(longrange),[1 2]));
            end
        end
    end 
end

% make copies
DdwPLIsr = DdwPLI;
DdwPLIlr = DdwPLI;
TdwPLIsr = TdwPLI;
TdwPLIlr = TdwPLI;
AdwPLIsr = AdwPLI;
AdwPLIlr = AdwPLI;
BdwPLIsr = BdwPLI;
BdwPLIlr = BdwPLI;
GdwPLIsr = GdwPLI;
GdwPLIlr = GdwPLI;

% short-range connectiviity
DdwPLIsr(repmat(~shortrange,1,1,length(files),3)) = nan;
TdwPLIsr(repmat(~shortrange,1,1,length(files),3)) = nan;
AdwPLIsr(repmat(~shortrange,1,1,length(files),3)) = nan;
BdwPLIsr(repmat(~shortrange,1,1,length(files),3)) = nan;
GdwPLIsr(repmat(~shortrange,1,1,length(files),3)) = nan;

% long-range connectivity
DdwPLIlr(repmat(~longrange,1,1,length(files),3)) = nan;
TdwPLIlr(repmat(~longrange,1,1,length(files),3)) = nan;
AdwPLIlr(repmat(~longrange,1,1,length(files),3)) = nan;
BdwPLIlr(repmat(~longrange,1,1,length(files),3)) = nan;
GdwPLIlr(repmat(~longrange,1,1,length(files),3)) = nan;

%% Plot of dwPLI spectrum
foi_hd = 2.^[linspace(0,log2(max(foi)),nfrq)];
[~,~,sr_pre_ci] = ttest(dwPLI_SR_Pre);
[~,~,sr_during_ci] = ttest(dwPLI_SR_During);
[~,~,sr_post_ci] = ttest(dwPLI_SR_Post);
[~,~,lr_pre_ci] = ttest(dwPLI_LR_Pre);
[~,~,lr_during_ci] = ttest(dwPLI_LR_During);
[~,~,lr_post_ci] = ttest(dwPLI_LR_Post);

%% Short-range

myfigure

plot(log2(foi_hd),interp1(foi,nanmean(dwPLI_SR_Pre)',foi_hd,'spline'),'Color','r','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,sr_pre_ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,sr_pre_ci(2,:),foi_hd,'spline'))],'r','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(dwPLI_SR_During)',foi_hd,'spline'),'Color','m','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,sr_during_ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,sr_during_ci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(dwPLI_SR_Post)',foi_hd,'spline'),'Color','b','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,sr_post_ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,sr_post_ci(2,:),foi_hd,'spline'))],'b','facealpha',0.3,'EdgeColor','none')

xlim([0 log2(50)])
ylim([-0.5 1])

xlabel('Frequency (Hz)','FontSize',12)
ylabel('dwPLI','FontSize',18)
title(sprintf('Short-range (%i - %i mm) connectivity',srb1,srb2),'FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
box off
legend({'Pre','Pre 95% CI','During','During 95% CI','Post','Post 95% CI'},'fontsize',16','location','northeastoutside')
legend boxoff
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 22)
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 22)
set(gca, 'TickDir', 'out');
makefighandsome
print('-dpng','./SR_dwPLI.png')

%% Long-range

myfigure

plot(log2(foi_hd),interp1(foi,nanmean(dwPLI_LR_Pre)',foi_hd,'spline'),'Color','r','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,lr_pre_ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,lr_pre_ci(2,:),foi_hd,'spline'))],'r','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(dwPLI_LR_During)',foi_hd,'spline'),'Color','m','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,lr_during_ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,lr_during_ci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(dwPLI_LR_Post)',foi_hd,'spline'),'Color','b','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,lr_post_ci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,lr_post_ci(2,:),foi_hd,'spline'))],'b','facealpha',0.3,'EdgeColor','none')

xlim([0 log2(50)])
ylim([-0.5 1])

xlabel('Frequency (Hz)','FontSize',12)
ylabel('dwPLI','FontSize',18)
title(sprintf('Long-range (> %i mm) connectivity',srb2),'FontSize',20)
set(h,'XTick',0:1:5.5,'XTickLabel',2.^(0:1:5.5))
box off
legend({'Pre','Pre 95% CI','During','During 95% CI','Post','Post 95% CI'},'fontsize',16','location','northeastoutside')
legend boxoff
set(gca,'linewidth',3)
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 22)
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 22)
set(gca, 'TickDir', 'out');
makefighandsome
print('-dpng','./LR_dwPLI.png')



%% DELTA

%[~,:] = sort(lay.pos(:,2),'descend'); % sort channels by y-coordinate
tstrs = {'PRE','DURING','POST'};


% Separate out conditions and perform ANOVAs
dSRpre = squeeze(nanmean(DdwPLIsr(:,:,:,1),[1 2]));
dSRduring = squeeze(nanmean(DdwPLIsr(:,:,:,2),[1 2]));
dSRpost = squeeze(nanmean(DdwPLIsr(:,:,:,3),[1 2]));
testme = [dSRpre dSRduring dSRpost];
ystr = '\delta SR dwPLI';
ANOVAviolin
print('-dpng','deltaSRdwPLI.png')

% pick z-axis limits based on short-range values
chi = round(max(abs(testme(:))));
clo = chi*-1;

dLRpre = squeeze(nanmean(DdwPLIlr(:,:,:,1),[1 2]));
dLRduring = squeeze(nanmean(DdwPLIlr(:,:,:,2),[1 2]));
dLRpost = squeeze(nanmean(DdwPLIlr(:,:,:,3),[1 2]));
testme = [dLRpre dLRduring dLRpost];
ystr = '\delta LR dwPLI';
ANOVAviolin
print('-dpng','deltaLRdwPLI.png')


for itrl = 1:size(DdwPLI,4)
    myfigure2
    mypcolor(nanmean(DdwPLI(:,:,:,itrl),3))
    xticks([1:nchan]+0.5)
    xticklabels(lay.label(:))
    yticks([1:nchan]+0.5)
    yticklabels(flipud(lay.label(:)))
    caxis([clo chi])
    axis([1 nchan+1 1 nchan+1])
    makefighandsome, axis square
    set(xAX,'FontSize', 11)
    xtickangle(90)
    set(yAX,'FontSize', 11)
    title(sprintf('delta dwPLI, %s',tstrs{itrl}),'fontsize',20)
    %colormap jet
    mycolorbar
    print('-dpng',sprintf('delta_dwPLI_%s.png',tstrs{itrl}))
end



%% THETA

tSRpre = squeeze(nanmean(TdwPLIsr(:,:,:,1),[1 2]));
tSRduring = squeeze(nanmean(TdwPLIsr(:,:,:,2),[1 2]));
tSRpost = squeeze(nanmean(TdwPLIsr(:,:,:,3),[1 2]));
testme = [tSRpre tSRduring tSRpost];
ystr = '\theta SR dwPLI';
ANOVAviolin
print('-dpng','thetaSRdwPLI.png')

% pick z-axis limits based on short-range values
chi = round(max(abs(testme(:))));
clo = chi*-1;

tLRpre = squeeze(nanmean(TdwPLIlr(:,:,:,1),[1 2]));
tLRduring = squeeze(nanmean(TdwPLIlr(:,:,:,2),[1 2]));
tLRpost = squeeze(nanmean(TdwPLIlr(:,:,:,3),[1 2]));
testme = [tLRpre tLRduring tLRpost];
ystr = '\theta LR dwPLI';
ANOVAviolin
print('-dpng','thetaLRdwPLI.png')

for itrl = 1:size(TdwPLI,4)
    myfigure2
    mypcolor(nanmean(TdwPLI(:,:,:,itrl),3))
    xticks([1:nchan]+0.5)
    xticklabels(lay.label(:))
    yticks([1:nchan]+0.5)
    yticklabels(flipud(lay.label(:)))
    caxis([clo chi])
    axis([1 nchan+1 1 nchan+1])
    makefighandsome, axis square
    set(xAX,'FontSize', 11)
    xtickangle(90)
    set(yAX,'FontSize', 11)
    title(sprintf('theta dwPLI, %s',tstrs{itrl}),'fontsize',20)
    %colormap jet
    mycolorbar
    print('-dpng',sprintf('theta_dwPLI_%s.png',tstrs{itrl}))
end

%% ALPHA

aSRpre = squeeze(nanmean(AdwPLIsr(:,:,:,1),[1 2]));
aSRduring = squeeze(nanmean(AdwPLIsr(:,:,:,2),[1 2]));
aSRpost = squeeze(nanmean(AdwPLIsr(:,:,:,3),[1 2]));
testme = [aSRpre aSRduring aSRpost];
ystr = '\alpha SR dwPLI';
ANOVAviolin
print('-dpng','alphaSRdwPLI.png')

% pick z-axis limits based on short-range values
chi = round(max(abs(testme(:))));
clo = chi*-1;

aLRpre = squeeze(nanmean(AdwPLIlr(:,:,:,1),[1 2]));
aLRduring = squeeze(nanmean(AdwPLIlr(:,:,:,2),[1 2]));
aLRpost = squeeze(nanmean(AdwPLIlr(:,:,:,3),[1 2]));
testme = [aLRpre aLRduring aLRpost];
ystr = '\alpha LR dwPLI';
ANOVAviolin
print('-dpng','alphaLRdwPLI.png')

for itrl = 1:size(AdwPLI,4)
    myfigure2
    mypcolor(nanmean(AdwPLI(:,:,:,itrl),3))
    xticks([1:nchan]+0.5)
    xticklabels(lay.label(:))
    yticks([1:nchan]+0.5)
    yticklabels(flipud(lay.label(:)))
    caxis([clo chi])
    axis([1 nchan+1 1 nchan+1])
    makefighandsome, axis square
    set(xAX,'FontSize', 11)
    xtickangle(90)
    set(yAX,'FontSize', 11)
    title(sprintf('alpha dwPLI, %s',tstrs{itrl}),'fontsize',20)
    %colormap jet
    mycolorbar
    print('-dpng',sprintf('alpha_dwPLI_%s.png',tstrs{itrl}))
end

%% BETA

bSRpre = squeeze(nanmean(BdwPLIsr(:,:,:,1),[1 2]));
bSRduring = squeeze(nanmean(BdwPLIsr(:,:,:,2),[1 2]));
bSRpost = squeeze(nanmean(BdwPLIsr(:,:,:,3),[1 2]));
testme = [bSRpre bSRduring bSRpost];
ystr = '\beta SR dwPLI';
ANOVAviolin
print('-dpng','betaSRdwPLI.png')

% pick z-axis limits based on short-range values
chi = round(max(abs(testme(:))));
clo = chi*-1;

bLRpre = squeeze(nanmean(BdwPLIlr(:,:,:,1),[1 2]));
bLRduring = squeeze(nanmean(BdwPLIlr(:,:,:,2),[1 2]));
bLRpost = squeeze(nanmean(BdwPLIlr(:,:,:,3),[1 2]));
testme = [bLRpre bLRduring bLRpost];
ystr = '\beta LR dwPLI';
ANOVAviolin
print('-dpng','betaLRdwPLI.png')

for itrl = 1:size(BdwPLI,4)
    myfigure2
    mypcolor(nanmean(BdwPLI(:,:,:,itrl),3))
    xticks([1:nchan]+0.5)
    xticklabels(lay.label(:))
    yticks([1:nchan]+0.5)
    yticklabels(flipud(lay.label(:)))
    caxis([clo chi])
    axis([1 nchan+1 1 nchan+1])
    makefighandsome, axis square
    set(xAX,'FontSize', 11)
    xtickangle(90)
    set(yAX,'FontSize', 11)
    title(sprintf('beta dwPLI, %s',tstrs{itrl}),'fontsize',20)
    %colormap jet
    mycolorbar
    print('-dpng',sprintf('beta_dwPLI_%s.png',tstrs{itrl}))
end

%% GAMMA

gSRpre = squeeze(nanmean(GdwPLIsr(:,:,:,1),[1 2]));
gSRduring = squeeze(nanmean(GdwPLIsr(:,:,:,2),[1 2]));
gSRpost = squeeze(nanmean(GdwPLIsr(:,:,:,3),[1 2]));
testme = [gSRpre gSRduring gSRpost];
ystr = '\gamma SR dwPLI';
ANOVAviolin
print('-dpng','gammaSRdwPLI.png')

% pick z-axis limits based on short-range values
chi = round(max(abs(testme(:))));
clo = chi*-1;

gLRpre = squeeze(nanmean(GdwPLIlr(:,:,:,1),[1 2]));
gLRduring = squeeze(nanmean(GdwPLIlr(:,:,:,2),[1 2]));
gLRpost = squeeze(nanmean(GdwPLIlr(:,:,:,3),[1 2]));
testme = [gLRpre gLRduring gLRpost];
ystr = '\gamma LR dwPLI';
ANOVAviolin
print('-dpng','gammaLRdwPLI.png')

for itrl = 1:size(GdwPLI,4)
    myfigure2
    mypcolor(nanmean(GdwPLI(:,:,:,itrl),3))
    xticks([1:nchan]+0.5)
    xticklabels(lay.label(:))
    yticks([1:nchan]+0.5)
    yticklabels(flipud(lay.label(:)))
    caxis([clo chi])
    axis([1 nchan+1 1 nchan+1])
    makefighandsome, axis square
    set(xAX,'FontSize', 11)
    xtickangle(90)
    set(yAX,'FontSize', 11)
    title(sprintf('gamma dwPLI, %s',tstrs{itrl}),'fontsize',20)
    %colormap jet
    mycolorbar
    print('-dpng',sprintf('gamma_dwPLI_%s.png',tstrs{itrl}))
end


%% Do plotting and stats

%%% Lempel-Ziv
[P,ANOVATAB,STATS] = anova1(LZC_all)
close

myfigure
myviolin(LZC_all,'medc',[])
xticks([1 2 3])
xticklabels({'Before','During','After'})
ylim([0.2 0.35])
xlim([0.5 3.5])
makefigpretty
ylabel('Lempel-Ziv complexity')
title({'Stroboscopic stimulation',sprintf('F(%i,%i) = %1.3f, p = %1.3f',...
    ANOVATAB{2,3},ANOVATAB{3,3},ANOVATAB{2,5},P)},'fontsize',18)
print('-dpng','LZC.png')

idx = ~isnan(sum(LZC_all,2));
d = cohen_d(LZC_all(idx,1),LZC_all(idx,3));
fprintf('Pre-post Lempel-Ziv effect size: d = %1.2f\n',abs(d))

%%% Lempel-Ziv (occipital channels)
[P,ANOVATAB,STATS] = anova1(LZC_allOcc)
close

myfigure
myviolin(LZC_allOcc,'medc',[])
xticks([1 2 3])
xticklabels({'Before','During','After'})
ylim([0.15 0.4])
xlim([0.5 3.5])
makefigpretty
ylabel({'Lempel-Ziv complexity','(Occipital channels)'})
title({'Stroboscopic stimulation',sprintf('F(%i,%i) = %1.3f, p = %1.3f',...
    ANOVATAB{2,3},ANOVATAB{3,3},ANOVATAB{2,5},P)},'fontsize',18)
print('-dpng','LZC_occpital.png')

idx = ~isnan(sum(LZC_allOcc,2));
d = cohen_d(LZC_allOcc(idx,1),LZC_allOcc(idx,3));
fprintf('Pre-post occipital Lempel-Ziv effect size: d = %1.2f\n',abs(d))

%% Lempel-Ziv topoplot
load easycapM11_59ch.mat
cfg = [];
cfg.parameter = 'avg';
cfg.zlim = 'maxabs';
cfg.layout = 'easycapM11_59ch.mat';
cfg.contournum = 3;
cfg.xlim = [0,0];
cfg.layout = lay;
cfg.gridscale = 200;
cfg.comment = 'no';
cfg.markersize = 8;
%cfg.style = 'fill';
cfg.style = 'straight';

lo = 0.26;
hi = 0.30;

myfigure
dummy.var = zeros(size(LZC_chans,2),1);
dummy.dof = ones(size(LZC_chans,2),1);
dummy.time = 0;
dummy.dimord = 'chan_time';
dummy.label = lay.label;
dummy.avg = nanmean(LZC_chans(:,:,1),1)';
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('Pre Lempel-Ziv complexity','fontsize',18)
print('-dpng','LZ_topo_pre.png')

myfigure
dummy.avg = nanmean(LZC_chans(:,:,2),1)';
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('During Lempel-Ziv complexity','fontsize',18)
print('-dpng','LZ_topo_during.png')

myfigure
dummy.avg = nanmean(LZC_chans(:,:,3),1)';
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('Post Lempel-Ziv complexity','fontsize',18)
print('-dpng','LZ_topo_post.png')

%% Spectral power--average across all channels and log-scale
Pre = squeeze(log10(mean(PWR_all(:,:,1,:))));
During = squeeze(log10(mean(PWR_all(:,:,2,:))));
Post = squeeze(log10(mean(PWR_all(:,:,3,:))));

% frequency interpolation for high-res plotting
foi_hd = 2.^[linspace(log2(min(foi)),log2(max(foi)),nfrq)];

% confidence intervals
[~,~,Preci] = ttest(Pre');
[~,~,Duringci] = ttest(During');
[~,~,Postci] = ttest(Post');

myfigure

plot(log2(foi_hd),interp1(foi,nanmean(Pre,2)',foi_hd,'spline'),'r','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,Preci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,Preci(2,:),foi_hd,'spline'))],'r','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(During,2)',foi_hd,'spline'),'m','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,Duringci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,Duringci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(Post,2)',foi_hd,'spline'),'b','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,Postci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,Postci(2,:),foi_hd,'spline'))],'b','facealpha',0.3,'EdgeColor','none')

xlabel('Frequency (Hz)')
ylabel('Power [log_{10}(\muV^{2}/Hz)]')
title('Stroboscopic stimulation (channel-average)','fontsize',18)
xlim([0 log2(50)])
legend({'Pre','Pre 95% CI','During','During 95% CI','Post','Post 95% CI'},'fontsize',16','location','northeastoutside')
legend boxoff
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefighandsome
print('-dpng','PSDs.png')

%% test all frequencies in exploratory manner
Pvals = nan(1,length(foi));
for ifrq = 1:length(foi)
    [~,p] = ttest(Pre(ifrq,:),Post(ifrq,:));
    Pvals(ifrq) = p;
end

% FDR correction
Q = mafdr(Pvals,'BHFDR',true);

myfigure
plot(log2(foi),-log10(Pvals),'linewidth',2)
plot(log2(foi),-log10(Q),'linewidth',2)
plot(log2(foi),ones(1,length(foi)).*-log10(0.05),'k--')
legend({'-log10(P), raw','-log10(P), FDR','P = 0.05'},'fontsize',16','location','northeastoutside')
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(P)')
title('All channel-average')
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefighandsome
print('-dpng','Pvals.png')


%% Extract alpha power (default 8 - 12 Hz)

f1=A1; f2 = A2;
idf = foi>= f1 & foi<f2;

alpha_Pre = log10(trapz(foi(idf),10.^Pre(idf,:)));
alpha_During = log10(trapz(foi(idf),10.^During(idf,:)));
alpha_Post = log10(trapz(foi(idf),10.^Post(idf,:)));

Alpha_all = [alpha_Pre' alpha_During' alpha_Post'];

[P,ANOVATAB,STATS] = anova1(Alpha_all)
close

myfigure
myviolin(Alpha_all,'medc',[])
xticks([1 2 3])
xticklabels({'Before','During','After'})
ylim([-2 3])
xlim([0.5 3.5])
makefigpretty
ylabel('Alpha power [log_{10}(\muV^{2})]')
title({'Stroboscopic stimulation (all channels)',sprintf('F(%i,%i) = %1.3f, p = %1.3f',...
    ANOVATAB{2,3},ANOVATAB{3,3},ANOVATAB{2,5},P)},'fontsize',18)
print('-dpng','Alpha.png')

idx = ~isnan(sum(Alpha_all,2));
d = cohen_d(Alpha_all(idx,1),Alpha_all(idx,3));
fprintf('Pre-post alpha power effect size: d = %1.2f\n',abs(d))

%% Topoplot of absolute alpha power 

plotmePre = nanmean(log10(trapz(foi(idf),squeeze(PWR_all(:,idf,1,:)),2)),3);
plotmeDuring = nanmean(log10(trapz(foi(idf),squeeze(PWR_all(:,idf,2,:)),2)),3);
plotmePost = nanmean(log10(trapz(foi(idf),squeeze(PWR_all(:,idf,3,:)),2)),3);

lo = 0.5; hi = 1.5;

myfigure
dummy.avg = plotmePre;
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('Pre, log absolute alpha power','fontsize',18)
print('-dpng','PreAbsPowTopo.png')

myfigure
dummy.avg = plotmeDuring;
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('During, log absolute alpha power','fontsize',18)
print('-dpng','DuringAbsPowTopo.png')

myfigure
dummy.avg = plotmePost;
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('Post, log absolute alpha power','fontsize',18)
print('-dpng','PostAbsPowTopo.png')

%% Spectral power--average across occipital channels and log-scale

PreO = squeeze(log10(mean(PWR_all(Oidx,:,1,:))));
DuringO = squeeze(log10(mean(PWR_all(Oidx,:,2,:))));
PostO = squeeze(log10(mean(PWR_all(Oidx,:,3,:))));

% frequency interpolation for high-res plotting
foi_hd = 2.^[linspace(log2(min(foi)),log2(max(foi)),nfrq)];

% confidence intervals
[~,~,PreOci] = ttest(PreO');
[~,~,DuringOci] = ttest(DuringO');
[~,~,PostOci] = ttest(PostO');

myfigure

plot(log2(foi_hd),interp1(foi,nanmean(PreO,2)',foi_hd,'spline'),'r','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,PreOci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,PreOci(2,:),foi_hd,'spline'))],'r','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(DuringO,2)',foi_hd,'spline'),'m','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DuringOci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,DuringOci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(PostO,2)',foi_hd,'spline'),'b','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,PostOci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,PostOci(2,:),foi_hd,'spline'))],'b','facealpha',0.3,'EdgeColor','none')

xlabel('Frequency (Hz)')
ylabel('Power [log_{10}(\muV^{2}/Hz)]')
title('Stroboscopic stimulation (occipital channel-average)','fontsize',18)
xlim([0 log2(50)])
legend({'PreO','PreO 95% CI','DuringO','DuringO 95% CI','PostO','PostO 95% CI'},'fontsize',16','location','northeastoutside')
legend boxoff
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefighandsome
print('-dpng','PSDs_occipital.png')

%% test all frequencies in exploratory manner
Pvals = nan(1,length(foi));
for ifrq = 1:length(foi)
    [~,p] = ttest(PreO(ifrq,:),PostO(ifrq,:));
    Pvals(ifrq) = p;
end

% FDR correction
Q = mafdr(Pvals,'BHFDR',true);

myfigure
plot(log2(foi),-log10(Pvals),'linewidth',2)
plot(log2(foi),-log10(Q),'linewidth',2)
plot(log2(foi),ones(1,length(foi)).*-log10(0.05),'k--')
legend({'-log10(P), raw','-log10(P), FDR','P = 0.05'},'fontsize',16','location','northeastoutside')
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(P)')
title('Occipital channel-average','fontsize',18)
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefighandsome
print('-dpng','Pvals_occipital.png')


%% Extract alpha power (default 8 - 12 Hz)

f1=A1; f2 = A2;
idf = foi>= f1 & foi<f2;

alpha_PreO = log10(trapz(foi(idf),10.^PreO(idf,:)));
alpha_DuringO = log10(trapz(foi(idf),10.^DuringO(idf,:)));
alpha_PostO = log10(trapz(foi(idf),10.^PostO(idf,:)));

Alpha_allO = [alpha_PreO' alpha_DuringO' alpha_PostO'];

[P,ANOVATAB,STATS] = anova1(Alpha_allO)
close

myfigure
myviolin(Alpha_allO,'medc',[])
xticks([1 2 3])
xticklabels({'Before','DuringO','After'})
ylim([-2 4])
xlim([0.5 3.5])
makefigpretty
ylabel('Alpha power [log_{10}(\muV^{2})]')
title({'Stroboscopic stimulation (occipital channels)',sprintf('F(%i,%i) = %1.3f, p = %1.3f',...
    ANOVATAB{2,3},ANOVATAB{3,3},ANOVATAB{2,5},P)},'fontsize',18)
print('-dpng','Alpha_occipital.png')

idx = ~isnan(sum(Alpha_allO,2));
d = cohen_d(Alpha_allO(idx,1),Alpha_allO(idx,3));
fprintf('Pre-Post occipital alpha power effect size: d = %1.2f\n',abs(d))

%% Relative (i.e., normalized) spectral power--average across all channels and log-scale
Pre = squeeze(log10(mean(PWR_all(:,:,1,:))));
During = squeeze(log10(mean(PWR_all(:,:,2,:))));
Post = squeeze(log10(mean(PWR_all(:,:,3,:))));

% Undo log-10 transform and integrate
PreTot = trapz(foi,10.^Pre);
DuringTot = trapz(foi,10.^During);
PostTot = trapz(foi,10.^Post);

% Relative power normalization
PreR = log10(squeeze(mean(PWR_all(:,:,1,:)))./PreTot);
DuringR = log10(squeeze(mean(PWR_all(:,:,2,:)))./DuringTot);
PostR = log10(squeeze(mean(PWR_all(:,:,3,:)))./PostTot);

% frequency interpolation for high-res plotting
foi_hd = 2.^[linspace(log2(min(foi)),log2(max(foi)),nfrq)];

% confidence intervals
[~,~,PreRci] = ttest(PreR');
[~,~,DuringRci] = ttest(DuringR');
[~,~,PostRci] = ttest(PostR');

myfigure

plot(log2(foi_hd),interp1(foi,nanmean(PreR,2)',foi_hd,'spline'),'r','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,PreRci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,PreRci(2,:),foi_hd,'spline'))],'r','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(DuringR,2)',foi_hd,'spline'),'m','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,DuringRci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,DuringRci(2,:),foi_hd,'spline'))],'m','facealpha',0.3,'EdgeColor','none')

plot(log2(foi_hd),interp1(foi,nanmean(PostR,2)',foi_hd,'spline'),'b','LineWidth',2)
patch([log2(foi_hd) fliplr(log2(foi_hd))], [interp1(foi,PostRci(1,:),foi_hd,'spline') ...
    fliplr(interp1(foi,PostRci(2,:),foi_hd,'spline'))],'b','facealpha',0.3,'EdgeColor','none')

xlabel('Frequency (Hz)')
ylabel('Relative power [log_{10}(1/Hz)]')
title('Stroboscopic stimulation (channel-average)','fontsize',18)
xlim([0 log2(50)])
legend({'PreR','PreR 95% CI','DuringR','DuringR 95% CI','PostR','PostR 95% CI'},'fontsize',16','location','northeastoutside')
legend boxoff
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefighandsome
print('-dpng','PSDs_relative.png')

%% test all frequencies in exploratory manner
Pvals = nan(1,length(foi));
for ifrq = 1:length(foi)
    [~,p] = ttest(PreR(ifrq,:),PostR(ifrq,:));
    Pvals(ifrq) = p;
end

% FDR correction
Q = mafdr(Pvals,'BHFDR',true);

myfigure
plot(log2(foi),-log10(Pvals),'linewidth',2)
plot(log2(foi),-log10(Q),'linewidth',2)
plot(log2(foi),ones(1,length(foi)).*-log10(0.05),'k--')
legend({'-log10(P), raw','-log10(P), FDR','P = 0.05'},'fontsize',16','location','northeastoutside')
legend boxoff
xlabel('Frequency (Hz)')
ylabel('-log_{10}(P)')
title('All channel-average (rel. power)')
set(h,'XTick',-1:1:5.5,'XTickLabel',2.^(-1:1:5.5))
makefighandsome
print('-dpng','Pvals_relative.png')


%% Extract alpha power (default 8 - 12 Hz)

f1=A1; f2 = A2;
idf = foi>= f1 & foi<f2;

alpha_PreR = log10(trapz(foi(idf),10.^PreR(idf,:)));
alpha_DuringR = log10(trapz(foi(idf),10.^DuringR(idf,:)));
alpha_PostR = log10(trapz(foi(idf),10.^PostR(idf,:)));

Alpha_allR = [alpha_PreR' alpha_DuringR' alpha_PostR'];

[P,ANOVATAB,STATS] = anova1(Alpha_allR)
close

myfigure
myviolin(Alpha_allR,'medc',[])
xticks([1 2 3])
xticklabels({'Before','During','After'})
ylim([-1.5 0.5])
xlim([0.5 3.5])
makefigpretty
ylabel('Alpha relative power')
title({'Stroboscopic stimulation (all channels, rel. power)',sprintf('F(%i,%i) = %1.3f, p = %1.3f',...
    ANOVATAB{2,3},ANOVATAB{3,3},ANOVATAB{2,5},P)},'fontsize',18)
print('-dpng','Alpha_relative.png')

idx = ~isnan(sum(Alpha_allR,2));
d = cohen_d(Alpha_allR(idx,1),Alpha_allR(idx,3));
fprintf('Pre-Post relative alpha power effect size: d = %1.2f\n',abs(d))

%% Topoplot of relative alpha power 

% Get total power
TotPre = nanmean(log10(trapz(foi,squeeze(PWR_all(:,:,1,:)),2)),3);
TotDuring = nanmean(log10(trapz(foi,squeeze(PWR_all(:,:,2,:)),2)),3);
TotPost = nanmean(log10(trapz(foi,squeeze(PWR_all(:,:,3,:)),2)),3);

plotmePreR = log10(10.^plotmePre./10.^TotPre); % total power normalization
plotmeDuringR = log10(10.^plotmeDuring./10.^TotDuring); % total power normalization
plotmePostR = log10(10.^plotmePost./10.^TotPost); % total power normalization

lo = -0.7; hi = -0.4;

myfigure
dummy.avg = plotmePreR;
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('Pre, log relative alpha power','fontsize',18)
print('-dpng','PreRelPowTopo.png')

myfigure
dummy.avg = plotmeDuringR;
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('During, log relative alpha power','fontsize',18)
print('-dpng','DuringRelPowTopo.png')

myfigure
dummy.avg = plotmePostR;
ft_topoplotER_JF(cfg,dummy);
caxis([lo hi])
mycolorbar
title('Post, log relative alpha power','fontsize',18)
print('-dpng','PostRelPowTopo.png')

%% Write results to CSV

% Extract subject IDs
SIDs = [];
for ifile = 1:length(files)
    SIDs = [SIDs; str2double(files(ifile).name(1:6))];
end

% build table
T = table(SIDs,'VariableNames',{'SubjectID'});
T2 = T;

for ich = 1:length(data.label)
    for ifrq = 1:length(foi)
        for icnd = 1:3
            switch icnd
                case  1
                    eval(sprintf('T.%s_PRE_%i_Hz_log_power = log10(squeeze(PWR_all(ich,ifrq,icnd,:)));',...
                        data.label{ich},ifrq));
                case 2
                    eval(sprintf('T.%s_DURING_%iHz_log_power = log10(squeeze(PWR_all(ich,ifrq,icnd,:)));',...
                        data.label{ich},ifrq));
                case 3
                    eval(sprintf('T.%s_POST_%i_Hz_log_power = log10(squeeze(PWR_all(ich,ifrq,icnd,:)));',...
                        data.label{ich},ifrq));
            end
        end
    end
    for icnd = 1:3
        switch icnd
            case  1
                eval(sprintf('T2.%s_PRE_LZC = LZC_chans(:,ich,icnd);',...
                    data.label{ich}));
            case 2
                eval(sprintf('T2.%s_DURING_LZC = LZC_chans(:,ich,icnd);',...
                    data.label{ich}));
            case 3
                eval(sprintf('T2.%s_POST_LZC = LZC_chans(:,ich,icnd);',...
                    data.label{ich}));
        end
    end
end

writetable(T,'./IASC_Stroboscope_EEG_Power.csv')
writetable(T2,'./IASC_Stroboscope_EEG_LZC.csv')

%% Now make a table with just channel-averaged EEG measures

T3 = table(Alpha_all(:,1),Alpha_all(:,2),Alpha_all(:,3),...
    Alpha_allO(:,1),Alpha_allO(:,2),Alpha_allO(:,3),...
    Alpha_allR(:,1),Alpha_allR(:,2),Alpha_allR(:,3),...
    LZC_all(:,1),LZC_all(:,2),LZC_all(:,3),...
    LZC_allOcc(:,1),LZC_allOcc(:,2),LZC_allOcc(:,3),...
    'VariableNames',{'LogAbsAlphaPower_PRE','LogAbsAlphaPower_DURING','LogAbsAlphaPower_POST',...
    'LogOccAbsAlphaPower_PRE','LogOccAbsAlphaPower_DURING','LogOccAbsAlphaPower_POST'...
    'LogRelAlphaPower_PRE','LogRelAlphaPower_DURING','LogRelAlphaPower_POST',...
    'LZC_PRE','LZC_DURING','LZC_POST',...
    'Occipital_LZC_PRE','Occipital_LZC_DURING','Occipital_LZC_POST'});

writetable(T3,'./IASC_Stroboscope_ch_avg_EEG_measures.csv')





    

