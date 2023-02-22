function [dat] = jf_interpchan_2(dat,cfg,lay)

% interpolate bad channels
%
% parameter / example
% dat .. [n_chan,n_time]
% cfg.badchannel = {'F8','CZ'}

%cfg.badchannel=upper(cfg.badchannel);
cfg.badchannel=upper(cfg.badchannel);

addfieldtrip    
if nargin<2, cfg = []; end
if ~isfield(cfg,'elec')
    load easycap_elec_3D_59ch elec
    elec.label=upper(elec.label);
    cfg.elec = elec;
end
if ~isfield(cfg,'method'), cfg.method = 'spline'; end

tmp.trial{1} = dat;
tmp.time{1}  = 1:size(dat,2); % dummy time
tmp.label    = lay.label;
tmp.dimord   = 'chan_time';
tmp.fsample  = 1; % dummy sample freq

tmp = ft_channelrepair_JF(cfg,tmp); %fixed bug in the fieldtrip code!
dat = tmp.trial{1};
dat = dat-repmat(mean(dat),[size(dat,1),1]); % re-compute average reference
