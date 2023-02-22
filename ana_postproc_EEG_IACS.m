clearvars
dbstop if error
files = dir('./EEG/imported/ICA/*.mat');
nch = 59; % target number of channels (after removing peripheral channels)
load BrainCap64chLayFile lay
lay64 = lay; clear lay
load 'easycapM11_59ch.mat' lay
lay59 = lay; clear lay
pchs = {'FT10','FT9','IZ','TP10','TP9'};
%myPool = parpool(6); % for dwPLI
winsize = 10; % length of the connectivity windows in seconds 
overlap = true; % do 50% overlapping windows for connectivity? 
        
% load the neighbour structure for 59 channels
n59 = load('easycap_neighb_59ch.mat');

start_ndx = 1;

if start_ndx ~= 1, warning('ATTENTION: the starting index is set to %i\n',start_ndx), end

for ifile = start_ndx:length(files)
    fprintf('Now loading %s ...\n',files(ifile).name)
    EEG = load(sprintf('./EEG/imported/ICA/%s',files(ifile).name));
    data = EEG.data; clear EEG
    
    % Apply ICA
    fulllbl = data.label; % save labels
    data.label = upper(data.label(data.cfg.ica.userawch)); % remove labels of channels not used for ICA
    idcompbad = [data.cfg.ica.ctype.muscle.num' data.cfg.ica.ctype.other.num'...
        data.cfg.ica.ctype.blink.num' data.cfg.ica.ctype.ocul.num'];
    
    % ica weights matrix
    weights = data.cfg.ica.a(:,idcompbad)*data.cfg.ica.w(idcompbad,:);
    
    % loop over trials
    for itrl = 1:length(data.trial)
        if itrl == 1, data.labelold = data.label; end

        ICA_cleaned = data.trial{itrl}(data.cfg.ica.userawch,:) - ...
            weights*data.trial{itrl}(data.cfg.ica.userawch,:);
        fprintf('     ICA subtracted \n')
        
        % interpolate bad channels
        data2 = data;
        data2.trial{itrl} = ICA_cleaned;
        data2.cfg.raw = data.cfg.raw;
        
        % fill in the missing channel as nan
        data3 = data2;
        data3.trial{itrl} = [];
        data3 = rmfield(data3,'label');
        for ielec=1:length(data3.elec.label)
            idx = find(strcmpi(data3.elec.label{ielec},data2.labelold));
            if ~isempty(idx)
                data3.trial{itrl}(ielec,:) = data2.trial{itrl}(idx,:);
                data3.label{ielec} = data3.elec.label{ielec};
            else
                data3.trial{itrl}(ielec,:) = nan(1,size(data2.trial{itrl},2));
                data3.label{ielec} = data3.elec.label{ielec};
            end
        end
        
        % Remove channels that aren't in lay file
        IS = setdiff(upper(data3.label),upper(lay64.label));
        idx = find(ismember(upper(data3.label),IS));
        data3.trial{itrl}(idx,:) = [];
        data3.elec.label(idx) = [];
        data3.elec.pnt(idx,:) = [];
        data3.label(idx) = []; 
        badch = data.cfg.raw.ctype.bad.label;
        if ~isempty(badch{1}) && length(setdiff(upper(badch),upper(pchs))) > sqrt(nch)
            warning('ATTENTION: This dataset contains a large number (n = %i) of bad channels, consider QC''ing. Press any key to continue\n',length(setdiff(upper(badch),upper(pchs))))
            pause
        end
        
        % overwrite data
        data = data3;
        
        % Remove peripheral channels 
        rmvidx = find(contains(data.label,pchs,'ignorecase',true));
        data.trial{itrl}(rmvidx,:) = [];
        data.label(rmvidx) = [];
        data.elec.pnt(rmvidx,:) = [];
        data.elec.label(rmvidx) = [];
                
        if ~isempty(badch)
            if itrl == 1
                clear cfg
                cfg.badchannel=badch;
                cfg.missingchannel = [];
                cfg2 = cfg;
                cfg2.method = 'template';
                cfg2.layout = 'easycapM11_59ch.mat';
                cfg2.elecfile = 'easycap_elec_3D_59ch.mat';
                cfg2.template = 'easycap_neighb_59ch.mat';
            end

            data4 = data;
            data4.trial = [];
            data4.time = [];
            data4.time{1} = data.time{itrl};
            data4.trial{1} = data.trial{itrl}; % so we don't get inconsistant number of channels error below
            neighbours = ft_prepare_neighbours(cfg2, data4);
                        
            cfg.layout = 'easycapM11_59ch.mat';
            cfg.method = 'spline';
            cfg.elecfile = 'easycap_elec_3D_59ch.mat';
            cfg.neighbours = neighbours;
            tmp = jf_interpchan_2(data.trial{itrl},cfg,lay59); % also re-computes average reference
            data.trial{itrl} = tmp; % overwrite with interpolated channels                
                       
            if itrl == 1
                data.cfg.raw.ctype.bad.num = [];
                data.cfg.raw.ctype.bad.label = cell(1,1);
                data.cfg.raw.ctype.good.num = 1:length(data.label);
                data.cfg.raw.ctype.good.label = data.label;
            end
        end
    end
    
    % Surface Laplacian/CSD
       
    % Surface laplacian (for EEG connectivity)
    cfg_sl.method = 'hjorth';
    cfg_sl.neighbours = n59.neighbours;
    data_sl = ft_scalpcurrentdensity(cfg_sl, data); % CSD/surface laplacian
    
    % copy over data fields
    for itrl = 1:size(data_sl.trial,2)
        data.CSD{itrl} = data_sl.trial{itrl};
    end
       
    % update other fields to reflect correct number of channels
    [~,rmvidx] = setdiff(upper(data.cfg.raw.alllabels),upper(lay59.label)); % find the channels to remove
    data.cfg.raw.alllabels(rmvidx) = [];
    data.cfg.raw.allpos(rmvidx,:) = [];

    assert(size(data.trial{1},1)==nch,'Wrong number of chanels')
    assert(length(data.label)==nch,'Wrong number of chanels')

    % Uncomment below to already save before computing EEG measures 
%     % save the ICA subtracted data
%     outname = sprintf('./EEG/imported/ICA/Cleaned/%s_cleaned.mat',files(ifile).name);
%     save(outname,'data','-v7.3')

    %% Compute EEG measures

    % Set bad data to NaNs for purpose of taking freq transform
    if isfield(data.cfg.dattype,'bigmusclemov')
        badart = [data.cfg.dattype.bigmusclemov; data.cfg.dattype.muscle; ...
            data.cfg.dattype.bad; data.cfg.dattype.sleep];
    else
        badart = [data.cfg.dattype.muscle; data.cfg.dattype.bad; data.cfg.dattype.sleep];
    end
    
    badart(badart==0) = 1; % do this so no 0 indicies that give bugs
    for itrl = 1:length(data.trial)
        for iart = 1:size(badart,1)
            idpt = badart(iart,1) - data.sampleinfo(itrl,1) + 1;
            if idpt > 0 && idpt+diff(badart(iart,:)) <= size(data.trial{itrl},2)
                data.trial{itrl}(:,idpt:idpt+diff(badart(iart,:))) = nan;
            end
        end
    end
        
    %%% LEMPEL-ZIV COMPLEXITY %%%
    fprintf('     Computing Lempel-Ziv complexity ...\n')
    tic
    LZC = cell(1,length(data.trial));
    for itrl = 1:length(data.trial)
        datamat = data.trial{itrl};
        n_win = data.fsample*10; % compute LZC in 10 s windows
        % pad with nans 
        datamat = [datamat nan(size(datamat,1),n_win-mod(size(datamat,2),n_win))];
        % reshape to channels x samples x windows
        datamat = reshape(datamat,size(datamat,1),n_win,size(datamat,2)/n_win);
        LZC{itrl} = nan(size(datamat,1),size(datamat,3)); % preallocation
        for iwn = 1:size(datamat,3)
            if ~any(isnan(datamat(:,:,iwn)),[1 2]) % if there are no NaNs in this window
                for ich = 1:size(datamat,1)
                    LZC{itrl}(ich,iwn) = LZ76(datamat(ich,:,iwn)>=median(datamat(ich,:,iwn))) ...
                        /(length(datamat(ich,:,iwn))/log2(length(datamat(ich,:,iwn)))); % normalization for data length
                end
            end
        end
    end
    
    data.LZC = LZC;
    fprintf('          %i seconds elapsed\n',round(toc))
    
    %%% SPECTRAL POWER %%%
    fprintf('     Computing spectral power (with Morlet wavelets) ...\n')
    tic
    for itrl = 1:length(data.trial)
        datamat = data.trial{itrl};
        cfg = [];
        cfg.fsample = data.fsample;
        cfg.verbose = false;
        cfg.foi_start = 1; % lowest frequency
        cfg.foi_end = 50; % highest frequency
        [pow,pow_full,pow_median,pow_var,n,unit,foi_target,foi_delta,foi,...
            cfg,w_lin2log] = ro_freq_wavelet_JF(datamat,cfg);
        if ~exist('PWR','var')
            PWR = nan(size(pow,1),size(pow,2),3);
        end
        PWR(:,:,itrl) = pow;
    end
    
    data.PWR = PWR;
    data.foi = foi; % frequency vector
    fprintf('          %i seconds elapsed\n',round(toc))
    
% UNCOMMENT BELOW FOR wSMI 
%%
%     %%% wSMI %%%
%     fprintf('     Computing functional connectivity with wSMI (weighted symbolic mutual information)\n')
%     tic
%     % partion into n sec windows and compute wSMI
%     
%     nwin = 5;
%     new_fsample = 125; % downsample to this frequency 
%     win_len = new_fsample*nwin; % n sec windows
%     
%     % BAND (tau): 32-80 Hz (tau = 4 ms), 16-40 Hz (tau = 8 ms), 8-20 Hz 
%     % (tau = 16 ms), 4-10 Hz (tau = 32 ms), 2-5 Hz (tau =  64 ms) and 1-2.5 Hz 
%     %(tau = 128 ms) [from   Bourdillon et al. 2019 BioRxiv) 
%     
%     for itrl = 1:size(data.trial,2)
%         datamat = data.CSD{itrl};
%         cfg = [];
%         cfg.sf = new_fsample;
%         cfg.verbose = false;
%         cfg.kernel = 3;
%         cfg.taus = round([8 16 32 64].*(new_fsample/1000)); % covert tau from ms to samples (e.g., 32 ms becomes 8 samples at fs = 250 Hz)
%         cfg.taus_ms = [8 16 32 64];
%         cfg.chan_sel = 1:size(datamat,1); % use all channels
%         cfg.data_sel = 1:win_len; % use all samples
% 
%         gcount = 0;
%         bcount = 0;
% 
%         % windowing parameters
%         % 1.0 --> no overlap; 0.5 --> 50 percent
%         window_shift = 1.0; % do x% overlap
%         n_shift = round(win_len*window_shift);
%         for isection = 1:n_shift:size(datamat,2)-win_len+1
%             section = double(datamat(:,isection:isection+win_len-1));
% 
%             if ~any(isnan(section(:))) % IF there are no NaNs in the window
%                 [cfg_out] = ro_wSMI(section,cfg);
% 
%                 for itau = 1:length(cfg.taus)
%                     % create variable if it doesn't already exist
%                     %%% wSMI %%%
%                     if ~exist(sprintf('wSMI_%i',cfg.taus(itau)),'var')
%                         eval(sprintf('wSMI_%i = [];',cfg.taus(itau)));
%                     end
%                     eval(sprintf('wSMI_%i = cat(3,cfg_out.wSMI{itau},wSMI_%i);',...
%                         cfg.taus(itau),cfg.taus(itau)))
%                 end
%                 gcount = gcount + 1;
%                 %fprintf('%2.2f percent complete \n',(isection/(size(datamat,2)-win_len+1))*100)
%                 %toc
%             else
%                 bcount = bcount + 1;
%                 %fprintf('Window contains NaNs, continuing to next window ...\n')
%             end
%         end
% 
%         wSMIout.cfg = cfg_out;
%         wSMIout.window_shift = window_shift; % window overlap
% 
%         for itau = 1:length(cfg.taus)
%             eval(sprintf('wSMIout.wSMI_%i = wSMI_%i;',cfg.taus(itau),cfg.taus(itau)));
%         end
% 
%         wSMIout.prc_win_used = gcount/(gcount+bcount)*100;
%         wSMIout.win_len_used = gcount;
%         wSMIout.win_len      = nwin;
% 
%         data.wSMI{itrl} = wSMIout;
%     end
%     toc
%     fprintf('          %i seconds elapsed\n',round(toc))     

%%
    %% dwPLI %%%
    fprintf('     Computing dwPLI ...\n')
    tic
    % start with large window and make smaller later
    n_win = data.fsample*winsize; % compute dwPLI in n s windows
    dwPLI = cell(1,length(data.CSD));
    for itrl = 1:length(data.CSD)
        fprintf('     Condition %i ...\n',itrl)
        datamat = data.CSD{itrl};
        % pad with nans
        datamat = [datamat nan(size(datamat,1),n_win-mod(size(datamat,2),n_win))];

        switch overlap
            case true
                datamat_3d = reshape(datamat,size(datamat,1),n_win,size(datamat,2)/n_win);

                % create second datamat to get 50% overlapping windows
                datamat2 = datamat(:,(n_win/2)+1:end-(n_win/2));
                datamat2_3d = reshape(datamat2,size(datamat2,1),n_win,size(datamat2,2)/n_win);

                % check that the second matrix has appropriate dimensions
                assert(size(datamat_3d,1)==size(datamat2_3d,1))
                assert(size(datamat_3d,2)==size(datamat2_3d,2))
                assert(size(datamat_3d,3)==size(datamat2_3d,3)+1)

                % interleave the two 3D arrays so we acheive 50% window overlap

                data_3d = nan(size(datamat_3d,1),size(datamat_3d,2),size(datamat_3d,3)*2-1);
                for ilay = 1:size(data_3d,3)
                    if mod(ilay,2) == 1 % if it's an odd numbered layer
                        data_3d(:,:,ilay) = datamat_3d(:,:,ceil(ilay/2));
                    else % if it's an even numbered layer
                        data_3d(:,:,ilay) = datamat2_3d(:,:,ilay/2);
                        % check that they really overlap
                            % if this and the last layer don't have nans
                        if ~any(isnan(data_3d(:,:,ilay)),[1 2]) && ~any(isnan(data_3d(:,:,ilay-1)),[1 2])
                            assert(all(data_3d(:,1:n_win/2,ilay) == data_3d(:,n_win/2+1:n_win,ilay-1),[1 2]))
                        end
                    end
                end

                [~,dwpli] = ro_wPLI(data_3d,data.fsample,data.foi);
            case false
                datamat_3d = reshape(datamat,size(datamat,1),n_win,size(datamat,2)/n_win);
                [~,dwpli] = ro_wPLI(datamat_3d,data.fsample,data.foi);
        end
        dwPLI{itrl} = dwpli;
    end
    data.dwPLI = dwPLI;
    fprintf('          %i seconds elapsed\n',round(toc))

%% 
    % save with EEG measures attached to data structure
    outname = sprintf('./EEG/imported/ICA/Cleaned/%s_cleaned_meas.mat',files(ifile).name(1:end-4));
    save(outname,'data','-v7.3')
    
end
    



