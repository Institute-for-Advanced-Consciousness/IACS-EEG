clearvars

% add folders to path
addpath('.\fieldtrip-20170827')
ft_defaults

% Condition to load, edit as necessary
condition = 'Pre Experience';
nchan = 64; % expected number of EEG channels

% duration of pre and post resting periods in *seconds* (
% can be used to set a timestamp if beginning or end is missing)
nsecpp = 300; 

hpfreq = 0.5; % highpass filter frequency in Hz
% lpfreq can also be set to 70 or 80 to look at high gamma, just need to
% add notch filter if you do this

lpfreq = 50; % lowpass filter frequency in Hz


% Find all XDF files in directory
FList = ReadXdfFileNames('./EEG/xdf');
start_ndx = 1;

if start_ndx ~= 1
    warning('start_ndx set to %i',start_ndx)
end

for ifile = start_ndx:length(FList)
    try
        tmp = strfind(FList(ifile),'1'); % find beginning of subject ID
        A = tmp{1}(1);
        B = tmp{1}(1) + 5;
        tmpstr = FList(ifile);
        SIDstr = tmpstr{1}(A:B);
        SID = str2double(tmpstr{1}(A:B));
        assert(~isempty(SID) & ~isnan(SID),'SID not found')
    catch
        tmp = strfind(FList(ifile),'P'); % find beginning of subject ID
        A = tmp{1}(1)+1;
        B = tmp{1}(1) + 3;
        tmpstr = FList(ifile);
        SIDstr = tmpstr{1}(A:B);
        SID = str2double(tmpstr{1}(A:B));
        assert(~isempty(SID) & ~isnan(SID),'SID not found')
    end
    fprintf('Now loading %s\n',FList{ifile})
    
    % Try to load the file, but skip if it's corrupt
    try
        [OUTDATA,fileheader] = load_xdf(FList{ifile}); % load XDF file
    catch
        fprintf('     This file couldn''t be loaded (probably corrupt), skipping this one ...')
        continue
    end
    
    Pre = [];
    Post = [];
    During = [];
    eegidx = nan;
    expidx = nan; 
    for icell = 1:length(OUTDATA)
        if iscell(OUTDATA{icell}.time_series)
            if contains(OUTDATA{icell}.time_series,'Pre Experience','IgnoreCase',true)
                assert(length(OUTDATA{icell}.time_stamps)==1,'More than one timestamp detected for Pre start/stop')
                    Pre = [Pre OUTDATA{icell}.time_stamps];
                    preidx = icell; % save this for later, in case we need to come back to find this
            elseif contains(OUTDATA{icell}.time_series,'Post Experience','IgnoreCase',true)
                assert(length(OUTDATA{icell}.time_stamps)==1,'More than one timestamp detected for Post start/stop')
                    Post = [Post OUTDATA{icell}.time_stamps];
                    postidx = icell; % save this for later, in case we need to come back to find this
            elseif contains(OUTDATA{icell}.time_series,'Meditation','IgnoreCase',true) || ...
                    contains(OUTDATA{icell}.time_series,'Minute Experience','IgnoreCase',true)
                assert(length(OUTDATA{icell}.time_stamps)==1,'More than one timestamp detected for During start/stop')
                    During = [During OUTDATA{icell}.time_stamps];
                    expidx = icell; % save this for later, in case we need to come back to find the duration
            end
        elseif strcmp(OUTDATA{icell}.info.type,'EEG') && str2double(OUTDATA{icell}.info.channel_count) >= nchan
            if isnan(eegidx)
                eegidx = icell;
            else
                error('Two cell array elements look like they correspond to EEG')
            end
        end      
    end
    
    if length(Pre) == 1 % if we only found one time stamp for Pre
        descript = OUTDATA{expidx}.time_series;
        assert(iscell(descript),'Wrong index found for the Pre condition')
        nsecpp = 300;
        if contains(descript,'Started') % if we have the start time stamp
            Pre = [Pre Pre(1)+nsecpp];
        elseif contains(descript,'Ended') % if we have the end time stamp
            Pre = [Pre(1)-nsecpp Pre(1)];
        else
            error('Only one timestamp found for Pre, can''t determine which one it is')
        end
    end
        
    
    if length(During) == 1 % if we only found one time stamp for During
        descript = OUTDATA{expidx}.time_series;
        assert(iscell(descript),'Wrong index found for the During condition')
        where = strfind(descript{1},'Minute');
        duration = str2num(descript{1}(1:where-2)); % duration of light flash/meditation
        nsec = duration*60; % number of seconds
        During = [During During(1)+nsec];
    end  
    
    if length(Post) == 1 % if we only found one time stamp for Post
        descript = OUTDATA{expidx}.time_series;
        assert(iscell(descript),'Wrong index found for the Post condition')
        if contains(descript,'Started') % if we have the start time stamp
            Post = [Post Post(1)+nsecpp];
        elseif contains(descript,'Ended') % if we have the end time stamp
            nsecpp = 300;
            Post = [Post(1)-nsecpp Post(1)];
        else
            error('Only one timestamp found for Post, can''t determine which one it is')
        end
    end
    

    assert(~isnan(eegidx),'Failed to find EEG data in the loaded data sturcture')
    
    try
        assert(length(Pre)==2,'Failed to return two Pre timestamps, one for start and one for stop')
        Pre = sort(Pre,'ascend'); % order timestamps chronologically
        assert(length(During)==2,'Failed to return two During timestamps, one for start and one for stop')
        During = sort(During,'ascend'); % order timestamps chronologically
        assert(length(Post)==2,'Failed to return two Post timestamps, one for start and one for stop')
        Post = sort(Post,'ascend'); % order timestamps chronologically
    catch
        fprintf('Error with timestamps, skipping this one ...\n')
        continue
    end
    
    % Create Fieldtrip data structure
    data.fsample  = str2num(OUTDATA{eegidx}.info.nominal_srate);
    % Pre
    [~,START]     = min(abs(OUTDATA{eegidx}.time_stamps - Pre(1)));
    [~,STOP]      = min(abs(OUTDATA{eegidx}.time_stamps - Pre(end)));
    data.trial{1} = double(OUTDATA{eegidx}.time_series(:,START:STOP));
    data.time{1}  = linspace(0,length(data.trial{1})/data.fsample,length(data.trial{1}));
    % During
    [~,START]     = min(abs(OUTDATA{eegidx}.time_stamps - During(1)));
    [~,STOP]      = min(abs(OUTDATA{eegidx}.time_stamps - During(end)));
    data.trial{2} = double(OUTDATA{eegidx}.time_series(:,START:STOP));
    data.time{2}  = linspace(0,length(data.trial{2})/data.fsample,length(data.trial{2}));
    % Post
    [~,START]     = min(abs(OUTDATA{eegidx}.time_stamps - Post(1)));
    [~,STOP]      = min(abs(OUTDATA{eegidx}.time_stamps - Post(end)));
    data.trial{3} = double(OUTDATA{eegidx}.time_series(:,START:STOP));
    data.time{3}  = linspace(0,length(data.trial{3})/data.fsample,length(data.trial{3}));
    
    % Get channel labels
    load elec64
    % Fix channel labels
    elec64.label{find(strcmpi(elec64.label,'Io'))} = 'Iz'; 
    elec64.label{find(strcmpi(elec64.label,'FT9'))} = 'TP9x'; % placeholder
    elec64.label{find(strcmpi(elec64.label,'FT10'))} = 'TP10x'; % placeholder
    elec64.label{find(strcmpi(elec64.label,'TP9'))} = 'FT9'; 
    elec64.label{find(strcmpi(elec64.label,'TP10'))} = 'FT10'; 
    elec64.label{find(strcmpi(elec64.label,'TP9x'))} = 'TP9'; 
    elec64.label{find(strcmpi(elec64.label,'TP10x'))} = 'TP10'; 
    data.elec = elec64;
    data.label = elec64.label;
    
    % Load 2D channel locations (topo projection)
    load BrainCap2DwLabels
        
    % Remove any channels that aren't in the lay file
    exclude = setdiff(upper(data.label),upper(BrainCapLabels));
    exidx = [];
    for ich = 1:length(exclude)
        exidx = [exidx find(strcmpi(data.label,exclude{ich}))];
    end
    
    data.elec.pnt(exidx,:) = [];
    data.elec.label(exidx) = [];
    data.label(exidx) = [];
    for itrl = 1:size(data.trial,1)
        data.trial{itrl}(exidx,:) = [];
    end
        
    % Sanity check
    assert(length(data.label)==length(BrainCapLabels),'Wrong number of channels')
    assert(isempty(setdiff(data.label,BrainCapLabels)),'Channels don''t match')
    NDX = sortref(BrainCapLabels,data.label); % sort A according to B
    if ~all(NDX==1:length(data.label))
        data.cfg.raw.alllabels = data.cfg.raw.alllabels(NDX);
        data.cfg.raw.allpos = data.cfg.raw.allpos(NDX);
    end   
        
    % Fill configuration field
    data.cfg.info = OUTDATA{eegidx}.info;
    data.cfg.raw.alllabels = BrainCapLabels;
    data.cfg.ctypelabels = {'good','bad','physio'};
    data.cfg.ctypemainlabel = 'good';
    data.cfg.unit = 'uV';
    data.cfg.raw.allpos = BrainCap2D;
    data.cfg.raw.haspos = true(length(data.label),1);
    
    %%%% Define filters (credit Pilar Garces for filter code) %%%%
    
    %%%%%% WANT TO ADD A NOTCH FILTER??? %%%%%%%
    % To add a notch filter, use this syntax
    % B = fir1(N,Wn,'stop') will design a bandstop filter.
    
    % Wn is filter frequency in terms of Nyquist freq, e.g., if fs = 250
    % Hz, then 60/(250/2) = 0.48 as your value for Wn
    
%      B = fir1(N,Wn) designs an N'th order lowpass FIR digital filter
%     and returns the filter coefficients in length N+1 vector B.
%     The cut-off frequency Wn must be between 0 < Wn < 1.0, with 1.0
%     corresponding to half the sample rate.  The filter B is real and
%     has linear phase.  The normalized gain of the filter at Wn is
%     -6 dB.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fNyquist = data.fsample/2;     

    % Highpass
    fsample = data.fsample;
    TRANSWIDTHRATIO = 0.25; % Constant in pop_eegfiltnew (the eeglab function)
    df = min([max([hpfreq * TRANSWIDTHRATIO 2]) hpfreq]);
    edgeArray = hpfreq;
    g.filtorder = 3.3 / (df / fsample); % Hamming window
    g.filtorder = ceil(g.filtorder / 2) * 2; % Filter order must be even.
    winArray = window('hamming', g.filtorder + 1);
    cutoffArray = edgeArray + -df / 2;
    bhp = fir1(g.filtorder, cutoffArray / fNyquist, 'high', winArray);
    
    % Lowpass
    fsample = data.fsample;
    TRANSWIDTHRATIO = 0.25; %Constant in pop_eegfiltnew
    df = min([max([lpfreq * TRANSWIDTHRATIO 2]) lpfreq]);
    edgeArray = lpfreq;
    g.filtorder = 3.3 / (df / fsample); % Hamming window
    g.filtorder = ceil(g.filtorder / 2) * 2; % Filter order must be even.
    winArray = window('hamming', g.filtorder + 1);
    cutoffArray = edgeArray + -df / 2;
    blp = fir1(g.filtorder, cutoffArray / fNyquist, 'low', winArray);
    
    % loop over data trials and filter
    for itrial = 1:length(data.trial)    
        tmp = data.trial{itrial};
        tmp2 = filtfilt(blp,1,double(tmp')); % lowpass filter
        tmp3 = filtfilt(bhp,1,tmp2)'; % highpass filter
        data.trial{itrial} = single(tmp3);
    end
    
    outname = sprintf('./EEG/imported/%s_imported.mat',SIDstr);
    save(outname,'data','-v7.3')
    
end
    