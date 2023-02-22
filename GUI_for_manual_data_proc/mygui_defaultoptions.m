function [dat] = mygui_defaultoptions(cfg,data)


dat = [];

if ~isfield(cfg,'raw')
    cfg.raw = [];
end

if ~isfield(cfg,'ica')
    cfg.ica = [];
end

if ~isfield(data,'cfg')
    data.cfg = [];
end

if isfield(cfg,'fileout')
    dat.fileout = cfg.fileout;
elseif isfield(data.cfg,'fileout')
    dat.fileout = data.cfg.fileout; 
else
    dat.fileout = 'mypreprocdata';
end

if isfield(data.cfg,'comments')
    dat.comments = data.cfg.comments;
end

%%  1 - EEG channel selection

if ~isfield(data.cfg,'raw')
    data.cfg.raw = [];
end

dat.raw.alllabels = data.label;
if isfield(data,'chaninfo') && isfield(data.chaninfo,'chantype')
    data.cfg.chaninfo = data.chaninfo;
end

% See what channel types are present - fill dat.raw.ctypelabels
if isfield(data.cfg.raw,'ctypelabels')
    dat.raw.ctypelabels = data.cfg.raw.ctypelabels;
elseif isfield(data.cfg,'chaninfo') && isfield(data.cfg.chaninfo,'chantype')
    if length(data.cfg.chaninfo.chantype) ~= length(data.label)
        error('data.(cfg).chaninfo.chantype should have the same length as data.label');
    elseif any(cellfun(@(x) isempty(x),data.cfg.chaninfo.chantype))
        error('data.(cfg).chaninfo.chantype cannot have empty elements');
    end
    idpossibleval = false(length(data.label),1); idpossibleval(1)= true;
    
    for ich = 2:length(data.label)
        if ~ismember(data.cfg.chaninfo.chantype{ich},data.cfg.chaninfo.chantype(idpossibleval))
            idpossibleval(ich) = true;
        end
    end
    dat.raw.ctypelabels = data.cfg.chaninfo.chantype(idpossibleval);
    dat.raw.ctypelabels = dat.raw.ctypelabels(:)';
else
    if isfield(cfg.raw,'ctypelabels')
        if ~iscell(cfg.raw.ctypelabels) || length(cfg.raw.ctypelabels)<1
            fprintf('cfg.raw.ctypelabels is not a cell of length>=1, using default raw channel types\n');
            dat.raw.ctypelabels ={'good','bad','physio'};
        else
            dat.raw.ctypelabels = cfg.raw.ctypelabels;
        end
    else
        dat.raw.ctypelabels ={'good','bad','physio'};
    end
end

% define which one is the default type to be included in main window
if isfield(data.cfg.raw,'ctypemainlabel')
    dat.raw.ctypemainlabel = data.cfg.raw.ctypemainlabel;
elseif isfield(cfg.raw,'ctypemainlabel')
    dat.raw.ctypemainlabel = cfg.raw.ctypemainlabel;
elseif any(strcmpi(dat.raw.ctypelabels,'good'))
    dat.raw.ctypemainlabel = 'good';
else
    dat.raw.ctypemainlabel = dat.raw.ctypelabels{1};
end
dat.raw.ctypemainid = find(strcmpi(dat.raw.ctypelabels,dat.raw.ctypemainlabel));

% define which channel types (if any) are type physio, to be plotted in upper window
if isfield(data.cfg.raw,'ctypephysiolabel')
    dat.raw.ctypephysiolabel = data.cfg.raw.ctypephysiolabel;
elseif isfield(cfg.raw,'ctypephysiolabel')
    dat.raw.ctypephysiolabel = cfg.raw.ctypephysiolabel;
elseif any(strcmpi(dat.raw.ctypelabels,'physio'))
    dat.raw.ctypephysiolabel = 'physio';
else
    dat.raw.ctypephysiolabel = [];
end
dat.raw.ctypephysioid = find(strcmpi(dat.raw.ctypelabels,dat.raw.ctypephysiolabel));


% classify data channels into the categories
if isfield(data.cfg,'chaninfo') && isfield(data.cfg.chaninfo,'chantype')
    for itype = 1:length(dat.raw.ctypelabels)
        num = find(ismember(data.cfg.chaninfo.chantype,dat.raw.ctypelabels{itype}));
        dat.raw.ctype.(dat.raw.ctypelabels{itype}).num = num;
        dat.raw.ctype.(dat.raw.ctypelabels{itype}).label = data.label(num);
    end
else
    if ~isfield(cfg.raw,'ctypepatterns')
        cfg.raw.ctypepatterns = [];
    end
    ctypeaux = zeros(1,length(data.label));
    for itype = setdiff(1:length(dat.raw.ctypelabels),dat.raw.ctypemainid)
        if isfield(data.cfg.raw,'ctype') && isfield(data.cfg.raw.ctype,dat.raw.ctypelabels{itype})
            ctypeaux(ismember(data.label, data.cfg.raw.ctype.(dat.raw.ctypelabels{itype}).label ) ) = itype;
        elseif isfield(cfg.raw.ctypepatterns, dat.raw.ctypelabels{itype}) && ~isempty(cfg.raw.ctypepatterns.(dat.raw.ctypelabels{itype}))
            ctypeaux(ismember(data.label, cfg.raw.ctypepatterns.(dat.raw.ctypelabels{itype}) ) ) = itype;
        elseif strcmpi(dat.raw.ctypelabels{itype},'eog')
            ctypeaux( cellfun(@(x) ~isempty(strfind(lower(x),eog)),data.label) ) = itype;
        end
        num = find(ctypeaux==itype);
        dat.raw.ctype.(dat.raw.ctypelabels{itype}).num = num;
        dat.raw.ctype.(dat.raw.ctypelabels{itype}).label = data.label(num);
    end
    num = find(ctypeaux==0);
    dat.raw.ctype.(dat.raw.ctypemainlabel).num = num;
    dat.raw.ctype.(dat.raw.ctypemainlabel).label = data.label(num);
end


% See if plot bad channels
if isfield(data.cfg.raw,'badchplot') && ~isfield(cfg.raw,'badchplot')
   cfg.raw.badchplot = data.cfg.raw.badchplot; 
end
if isfield(cfg.raw,'badchplot')
    if ~islogical(cfg.raw.badchplot)
        dat.raw.badchplot = false;
        fprintf('cfg.raw.badchplot is not a logical - using default value (false)\n');
    else
        dat.raw.badchplot = cfg.raw.badchplot;
    end
else
    dat.raw.badchplot = false;
end
if dat.raw.badchplot
    if isfield(cfg.raw,'badchcolor') && isnumeric(cfg.raw.badchcolor) ...
            && length(cfg.raw.badchcolor)==3
        dat.raw.badchcolor = cfg.raw.badchcolor;
        dat.raw.badchcolor(dat.raw.badchcolor>1 | dat.raw.badchcolor<0 ) = 0;
    else
        dat.raw.badchcolor = 0.9*[1 1 1];
    end
end

% determine default number of channels to be displayed simultaneously
if isfield(cfg.raw,'nsigblock');
    dat.raw.nsigblock = min( max(round(cfg.raw.nsigblock),1), length(data.label));
    if isnan(dat.raw.nsigblock), dat.raw.nsigblock = length(data.label); end
elseif isfield(data.cfg.raw,'nsigblock');
    dat.raw.nsigblock = data.cfg.raw.nsigblock;
else
    dat.raw.nsigblock = length(data.label);
end

% Plot online average referencing?
if isfield(cfg.raw,'plotavgref') && islogical(cfg.raw.plotavgref)
    dat.raw.plotavgref = cfg.raw.plotavgref;
elseif isfield(data.cfg.raw,'plotavgref')
    dat.raw.plotavgref = data.cfg.raw.plotavgref;
else
    dat.raw.plotavgref = false;
end

% Plot clinical 'banana' montage?
dat.raw.plotbanana = false; 

% Channel units
if isfield(data.cfg,'chaninfo') && isfield(data.cfg.chaninfo,'chanunit')
    dat.raw.unit = data.cfg.chaninfo.chanunit{dat.raw.ctype.(dat.raw.ctypemainlabel).num(1)};
elseif isfield(cfg.raw,'unit')
    dat.raw.unit = cfg.raw.unit;
elseif isfield(data.cfg.raw,'unit')
    dat.raw.unit = data.cfg.raw.unit;
else
    dat.raw.unit = 'uV'; % EEG: uV, V; MEG: T, fT
end

% Default amplitude range to plot 
if isfield(cfg.raw,'amplitude_range_default')
    dat.raw.amplitude_range_default = cfg.raw.amplitude_range_default;
elseif isfield(data.cfg.raw,'amplitude_range_default')
    dat.raw.amplitude_range_default = data.cfg.raw.amplitude_range_default;
else
    dat.raw.amplitude_range_default = 50; % in units of cfg.raw.unit
end
dat.raw.amplitude_range = dat.raw.amplitude_range_default;

% Get channel positions
if isfield(data.cfg.raw,'allpos')
    cfg.raw.allpos = data.cfg.raw.allpos;
end
if isfield(cfg.raw,'allpos')
    if size(cfg.raw.allpos,1)~=length(dat.raw.alllabels) || size(cfg.raw.allpos,2) ~= 2
        error('cfg.raw.allpos must be Nchannels x 2 ')
    else
        dat.raw.allpos = cfg.raw.allpos;
        dat.raw.haspos = ~isnan(cfg.raw.allpos(:,1));
    end
else
    if isfield(cfg.raw, 'layout')
        if exist(cfg.raw.layout,'file')==2
            dat.raw.layout = cfg.raw.layout;
        else
            error('cannot find the requested layout (%s)',cfg.raw.layout)
        end
    else
        dat.raw.layout = 'elec1005.lay';
    end
    fprintf('Reading layout from %s \n',dat.raw.layout);
    if strcmpi(dat.raw.layout(end-3:end),'.mat')
        lay = load(dat.raw.layout);
        aux = fieldnames(lay);
        if any(strcmpi(aux,'lay'))
            lay = lay.lay;
        else
            lay = lay.(aux{1});
        end
        %%% JF DEBUG
        try
            pos = lay.pos;
        catch
            pos = lay.chanpos;
        end
        %%%
        lab = lay.label;
    else
        fileID = fopen(dat.raw.layout,'r');
        dataArray = textscan(fileID, '%*s%f%f%*s%*s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
        fclose(fileID);
        pos = [dataArray{1:2}];
        lab = dataArray{3};
    end
    [found,id] = ismember(lower(data.label),lower(lab));
    dat.raw.allpos = nan(length(dat.raw.alllabels),2);
    dat.raw.allpos(found,:) = pos(id(found),:);
    dat.raw.haspos = found;
end

%% MAIN WINDOW PARAMETERS

if isfield(cfg,'win') && isfield(cfg.win,'tdefault') &&  ...
        isnumeric(cfg.win.tdefault) && cfg.win.tdefault>0 
    dat.win.tdefault = cfg.win.tdefault;
else
    dat.win.tdefault = 10; % [s]
end



%% DATTYPE DEFINITION

if isfield(data.cfg,'dattype')
    dat.dattype.labels = fieldnames(data.cfg.dattype);
    for iaty = 1:length(dat.dattype.labels)
        aux = data.cfg.dattype.(dat.dattype.labels{iaty});
        if isempty(aux)
            dat.dattype.(dat.dattype.labels{iaty}) = zeros(0,2);
        else
            dat.dattype.(dat.dattype.labels{iaty}) = data.cfg.dattype.(dat.dattype.labels{iaty});
        end
    end
    if isfield(cfg,'dattype') && isfield(cfg.dattype,'labels') && iscell(cfg.dattype.labels)
        for iaty = 1:length(cfg.dattype.labels)
            if ~ismember(cfg.dattype.labels{iaty},dat.dattype.labels)
                dat.dattype.(cfg.dattype.labels{iaty}) = zeros(0,2);
                dat.dattype.labels = [dat.dattype.labels(:); cfg.dattype.labels{iaty}];
            end
        end
    end
else
    if isfield(cfg,'dattype') && isfield(cfg.dattype,'labels')
        if ischar(cfg.dattype.labels)
            cfg.dattype.labels = {cfg.dattype.labels};
        elseif iscell(cfg.dattype.labels)
            isok = cellfun(@(x) ischar(x),cfg.dattype.labels);
            cfg.dattype.labels = cfg.dattype.labels(isok);
        end 
        if ~iscell(cfg.dattype.labels) || isempty(cfg.dattype.labels)
            dat.dattype.labels = {'other','bigmusclemov','muscle','ocular'};
        else
            dat.dattype.labels = cfg.dattype.labels;
        end
    else
        dat.dattype.labels = {'other','bigmusclemov','muscle','ocular'};
    end
    for iaty = 1:length(dat.dattype.labels)
        dat.dattype.(dat.dattype.labels{iaty}) = zeros(0,2); %initializes with empty artifacts. 
    end
end
dat.dattype.labels = dat.dattype.labels(:)';
if length(dat.dattype.labels)>9
    dat.dattype.colors = jet(length(dat.dattype.labels));
else
    dat.dattype.colors = [0.9686 0.7608 0.7686; 0.7529 0.7098 0.9647; 0.7373 0.9725 0.6824;0.9725 0.6745 0.4784; 0.9765 0.9176 0.5686; 0.6863 1 1; 1 0.6863 1; 0 1 0.6000];
end
dat.dattype.current = length(dat.dattype.labels);

if isfield(cfg,'dattype') && isfield(cfg.dattype,'send2background')
    if ~iscell(cfg.dattype.send2background)
        fprintf('cfg.dattype.send2background is expected to be a cell - ignoring \ n')
        dat.dattype.background = [];
    else
        send2background = false(1,length(dat.dattype.labels));
        for iaty = 1:length(cfg.dattype.send2background)
            iddattype = find(strcmpi(dat.dattype.labels,cfg.dattype.send2background{iaty}));
            if isempty(iddattype)
                fprintf('Ignoring dattype2background %s - not found in cfg.dattype.labels \n',cfg.dattype.send2background{iaty});
            else
                send2background(iddattype) = true;
            end
        end
        dat.dattype.background = find(send2background);
    end
else
    dat.dattype.background = [];
end



%% ICA PARAMETERS

if isfield(data.cfg,'ica')
    if ~all(isfield(data.cfg.ica,{'a','w'}))
        fprintf('Ignoring data.cfg.ica - a or w matrices not present\n')
        dat.ica = [];
    elseif any(size(data.cfg.ica.a)~=size(data.cfg.ica.w'))
        fprintf('Ignoring data.cfg.ica - a and transpose(w) have different sizes\n')
        dat.ica = [];
    else
        dat.ica = data.cfg.ica;
        if ~isfield(dat.ica,'alllabels')
            for ic = 1:size(dat.ica.a,2)
                dat.ica.alllabels{ic} = sprintf('comp%i',ic);
            end
        end
    end
else
    dat.ica = [];
end

% set number of IC
if isfield(dat.ica,'a')
    dat.ica.numofic = size(dat.ica.a,2);
elseif isfield(cfg.ica,'numofic')
    dat.ica.numofic = cfg.ica.numofic ;
else
    dat.ica.numofic = length(dat.raw.ctype.(dat.raw.ctypemainlabel).num);
end

% Set default data types to exclude from ICA computation
if isfield(dat.ica,'discarddattypeid') 
    if isnumeric(dat.ica.discarddattypeid) && all(dat.ica.discarddattypeid <= length(dat.dattype.labels))
        cfg.ica.discarddattypetype = dat.dattype.labels(dat.ica.discarddattypeid);
    end
end
if isfield(dat.ica,'discarddattype') && ~isfield(cfg.ica,'discarddattype')
   cfg.ica.discarddattype = dat.ica.discarddattype;
end
if isfield(cfg.ica,'discarddattype')
    if strcmpi(cfg.ica.discarddattype,'all')
        cfg.ica.discarddattype = dat.dattype.labels;
    elseif ischar(cfg.ica.discarddattype)
        cfg.ica.discarddattype = {cfg.ica.discarddattype};
    end
    oktype= ismember(cfg.ica.discarddattype, dat.dattype.labels);
    if ~all(oktype)
        fprintf('Data type in cfg.ica.discarddattype do not correspond to data type labels - using default = all\n');
        dat.ica.discarddattypeid = 1:length(dat.dattype.labels);
    else
        dat.ica.discarddattypeid = find(ismember (dat.dattype.labels,cfg.ica.discarddattype));
    end
else
    dat.ica.discarddattypeid = 1:length(dat.dattype.labels);
end

% Number of IC to plot
if isfield(cfg.ica,'nsigblock');
    dat.ica.nsigblock = max(round(cfg.raw.nsigblock),1);
    if isnan(dat.ica.nsigblock)
        dat.ica.nsigblock = dat.ica.numofic;
    end
else
    dat.ica.nsigblock = dat.ica.numofic;
end

% IC labels
if all(isfield(dat.ica,{'ctype','ctypelabels','ctypemainlabel','ctypemainid'})) %if ICa already computed and classification already present
    if isfield(cfg.ica,'ctypelabels') && iscell(cfg.ica.ctypelabels) %if another type is requested in the cfg, add it with 0 ICs
        for ity = 1:length(cfg.ica.ctypelabels)
            if ~ismember(data.cfg.ica.ctypelabels, cfg.ica.ctypelabels{ity})
                dat.ica.ctypelabels{end+1} = cfg.ica.ctypelabels{ity};
                if isfield(cfg.ica,'alllabels') % ica has already been computed
                    dat.ica.ctype.(cfg.ica.ctypelabels{ity}).num = [];
                    dat.ica.ctype.(cfg.ica.ctypelabels{ity}).label = cell(1,0);
                    dat.ica.ctypelabels{end+1} = cfg.ica.ctypelabels{ity};
                end
            end
        end
    end
else
    if all(isfield(dat.ica,{'a','ctype'}))
        dat.ica.ctypelabels = fieldnames(dat.ica.ctype);
    elseif isfield(cfg.ica,'ctypelabels')
        if ~iscell(cfg.ica.ctypelabels) || length(cfg.ica.ctypelabels)<1
            fprintf('cfg.ica.ctypelabels is not a cell of length>=1, using default raw channel types\n');
            dat.ica.ctypelabels ={'good','bad'};
        else
            dat.ica.ctypelabels = cfg.ica.ctypelabels;
        end
    else
        dat.ica.ctypelabels ={'good','bad'};
    end
    % define which one is the default type to be included in main window
    if isfield(cfg.ica,'ctypemainlabel')
        dat.ica.ctypemainlabel = cfg.ica.ctypemainlabel;
    elseif any(strcmpi(dat.ica.ctypelabels,'good'))
        dat.ica.ctypemainlabel = 'good';
    else
        dat.ica.ctypemainlabel = dat.ica.ctypelabels{1};
    end
    dat.ica.ctypemainid = find(strcmpi(dat.ica.ctypelabels,dat.ica.ctypemainlabel));
end

% Add other plotting options that are required if ICA is already computed
if isfield(dat.ica,'a') && ~isfield(dat.ica,'userawch')
    
    if size(dat.ica.a,1) == length(dat.raw.ctype.(dat.raw.ctypemainlabel).num)
        dat.ica.userawch = dat.raw.ctype.(dat.raw.ctypemainlabel).num;
    elseif size(dat.ica.a,1) == length(dat.raw.alllabels)
        dat.ica.userawch = 1:length(dat.raw.alllabels);
    else
        fprintf('Discarding previously computed ICA - cannot find with raw channels were employed to compute it\n')
        dat.ica = rmfield(dat.ica,{'a','w'});
    end
end
if isfield(dat.ica,'a') 

    if ~isfield(dat.ica,'mainplot')
        dat.ica.mainplot.allchlabel = dat.ica.alllabels;
        dat.ica.mainplot.allchid  = 1:1:dat.ica.numofic; %use all components in the main plot
        dat.ica.mainplot.allchcolor = lines(length(dat.ica.mainplot.allchid));
        idbad = setdiff(1:dat.ica.numofic,dat.ica.ctype.(dat.ica.ctypemainlabel).num);
        dat.ica.mainplot.allchcolor(idbad,:) = repmat(0.6*[1 1 1],length(idbad),1);
        dat.ica.mainplot.chhighlight = [];
        dat.ica.mainplot.showchid = dat.ica.nsigblock : -1:1;
    end
    
    if ~isfield(dat.ica,'amplitude_range')
        dat.ica.amplitude_range = 0.5*max(max(dat.ica.w  *  data.trial{1}(dat.ica.userawch,:) ,[], 2));
        dat.ica.amplitude_range_default = dat.ica.amplitude_range;
    end
    
    if ~isfield(dat.ica,'ctype')
        for ity = 1:length(dat.ica.ctypelabels)
            dat.ica.ctype.(dat.ica.ctypelabels{ity}).label = cell(0,1);
            dat.ica.ctype.(dat.ica.ctypelabels{ity}).num = nan(0,1);
        end
        dat.ica.ctype.(dat.ica.ctypemainlabel).label = dat.ica.alllabels;
        dat.ica.ctype.(dat.ica.ctypemainlabel).num = 1:1:dat.ica.numofic;
    else
        for ity = 1:length(dat.ica.ctypelabels)
            % check if that label type is even stored
            if ~isfield(dat.ica.ctype,dat.ica.ctypelabels{ity})
                dat.ica.ctype.(dat.ica.ctypelabels{ity}).num = []; % create subfield
                dat.ica.ctype.(dat.ica.ctypelabels{ity}).label = dat.ica.alllabels(dat.ica.ctype.(dat.ica.ctypelabels{ity}).num);
            elseif ~isfield(dat.ica.ctype.(dat.ica.ctypelabels{ity}),'label')
                dat.ica.ctype.(dat.ica.ctypelabels{ity}).label = dat.ica.alllabels(dat.ica.ctype.(dat.ica.ctypelabels{ity}).num);
            end
        end
    end
    
    
end





%% TRIAL INFORMATION
if isfield(data.cfg,'trl')
    dat.trl.trl = data.cfg.trl.trl;
    if isfield(data.cfg.trl,'nextra')
        dat.trl.nextra = data.cfg.trl.nextra;
        trlaux = data.cfg.trl.trl;
        trlaux(:,1) = trlaux(:,1) + data.cfg.trl.nextra;
        trlaux(:,2) = trlaux(:,2) - data.cfg.trl.nextra;
        if size(trlaux,2)>2
            trlaux(:,3) = trlaux(:,3) + data.cfg.trl.nextra; % or not have 3 columns if sampleinfo
        end
        dat.trl.trl = trlaux; % This trl is just used for plotting / computing ICA. This way the padding data will not be included
    end
    if size(dat.trl.trl,1)~=length(data.trial)
        mystring = sprintf('using data.cfg.trl.trl with %i trials - but %i trials present in data.trial',size(data.cfg.trl.trl,1),length(data.trial));
        error(mystring);
    end
else
    dat.trl.trl = [data.sampleinfo, zeros(size(data.sampleinfo,1),1)]; 
end
dat.trl.n_trl    = size(dat.trl.trl,1);
if ~isfield(dat.trl,'nextra')
    dat.trl.nextra = 0;
end
%% EVENT INFORMATION

if isfield(data.cfg,'event') && ( ~isfield(cfg,'plotevents') || cfg.plotevents )
    if ~all(isfield(data.cfg.event,{'value','sample'})) || ~isstruct(data.cfg.event)
        fprintf('Ignoring data.cfg.event - fields value or sample not found\n')
    elseif ~isnumeric(data.cfg.event(1).sample)
        fprintf('Ignoring data.cfg.event - sample field is not numeric\n')
    else
        auxsample = nan(length(data.cfg.event),1);
        auxval = cell(size(auxsample));
        for iev = 1: length(data.cfg.event)
            if ~isempty(data.cfg.event(iev).value)
                auxsample(iev) = data.cfg.event(iev).sample;
                if isnumeric(data.cfg.event(iev).value)
                   auxval{iev} = num2str(data.cfg.event(iev).value); 
                else
                auxval{iev} = data.cfg.event(iev).value;
                end
            end
        end
        idok = ~isnan(auxsample);
        
        if any(idok)
            dat.event.sample = auxsample(idok); dat.event.value = auxval(idok);
            idpossibleval = false(sum(idok),1); idpossibleval(1)= true;
            for iev = 2:length(idpossibleval)
                if ~ismember(dat.event.value{iev},dat.event.value(idpossibleval))
                   idpossibleval(iev) = true; 
                end
            end
            possibleval = dat.event.value(idpossibleval);
            whatval = nan(size(idpossibleval));
            for ival = 1:length(possibleval)
               whatval(ismember(dat.event.value,possibleval{ival})) = ival;
            end
            bascolor = jet(length(possibleval));
            dat.event.color = bascolor(whatval,:);
        end
        
    end
end
%% DISPLAY OPTIONS

dat.disp.alllabel = {'raw','ica','cleanraw','badproj'};
dat.disp.label = 'raw';
dat.disp.icur = find(strcmpi(dat.disp.alllabel,dat.disp.label));
dat.disp.disptype = 'raw';


%% VISUALIZATION FILTERING PARAMETERS

if ~isfield(cfg.raw,'visfilt')
    dat.raw.visfilt = [];
else
    if isfield(dat.raw.visfilt,'lpfreq')
        dat.raw.visfilt.lpfreq = cfg.raw.visfilt.lpfreq;
        dat.raw.visfilt.filter = true;
    end
    if isfield(dat.raw.visfilt,'hpfreq')
        dat.raw.visfilt.hpfreq = cfg.raw.visfilt.hpfreq;
        dat.raw.visfilt.filter = true;
    end
end

%% ADD OTHER NECESSARY PARAMETERS
dat.fsample = data.fsample;
dat.trl.last    = 1; % last trial
dat.trl.current = 1; % current trial
dat.stop     = 0;
dat.changes  = 0; %if there has been changes in artifact information and something has to be saved
dat.back     = 0;

dat.win.t = dat.win.tdefault;
dat.win.n = dat.win.t*dat.fsample;
dat.win.begs = 1;
dat.win.ends = min([diff(dat.trl.trl(dat.trl.current,1:2))+1, dat.win.n]);




%% MAIN PLOT 
% EEG channels
auxmainlabel = dat.raw.ctypemainlabel;
dat.raw.mainplot.allchid = dat.raw.ctype.(auxmainlabel).num(:);
dat.raw.mainplot.allchcolor = lines(length(dat.raw.mainplot.allchid));
if dat.raw.badchplot 
    auxnum = dat.raw.ctype.bad.num(:);
    dat.raw.mainplot.allchid = [dat.raw.mainplot.allchid; auxnum];
    dat.raw.mainplot.allchcolor = [dat.raw.mainplot.allchcolor; dat.raw.badchcolor(ones(1,length(auxnum)),:)];
end
[dat.raw.mainplot.allchid, goodorder] = sort(dat.raw.mainplot.allchid,'ascend');
dat.raw.mainplot.allchcolor = dat.raw.mainplot.allchcolor(goodorder,:);
dat.raw.mainplot.allchlabel = dat.raw.alllabels(dat.raw.mainplot.allchid);

dat.raw.mainplot.showchid = min(dat.raw.nsigblock, length(dat.raw.mainplot.allchid)):-1:1;
dat.raw.nsigblock = length(dat.raw.mainplot.showchid);
dat.raw.mainplot.chhighlight = [];








