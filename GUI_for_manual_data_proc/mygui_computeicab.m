function mygui_computeicab(h,~,data)
cfg=guidata(h);

display(' ')
display('this is the edited function')
display(' ')


% Computes ICA from the raw data with fastica
% Parameters:
% - number of IC wanted (cfg.ica.numofic, can be modified interactively)
% - artifacts types that will not be included for ICA computation
%       (cfg.ica.discarddattypeid, can be modified interactively)
% - raw channels that are included in the ICA: the raw channels having the
%       same label as cfg.raw.ctypemainlabel (i.e. 'good')
% Output: cfg variable passed in guidata, with fields .ica. ...
%   numofic  = number of IC
%   a = mixing matrix
%   w = separating maxtrix
%   discarddattypeid  = artifact samples not included in ICA
%   userawch = channels indices from all data that were used in the ICA


%% 1 - Get Current parameters
param = [];
goodlabel = cfg.raw.ctypemainlabel;
param.nrawch = length(cfg.raw.ctype.(goodlabel).label);

if isfield(cfg.ica,'numofic')
    tgnum = round(cfg.ica.numofic);
    if tgnum<0 || tgnum> param.nrawch || isnan(tgnum)
        cfg.ica.numofic = param.nrawch;
        fprintf('cannot use specified number of components - using %i ICs \n',cfg.ica.numofic);
    end
else
    cfg.ica.numofic = param.nrawch;
end
param.numofic = cfg.ica.numofic;
param.discarddattypeid = cfg.ica.discarddattypeid;

%% 2- Create popup window

pos      = get(0,'DefaultFigurePosition');
pos(3:4) = [300 400];
% dlg      = dialog('Name', 'ICA computation', 'Position', pos);
dlg = figure('Position', pos);
set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears

%% 3 - Create uicontrol elements in the window
uicontrol(dlg, 'style', 'pushbutton', 'position', [ 55  10    80  20], 'string', 'Cancel',       'callback', 'close');
uicontrol(dlg, 'style', 'pushbutton', 'position', [155  10    80  20], 'string', 'OK',           'callback', 'uiresume');

uicontrol(dlg,'units','pixels','position', [20 320 100 30],'Style', 'text', 'String', 'number of ICs');
param.boxnumofic = uicontrol(dlg,'units','pixels','position',[20 290 100 30], ...
    'Style','edit','Tag','icnumber','String',num2str(cfg.ica.numofic),'Callback',@icnumber,'HorizontalAlignment','center');

uicontrol(dlg,'units','pixels','position', [150 320 100 30],'Style', 'text', 'String', 'Discard samples of type');
param.boxusdattypediscard = uicontrol(dlg,'units','pixels','position',[150 220 100 100],'min',0,'max',2, ...
    'Style','listbox','String',cfg.dattype.labels,'Value',param.discarddattypeid,'Callback',@artdiscard,'HorizontalAlignment','center');


% Start with previous answer.
if isfield(cfg.ica,'a')
    oldA = cfg.ica.a;
    oldW = cfg.ica.w;
    if size(cfg.ica.a,1)~= length(cfg.raw.ctype.(cfg.raw.ctypemainlabel).num)
        fprintf('Cannot use the same components as initial guesses because the employed channels have changed\n');
        param.boxreuse.Value = 0;
    else
        param.textreuse = uicontrol(dlg,'units','pixels','position', [50 120 200 50],'Style', 'text', 'String', 'Use old ICA as initial guess');
        param.boxreuse = uicontrol(dlg,'units','pixels','position',[30 150 20 20], ...
            'Style','checkbox','Value',0,'HorizontalAlignment','left'); %,'Callback',@artdiscard);
    end
    oldorganization = nan(size(oldW,1),1);
    for ity = 1:length(cfg.ica.ctypelabels)
        oldorganization(cfg.ica.ctype.(cfg.ica.ctypelabels{ity}).num) = ity;
    end
else
    oldA = [];
    param.boxreuse.Value = 0;
end

set(dlg,'KeyPressFcn',@keypress);
%% loop until user presses cancel, ok, or closes

dlg.UserData = param;
uiwait(dlg);
if ishandle(dlg)
    param = dlg.UserData;
    
    
    %% Get the values of the final parameters
    %(in case modified interactively)
    cfg.ica.numofic = param.numofic;
    cfg.ica.discarddattypeid = param.discarddattypeid;
    fprintf('Computing %i ICs from %i EEG channels\n',cfg.ica.numofic,param.nrawch)
    
    boxstatus = uicontrol(dlg,'units','pixels','position', [80 50 150 60],'Style', 'text', ...
        'String', sprintf('Computing %i ICs from %i EEG channels\n',cfg.ica.numofic,param.nrawch));
    drawnow
    
    %% Compute ICA
    
    % Channels and samples to include
    idch = cfg.raw.ctype.(cfg.raw.ctypemainlabel).num;
    usesample = true(1,max(cfg.trl.trl(:,2)));
    for iaty = cfg.ica.discarddattypeid
        aux =  cfg.dattype.(cfg.dattype.labels{iaty});
        for iart = 1:size(aux,1)
            usesample(aux(iart,1):aux(iart,2)) = false;
        end
    end
    
    % Concatenate all wanted data
    alldat = cell(1,length(data.trial));
    alldateog = cell(1,length(data.trial));
    icheogv = find(strcmpi(data.label,'EOGV'));
    for itrial = 1:length(data.trial)
        usesampletrial = usesample(cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2));
        if isfield(cfg.trl,'nextra')
            aux = double(data.trial{itrial}(idch,1+cfg.trl.nextra:end-cfg.trl.nextra));
            auxeog = double(data.trial{itrial}(icheogv,1+cfg.trl.nextra:end-cfg.trl.nextra));
        else
            aux = double(data.trial{itrial}(idch,:));
            auxeog = double(data.trial{itrial}(icheogv,:));
        end
        alldat{itrial} = aux(:,usesampletrial);
        alldateog{itrial} = auxeog(:,usesampletrial);
    end
    alldat = cat(2,alldat{:});
    alldateog = cat(2,alldateog{:});
    
    alldat = cat(1, alldat,alldateog);
    
    % Compute ICA
    
    if param.boxreuse.Value == 0
        [A,W] = fastica(alldat,...
            'approach','defl',...
            'g','pow3',... % ,'tanh'
            'numOfIC',cfg.ica.numofic,...
            'stabilization' ,'on',...
            'verbose','on',... %'sampleSize',0.25,...
            'displayMode','off'); %,'interactivepca','on');
    else
        if cfg.ica.numofic > size(oldA,2)
            initguess = [oldA rand(size(oldA,1), cfg.ica.numofic-size(oldA,2) )];
        else
            initguess = oldA(:,1:cfg.ica.numofic);
        end
        [A,W] = fastica(alldat,...
            'approach','defl',...
            'g','pow3',... % ,'tanh'
            'numOfIC',cfg.ica.numofic,...
            'stabilization' ,'on',...
            'verbose','on',... %'sampleSize',0.25,...
            'displayMode','off', ...
            'initGuess', initguess);
    end
    A = A(1:end-length(icheogv),:);
    W = W(:,1:end-length(icheogv));
    alldat = alldat(1:end-length(icheogv),:);
    % sort by power
    [~,index] = sort(-sum(A.^2));
    
    
    
    if size(W,1) ~= cfg.ica.numofic
        cfg.ica.numofic = size(W,1);
        boxstatus.String =  sprintf('Only %i ICs were extracted\n',cfg.ica.numofic);
        drawnow
    end
       
    
    %% Save results
    
    cfg.ica.userawch = idch;
    cfg.ica.alllabels = cell(cfg.ica.numofic,1);
    for ic = 1:cfg.ica.numofic
        cfg.ica.alllabels{ic} = sprintf('comp%i',ic);
    end
    
    
    
    % Update mainplot information
    cfg.ica.mainplot.allchlabel = cfg.ica.alllabels;
    cfg.ica.mainplot.allchid  = 1:1:cfg.ica.numofic; %use all components in the main plot
    cfg.ica.mainplot.allchcolor = lines(length(cfg.ica.mainplot.allchid));
    cfg.ica.mainplot.chhighlight = [];
    cfg.ica.nsigblock = min(cfg.ica.nsigblock,cfg.ica.numofic);
    cfg.ica.mainplot.showchid = cfg.ica.nsigblock : -1:1;
    
    cfg.ica.amplitude_range = 0.5*max(max(W  *  alldat(:,1:min(10000,size(alldat,2))) ,[], 2));
    cfg.ica.amplitude_range_default = cfg.ica.amplitude_range;

    
    
    % change display to ICA components
    cfg.disp.label = 'ica';
    cfg.disp.icur = find(strcmpi(cfg.disp.alllabel,'ica'));
    cfg.disp.disptype = 'ica';
    cfg.disp.guichanselect.String = 'components';
    cfg.win.mainplotbegch.String = num2str(min(cfg.ica.mainplot.showchid));
    cfg.win.mainplotendch.String = num2str(max(cfg.ica.mainplot.showchid));
    cfg.disp.guiselect.Value = cfg.disp.icur;
    
    
    cfg.ica.a=A(:,index);
    cfg.ica.w=W(index,:);
    
    
    % Erase component classfication
    for ity = 1:length(cfg.ica.ctypelabels)
        cfg.ica.ctype.(cfg.ica.ctypelabels{ity}).label = cell(0,1);
        cfg.ica.ctype.(cfg.ica.ctypelabels{ity}).num = nan(0,1);
    end
    cfg.ica.ctype.(cfg.ica.ctypemainlabel).label = cfg.ica.alllabels;
    cfg.ica.ctype.(cfg.ica.ctypemainlabel).num = 1:1:cfg.ica.numofic;
    
    %% Explore relation of ICa components with artifacts
    
    % time series of presence of artifacts
    sampletype = zeros(1,max(cfg.trl.trl(:)));
    for iaty = 1:length(cfg.dattype.labels)
        aux =  cfg.dattype.(cfg.dattype.labels{iaty});
        for iart = 1:size(aux,1)
            sampletype(aux(iart,1):aux(iart,2)) = iaty;
        end
    end
    
    alldattype = cell(1,length(data.trial));
    physioch = cell(1,length(data.trial));
    idchphysio = cfg.raw.ctype.(cfg.raw.ctypephysiolabel).num;
    for itrial = 1:length(data.trial)
        idx = cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2);
        usesampletrial = usesample(cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2));
        alldattype{itrial} = sampletype(idx(usesampletrial));
        if isfield(cfg.trl,'nextra')
            aux = data.trial{itrial}(idchphysio,1+cfg.trl.nextra:end-cfg.trl.nextra);
        else
            aux = data.trial{itrial}(idchphysio,:);
        end
        physioch{itrial} = aux(:,usesampletrial);
    end
    alldattype = cat(2,alldattype{:});
    physioch = cat(2,physioch{:});
    
    icatimeseries = cfg.ica.w*alldat;
    for iaty = setdiff(1:length(cfg.dattype.labels),cfg.ica.discarddattypeid)
        cfg.ica.artifactpower.(cfg.dattype.labels{iaty}) = ...
            mean(abs(icatimeseries(:,alldattype==iaty)).^2,2) ./ ...
            mean(abs(icatimeseries(:,alldattype==0)).^2,2);
    end
    
    if cfg.ica.numofic>=3
        for ich = 1:length(idchphysio)
            mycorr = corr(physioch(ich,:)',icatimeseries');
            [~,idx] = sort(mycorr,'descend');
            fprintf('Highest correlation with %s - %.1f (comp%i) - %.1f (comp%i) - %.1f (comp%i)\n',...
                cfg.raw.ctype.(cfg.raw.ctypephysiolabel).label{ich},mycorr(idx(1)),idx(1), ...
                mycorr(idx(2)),idx(2),mycorr(idx(3)),idx(3));
        end
    end
    

    
    %% Reuse component classification from previous ICA computation
    if param.boxreuse.Value ~= 0
        
        % id interest
        id = find(oldorganization ~= cfg.ica.ctypemainid);
        if ~isempty(id)
            corroldnew = abs(corr((oldW(id,:)*alldat)',icatimeseries'));
            [maxcor,idmaxcor] = max(corroldnew,[],2);
            reuseclassif = maxcor>0.9;
            auxlabels = cfg.ica.ctypelabels(oldorganization(id));
            d = [num2cell(id(:)), auxlabels(:), ...
                num2cell(idmaxcor),num2cell(maxcor),num2cell(reuseclassif)];
            % figclas  = dialog('Name', 'Reuse component classification'); %, 'Position', pos);
            figclas = figure;
            t = uitable('Data', d,'Units','normalized','Position',[0.05 0.25 0.90 0.70],...
                'ColumnName', {'old component','type','new component','correlation','assign label'},...
                'ColumnFormat', {'numeric','char','numeric','numeric','logical'},...
                'ColumnEditable', [false false false false true],...
                'RowName',[]);
            uicontrol(figclas, 'style', 'pushbutton', 'units','normalized', ...
                'position', [ 0.1  0.1    0.3  0.1], 'string', 'Cancel',       'callback', 'close');
            uicontrol(figclas, 'style', 'pushbutton','units','normalized', ...
                'position', [ 0.5  0.1    0.3  0.1], 'string', 'OK',           'callback', 'uiresume');
            uiwait(figclas);
            if ishandle(figclas)
                neworganization = cfg.ica.ctypemainid * ones(1,cfg.ica.numofic);
                for ic = 1:length(id)
                   if t.Data{ic,5} 
                       neworganization(t.Data{ic,3}) = ...
                           find(strcmpi(cfg.ica.ctypelabels,t.Data{ic,2}));
                   end
                end
                for ity = 1:length(cfg.ica.ctypelabels)
                    aux = find(neworganization == ity);
                    cfg.ica.ctype.(cfg.ica.ctypelabels{ity}).num = aux;
                    cfg.ica.ctype.(cfg.ica.ctypelabels{ity}).label = cfg.ica.alllabels(aux);
                end
                idbad = setdiff(1:cfg.ica.numofic,cfg.ica.ctype.(cfg.ica.ctypemainlabel).num);
                cfg.ica.mainplot.allchcolor(idbad,:) = repmat(0.6*[1 1 1],length(idbad),1);
                
                
            end
            close(figclas)
        end
    end
    
    %%  power in 40-80Hz range
    alldat = nan(length(cfg.ica.userawch),max(cfg.trl.trl(:,2)));
    for itrial = 1:length(data.trial)
        if isfield(cfg.trl,'nextra')
            aux = data.trial{itrial}(idch,1+cfg.trl.nextra:end-cfg.trl.nextra);
        else
            aux = data.trial{itrial}(idch,:);
        end
        alldat(:,cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2)) = aux;
    end
    alldat(:,~usesample) = nan;
    cfg.ica.artifactpower.pow25_80Hz = nan(cfg.ica.numofic,1);
    optfreq = []; optfreq.foilim = [25 80];
    [foi,auxpowspctrm] = my_compute_powerspctrm(optfreq,cfg.fsample,cfg.ica.w*alldat);
    cfg.ica.artifactpower.pow25_80Hz = sum(auxpowspctrm,1)*mean(diff(foi));
    
    
    close(dlg);
end
cfg.redraw.main = true; cfg.redraw.topo = true; cfg.redraw.dattypelines = true;
guidata(h,cfg);
uiresume;
end


function icnumber(dlg, eventdata)
dlg = get(dlg, 'parent');
param = get(dlg, 'userdata');
tgnum = str2double(eventdata.Source.String);
if isempty(eventdata.Source.String)
    param.boxnumofic.String = num2str(param.numofic);
elseif tgnum<0 || tgnum> param.nrawch || isnan(tgnum)
    param.boxnumofic.String = num2str(param.numofic);
    fprintf('cannot use %s number of components - using %i ICs \n',eventdata.Source.String,param.numofic);
else
    param.numofic = tgnum;
end
set(dlg, 'userdata', param);
end

function artdiscard(dlg, eventdata)
dlg = get(dlg, 'parent');
param = get(dlg, 'userdata');
tgval =  eventdata.Source.Value;
if tgval == param.discarddattypeid %if repeated click deselect
    tgval = [];
    param.boxusdattypediscard.Value = [];
end
param.discarddattypeid = tgval;
set(dlg, 'userdata', param);
end


% Allow key - controll
function keypress(h, ~)
key=double(get(gcbf, 'CurrentCharacter'));
if  key
    switch key
        case 13     % arrow left
            uiresume(h);
            
    end
end
end