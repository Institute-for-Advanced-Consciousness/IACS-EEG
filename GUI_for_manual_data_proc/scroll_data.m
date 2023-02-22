
function [] = scroll_data(cfg,data)
%% STARTS GUI TO SCROLL DATA,
% INPUT RAW DATA either as the second argument or in cfg.datafile = full
% path to the data to load.
% - If a cfg.preprocfun exists with a valid handle to a function in the path
%   exists, it will be used to load the data (i.e. if cfg.preprocfun =
%   @mypreprocfun_projectx, then data will be loaded as data = mypreprocfun_projectx(cfg)  )
% - Else cfg.datafile must store a single variable with fieldtrip data
%   structure

% OUTPUT DATA
% If clicking in Done, data will be saved in data.cfg.fileout (if existent),
% or cfg.fileout (default = mypreprocdata.mat)
% The new file may contain modifications in the following fields only:
%   - data.sampleinfo (if .sampleinfo field did not exist before, it is
%   added)

% Overview of functionalities:
% DEFINE DATA TYPES,
% - It loads previous data types definition from data.cfg.dattype if existent
% - Else artifact types can be defined (initialized empty) in cfg.dattypes.labels 
% - Else it creates default data types: {'other','bigmusclemov','muscle','ocular'};
% CHANNEL SELECTION,
% - Loaded from data.cfg.raw if existent
% - Else channel types and patterns can be defined in cfg.raw
%    cfg.raw.ctypelabels = raw channel types (default =  {'good','bad','physio'}
%    cfg.raw.ctypepatterns.typex = labels that will be considered as typex
%           (default = [])
%    cfg.raw.ctypemainlabel = main raw channel type, which will be plotted
%       in the main axes, and will be employed for ICA (default = 'good')
%    cfg.raw.ctypephysiolabel = channel type that will be plotted in the
%       upper plot (default = 'physio')
% ICA COMPUTATION
% - Loaded from data.cfg.ica if existent
% - Else reads from cfg.ica
%       cfg.ica.numofic = number of independent components (default =
%           number of good channels)
%       cfg.ica.discarddattype = data types that should not be included
%           in ICA computation (default = 'all')
%       cfg.ica.ctypelabels = component types (default = {'good','bad'})
%       cfg.ica.ctypemainlabel = component type that is considered as good
%           (not noisy) (default = 'good')
% (Note: ICA will use all raw channels that are labeled as ctypemainlabel)
% TRIAL INFORMATION. is read from data.cfg.trl, which should contain
% either data.cfg.trl.trl with [begsample endsample offset (trialinfo)] for
% each trial, or data.sampleinfo with [begsample endsample] informatoin for
% each trial
% Optional: data.cfg.trl.nextra can contain the amount of samples that are
%   used for padding in the input data but should not be displayed or used
%   in ICA



% other options:
% DISPLAY OPTIONS
% cfg.win.tdefault = data length (seconds) displayed in the main plot
%       (default = 10)
% cfg.raw.nsigblock = maximum number of channels to display simultaneously
%       in the main window (default = plots all channels simultaneously)
% cfg.raw.plotavgref = true if data of channels type ctypemainlabel should
%       be averaged re-referenced for plotting (default = false)
% cfg.raw.badchplot = if true plot badchannels in light gray, if false
%       do not plot them (default = false)
% cfg.raw.unit = units of the raw data (default = 'uV')
% cfg.raw.amplitude_range_default = limits of the verical axes plotted for
%       each channel in the main plot (in units of cfg.raw.unit) (default = 50)
% RAW CHANNEL POSITIONS. Necessary for topoplots. Can be added either with
%      cfg.raw.allpos = Nchannel x 2 matrix with coordinates, or
%      cfg.raw.layout = full path to target layout (default 'elec1005.lay')
%

%% LOAD DATA

dbstop if error

if nargin <2
    if isfield(cfg,'preprocfun')
        f = cfg.preprocfun;
        data = f(cfg);
    else
        data = load(cfg.filename);
        aux = fieldnames(data);
        if length(aux) == 1
            data = data.(aux{1});
        elseif any(strcmpi(aux,'data'))
            data = data.data;
        elseif length(aux)>1
            data = data.(aux{1});
        end
    end
end

% Add sampleinfo if not present
if ~isfield(data,'sampleinfo')
    if isfield(data,'cfg') && isfield(data.cfg,'trl')
        data.sampleinfo = data.cfg.trl.trl(:,1:2);
    else
        fprintf('No trial or sampleinfo present - creating trl info assumming consecutive and adjacent trials\n')
        data.sampleinfo = nan(length(data.trial),2);
        data.sampleinfo(1,:) = [1 length(data.time{1})];
        for itrial = 2:length(data.trial)
            data.sampleinfo(itrial,1) = data.sampleinfo(itrial-1,2) + 1;
            data.sampleinfo(itrial,2) = data.sampleinfo(itrial,1) + length(data.time{itrial})-1;
        end
    end
end

%% GET CONFIGURATION OPTIONS


cfg = mygui_defaultoptions(cfg,data);


%% CREATE VARIABLE WITH DATA TO PLOT
% copy from the original data but will be modified to plot (filtering, ...)
plotdata.trial = cell(size(data.trial));
plotdata.time = data.time;
for itrial = 1:length(data.trial)
    if isfield(cfg.trl,'nextra')
        plotdata.trial{itrial}=double(data.trial{itrial}(:,1+cfg.trl.nextra:end-cfg.trl.nextra));
        plotdata.time{itrial}=data.time{itrial}(1+cfg.trl.nextra:end-cfg.trl.nextra);
    else
        plotdata.trial{itrial}=double(data.trial{itrial});
    end
end



%% CREATE FIGURE AND MAIN AXES


fig = figure('Color','w','Renderer','zbuffer'); % changer render method to 'painters' to restore orginal setting
set(fig,'KeyPressFcn',@keypress);
set(fig,'WindowButtonDownFcn',@mouse_click);
figure(fig)

cfg.figlay.maxmain = [0.05  0.1 0.80  0.7];
cfg.figlay.maxphysio = [0.05 0.80 0.80 0.15];
cfg.figlay.topowidth = 0.05;
aux = cfg.figlay.maxmain; aux(3) = cfg.figlay.topowidth ;
cfg.figlay.topopos = aux;
cfg.h_main  = axes('Parent',fig,'Position',cfg.figlay.maxmain); hold on
h_selectionlines = plot([],[]);
cfg.h_physio  = axes('Parent',fig,'Position',cfg.figlay.maxphysio); box on; set(cfg.h_physio,'XTick',[]);
cfg.h_topo = axes('Parent',fig,'Position',cfg.figlay.topopos); axis off; hold on;


axes('Position',[0.05  0.965 0.6  0.04]); axis off
cfg.h_title = text(0,0,'Loading ...','Interpreter','none');
drawnow

colormap(fig,'jet')
%% CREATE UI BUTTONS

% done
uicontrol(fig,'units','normalized','position',[0.005 0.005  0.03  0.03],'String','done','Callback',@stop);

% filtering
uicontrol(fig,'units','normalized','position',[0.06 0.005  0.06  0.03],'String','filter','Callback',@mygui_definefiltering);

% toggle montage 
uicontrol(fig,'units','normalized','position',[0.14 0.005  0.06  0.03],'String','clin. montage','Callback',@mygui_togglemontage);

% change to previous/next trial
uicontrol(fig,'units','normalized','position',[0.205 0.03 0.055 0.03],'Style', ...
    'text','String', 'trial','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
cfg.trl.currentbox = uicontrol(fig,'units','normalized','position',[0.205 0.005 0.055 0.03], ...
    'Style','edit','String',num2str(cfg.trl.current),'Callback',@trialnumber,'HorizontalAlignment','center');
uicontrol(fig,'units','normalized','position',[0.18  0.005  0.025  0.03],'String','<','Callback',@prevtrial);
uicontrol(fig,'units','normalized','position',[0.26  0.005  0.025  0.03],'String','>','Callback',@nexttrial);


% change to previous/next window
uicontrol(fig,'units','normalized','position',[0.33  0.0305  0.17  0.03],'Style', ...
    'text','String', 'time','BackgroundColor',[1 1 1],'HorizontalAlignment','center')
uicontrol(fig,'units','normalized','position',[0.305  0.005  0.025  0.03],'String','<','Callback',@prevwindow);
uicontrol(fig,'units','normalized','position',[0.50  0.005  0.025  0.03],'String','>','Callback',@nextwindow);
cfg.h_bar   = axes('Parent',fig,'Position',[0.3295  0.01  0.17  0.03]); axis off; box on; hold on;


% window length
uicontrol(fig,'units','normalized','position',[0.54 0.03 0.055 0.03],'Style', ...
    'text','String', 'window length','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
cfg.win.tbox = uicontrol(fig,'units','normalized','position',[0.54 0.005 0.055 0.03], ...
    'Style','edit','String',num2str(cfg.win.t),'Callback',@changewindowt,'HorizontalAlignment','center');

% amplitude
uicontrol(fig,'units','normalized','position',[0.64 0.03 0.055 0.03],'Style', ...
    'text','String', 'vert scale','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
cfg.win.verts = uicontrol(fig,'units','normalized','position',[0.64 0.005 0.055 0.03], ...
    'Style','edit','String',num2str(cfg.raw.amplitude_range),'Callback',@changeverts,'HorizontalAlignment','center');
uicontrol(fig,'units','normalized','position',[0.615  0.005  0.025  0.03],'String','<','Callback',@changeverts);
uicontrol(fig,'units','normalized','position',[0.695  0.005  0.025  0.03],'String','>','Callback',@changeverts);


% channelindex
uicontrol(fig,'units','normalized','position',[0.76 0.03 0.07 0.03],'Style', ...
    'text','String', 'displaych','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
cfg.win.mainplotbegch = uicontrol(fig,'units','normalized','position',[0.76 0.005 0.035 0.03], ...
    'Style','edit','String',num2str(min(cfg.raw.mainplot.showchid)),'Tag','begch','Callback',@changechanind,'HorizontalAlignment','center');
cfg.win.mainplotendch = uicontrol(fig,'units','normalized','position',[0.795 0.005 0.035 0.03], ...
    'Style','edit','String',num2str(max(cfg.raw.mainplot.showchid)),'Tag','endch','Callback',@changechanind,'HorizontalAlignment','center');
uicontrol(fig,'units','normalized','position',[0.735  0.005  0.025  0.03],'Tag','<','String','<','Callback',@changechanind);
uicontrol(fig,'units','normalized','position',[0.83  0.005  0.025  0.03],'Tag','>','String','>','Callback',@changechanind);


% % GET RID OF OLD DATA TYPES (just clutters the GUI)
% 
% for icell = 1:length(cfg.dattype.labels)
%     switch cfg.dattype.labels{icell}
%         case 'notched_delta'
%             cfg.dattype.labels{icell} = [];
%         case 'AS_theta'
%             cfg.dattype.labels{icell} = [];
%         case 'sawtooth_delta'
%             cfg.dattype.labels{icell} = [];
%     end
% end
% 
% % remove empty cells
% keep = ~cellfun('isempty', cfg.dattype.labels);
% cfg.dattype.labels = cfg.dattype.labels(keep);

% datatype definition
cfg.dattype.buttontype = cell(1,length(cfg.dattype.labels));
%We have so many labels that we need more colors now 
HSV = hsv; % get the hsv colormap
chrome = HSV(1:50:end,:)./3 + 0.66;
cfg.dattype.colors = chrome(end:-1:1,:);
for iaty = 1:length(cfg.dattype.labels)
    cfg.dattype.buttontype{iaty} = uicontrol(fig,'units','normalized','position', ...
        [0.87  0.4+iaty*0.025  0.13  0.025],'Style', 'togglebutton', 'String', cfg.dattype.labels{iaty},...
        'BackgroundColor',cfg.dattype.colors(iaty,:),...
        'Callback', @dattype_buttonpush);

end
cfg.dattype.buttontype{cfg.dattype.current}.Value = 1;
cfg.dattype.buttontype{cfg.dattype.current}.FontWeight= 'bold';

cfg.dattype.buttonadd = uicontrol(fig,'units','normalized','position', ...
    [0.87  0.39  0.043  0.04],'Style', 'pushbutton', 'String', 'add',...
    'BackgroundColor',0.7*[1 1 1],...
    'Callback', @button_dattype_change);
cfg.dattype.buttondel = uicontrol(fig,'units','normalized','position', ...
    [0.913  0.39  0.043  0.04],'Style', 'pushbutton', 'String', 'del',...
    'BackgroundColor',0.7*[1 1 1],...
    'Callback', @button_dattype_change);
uicontrol(fig,'units','normalized','position', ...
    [0.956  0.39  0.043  0.04],'Style', 'pushbutton', 'String', 'topo',...
    'BackgroundColor',0.7*[1 1 1],'Callback', {@topo_segment,plotdata});



cfg.dattype.boxstart = uicontrol(fig,'Style','edit','units','normalized','position', ...
    [0.87  0.32  0.07  0.05],'String','-','Callback',@start_datasel);
cfg.dattype.boxstop  = uicontrol(fig,'Style','edit','units','normalized','position',...
    [0.94  0.32  0.07  0.05],'String','-','Callback',@stop_datasel);

% EEG channelselection
cfg.disp.guichanselect = uicontrol(fig,'units','normalized','position', ...
    [0.87  0.25  0.13  0.05],'Style', 'pushbutton', 'String', 'channel/comp',...
    'Callback', @channelselect);

% channelplots
uicontrol(fig,'units','normalized','position', ...
    [0.87  0.20  0.13  0.05],'Style', 'pushbutton', 'String', 'c. plots',...
    'Callback', {@mygui_cplots,plotdata});

% computeica
uicontrol(fig,'units','normalized','position', ...
    [0.87  0.15  0.13  0.05],'Style', 'pushbutton', 'String', 'compute ICA',...
    'Callback', {@mygui_computeica,data});

% gamma thresholding
cfg.figlay.boxtopo = uicontrol(fig,'units','normalized','position', ...
    [0.87  0.1  0.13  0.05],'Style', 'pushbutton', 'String', 'set threshold',...
    'Value',0,'Callback',{@mygui_thresh,data});

% ICA overview
cfg.figlay.boxtopo = uicontrol(fig,'units','normalized','position', ...
    [0.87  0.05  0.13  0.05],'Style', 'pushbutton', 'String', 'ICA overview',...
    'Value',0,'Callback',{@mygui_ICA_overview,plotdata});

% showtopos
cfg.figlay.boxtopo = uicontrol(fig,'units','normalized','position', ...
    [0.87  0.00  0.13  0.05],'Style', 'togglebutton', 'String', 'ICA minitopo',...
    'Value',0,'Callback',@ICA_minitopo);

% Display
uicontrol(fig,'units','normalized','position', [0.87  0.75  0.13  0.05],'Style', ...
    'text', 'String', 'Display')
cfg.disp.guiselect = uicontrol(fig,'units','normalized','position', [0.87  0.7  0.13  0.05],'Style', ...
    'popupmenu', 'String', cfg.disp.alllabel,'Value',cfg.disp.icur,...
    'Callback', @selectdisplaymag);

% Comments
uicontrol(fig,'units','normalized','position', ...
    [0.87  0.88  0.13  0.05],'Style', 'pushbutton', 'String', 'Comments',...
    'Callback', @addcomments);
cfg.comments = [];

% save
uicontrol(fig,'units','normalized','position', ...
    [0.87  0.82  0.13  0.05],'Style', 'pushbutton', 'String', 'save',...
    'Callback', {@mysave,data});

% Dataout
uicontrol(fig,'units','normalized','position', ...
    [0.87  0.94  0.13  0.05],'Style', 'pushbutton', 'String', 'Output file',...
    'Callback', @changeoutfile);
cfg.outfiletext = uicontrol(fig,'units','normalized','position', ...
    [0.47  0.95  0.40  0.05],'Style', 'text', 'String', cfg.fileout, ...
    'BackgroundColor',[1 1 1]); %,'HorizontalAlignment','left');
%%  MAIN SCROLL LOOP


cfg.redraw.physio = true;
cfg.redraw.main = true;
cfg.redraw.dattypelines = true;
cfg.redraw.bar = true;
cfg.redraw.topo = true;

while(ishandle(fig))
    
    %% filter EEG if needed
    if ~isempty(cfg.raw.visfilt) && cfg.raw.visfilt.filter
        plotdata.trial = mygui_dofiltering(cfg,data);
        cfg.raw.visfilt.filter = false;
        cfg.redraw.main = true; 
    end
    
    %% Display banana montage
     if cfg.raw.plotbanana
         plotdata.avgref = plotdata.trial;
         cfg.raw.mainplot.oldLabels = cfg.raw.mainplot.allchlabel;
         cfg.raw.oldLabels = cfg.raw.alllabels;
         [plotdata.trial,newLabels] = mygui_banana(cfg,data);
         cfg.raw.alllabels = newLabels;
         cfg.raw.mainplot.allchlabel = newLabels;
         cfg.raw.mainplot.showchid = length(newLabels):-1:1;
%          for ilb = 1:length(newLabels)
%              cfg.raw.alllabels{ilb} = newLabels{ilb};
%              cfg.raw.mainplot.allchlabel{ilb} = newLabels{ilb};
%          end
%          cfg.raw.alllabels(ilb+1:end) = [];
%          cfg.raw.mainplot.allchlabel(ilb+1:end) = [];
         cfg.redraw.main = true;
     else
         if isfield(plotdata,'avgref')
            plotdata.trial = plotdata.avgref; % revert back 
            cfg.raw.mainplot.allchlabel = cfg.raw.mainplot.oldLabels;
            cfg.raw.alllabels = cfg.raw.oldLabels
            cfg.raw.mainplot.allchlabel = cfg.raw.oldLabels;
            cfg.raw.mainplot.showchid = length(cfg.raw.oldLabels):-1:1;
%             for ilb = 1:length(cfg.raw.oldLabels)
%                 cfg.raw.alllabels{ilb} = cfg.raw.oldLabels{ilb};
%                 cfg.raw.mainplot.allchlabel{ilb} = cfg.raw.oldLabels{ilb};
%             end
            cfg.redraw.main = true;
         end
     end
    
    %% Selects time of interest
    
    % 1 - selects time of interest
    idx_current_view = cfg.win.begs:cfg.win.ends; %samples of the current trial to be displayed.
    
    %%% PATCH ADDED BY JF 08.18.19
    % Sometimes idx_current_view is one element too long, this should fix
    if idx_current_view(end) == length(plotdata.time{cfg.trl.current}) + 1
        idx_current_view = idx_current_view(1:end-1);
    end
    %%%
    
    t    = plotdata.time{cfg.trl.current}(idx_current_view);
    cfg.win.tlim = [min(t) min(t)+cfg.win.t];
    cfg.win.tend = t(end);
    
    
    %% plot physio time series in upper plot
    
    if cfg.redraw.physio
        if ~isempty(cfg.raw.ctypephysiolabel) && ~isempty(cfg.raw.ctype.(cfg.raw.ctypephysiolabel).num)
            % Prepare physio time series
            nch = length(cfg.raw.ctype.(cfg.raw.ctypephysiolabel).num);
            physio = nan(nch,length(t));
            auxid = cfg.raw.ctype.(cfg.raw.ctypephysiolabel).num;
            if length(idx_current_view)>10
                for i=1:nch
                    aux = double(plotdata.trial{cfg.trl.current}(auxid(i),idx_current_view) ) ;
                    [B, A] = butter(5, 15/(data.fsample/2));
                    aux = filtfilt(B, A, (aux-mean(aux))');
                    aux = 0.5*aux/prctile(aux,99); % ~[0.5 0.5]
                    physio(i,:) = aux;
                end
            end
            
            
            cla(cfg.h_physio,'reset'); hold(cfg.h_physio,'on')
            set(cfg.h_physio,'XLim',cfg.win.tlim);
            set(cfg.h_physio,'YLim',[0.5 size(physio,1)+0.5] );
            for icnt=1:size(physio,1)
                plot(cfg.h_physio,t,physio(icnt,:)+icnt,'k')
            end
            
            
            
            set(cfg.h_physio,'YTickLabel',cfg.raw.ctype.(cfg.raw.ctypephysiolabel).label)
            set(cfg.h_physio,'XTick',''); 
        end
        cfg.redraw.physio = false;
    end
    %% main plot - add data types in the background
    
    if cfg.redraw.main
        
        subplot(cfg.h_main), cla,  
        auxYLimmain = [0.5, cfg.(cfg.disp.disptype).nsigblock+0.5];
        set(cfg.h_main,'YLim',auxYLimmain);
        set(cfg.h_main,'XLim',cfg.win.tlim);
        xlabel('Time [s]')
        
        
        % plot data types in the background
        aux = [cfg.dattype.background , setdiff(1:length(cfg.dattype.labels),cfg.dattype.background)];
        if ~ismember(cfg.dattype.current,cfg.dattype.background)
            aux = [setdiff(aux,cfg.dattype.current), cfg.dattype.current];
        end
        for iaty = aux
            % convert to samples of the current window.
            ibegwindow_rawdata = cfg.trl.trl(cfg.trl.current,1) + ...
                cfg.win.begs -1; %initial sample of the window in sample number of the raw continuous eeg data
            limart = cfg.dattype.(cfg.dattype.labels{iaty}) - ibegwindow_rawdata; %samples relative to window. = 0 if sample at the begin of the window
            tart = limart/data.fsample + t(1);
            iddiscard = tart(:,2) < t(1) | tart(:,1) > t(end);
            tart = tart(~iddiscard,:);
            
            for iart = 1:size(tart,1)
                fill([tart(iart,1) tart(iart,2)*[1 1] tart(iart,1)], ...
                    [auxYLimmain(1)*[1 1],auxYLimmain(2)*[1 1]],cfg.dattype.colors(iaty,:), ...
                    'EdgeColor','none','FaceAlpha',0.8)
            end
            
        end
        
        %% main plot - plot raw/ICA time series
        
        % Plot small 10Hz wave
        auxlims = round(length(t)*[0.051 0.19]);
        if length(auxlims)>1
            auxt = t(max(auxlims(1),1):auxlims(2));
            plot(auxt,0.5*cos(2*pi*10*auxt)+1.5,'k'); %,'LineWidth',1);
        end
        
        % Add units for raw data
        if strcmpi(cfg.disp.disptype,'raw')
            plot(0.05*t(end)+0.95*t(1)*[1,1],[1 2],'k','LineWidth',2)
            text(double(0.06*t(end)+0.94*t(1)),0,sprintf('%.0f %s',cfg.raw.amplitude_range,cfg.raw.unit))
        end
        
        %  prepare raw/ICA time series
        dataplot =  1 / cfg.(cfg.disp.disptype).amplitude_range * ...
            double(plotdata.trial{cfg.trl.current}(:,idx_current_view) ); %nallrawch x ntime
        
        switch cfg.disp.label
            case 'ica'
                dataplot = cfg.ica.w*dataplot(cfg.ica.userawch,:); %nallcomp x ntime
            case 'cleanraw'
                idcompbad = setdiff(1:cfg.ica.numofic,cfg.ica.ctype.(cfg.ica.ctypemainlabel).num);
                dataplot(cfg.ica.userawch,:) = dataplot(cfg.ica.userawch,:) - ...
                    cfg.ica.a(:,idcompbad)*cfg.ica.w(idcompbad,:)*dataplot(cfg.ica.userawch,:);
                dataplot(setdiff(1:length(cfg.raw.alllabels),cfg.ica.userawch),:) = nan;
            case 'badproj'
                idcompbad = setdiff(1:cfg.ica.numofic,cfg.ica.ctype.(cfg.ica.ctypemainlabel).num);
                dataplot(cfg.ica.userawch,:) = cfg.ica.a(:,idcompbad)*cfg.ica.w(idcompbad,:)*dataplot(cfg.ica.userawch,:);
                dataplot(setdiff(1:length(cfg.raw.alllabels),cfg.ica.userawch),:) = nan;
        end
        if cfg.raw.plotavgref && ~strcmpi(cfg.disp.label,'ica')
            idch = cfg.raw.ctype.(cfg.raw.ctypemainlabel).num;
            dataplot(idch,:) = bsxfun(@minus, dataplot(idch,:), mean(dataplot(idch,:),1));
        end
        
        

        auxcolor = cfg.(cfg.disp.disptype).mainplot.allchcolor(cfg.(cfg.disp.disptype).mainplot.showchid,:);
        idch_alldata = cfg.(cfg.disp.disptype).mainplot.allchid(cfg.(cfg.disp.disptype).mainplot.showchid);
        for i=1:length(idch_alldata)
            plot(cfg.h_main,t,dataplot(idch_alldata(i),:) + i,'Color',auxcolor(i,:))
        end
        
        
        
        % Highlight channel if any has been selected
        if ~isempty(cfg.(cfg.disp.disptype).mainplot.chhighlight)
            ich = find(cfg.(cfg.disp.disptype).mainplot.showchid == ...
                cfg.(cfg.disp.disptype).mainplot.chhighlight);
            if ~isempty(ich)
                fill([t(1) t(end)*[1 1] t(1)], ...
                    [(ich - 0.5)*[1 1],(ich + 0.5)*[1 1]],[1 0 0], ...
                    'EdgeColor','none','FaceAlpha',0.8)
                plot(t,dataplot(idch_alldata(ich),:) + ich,'Color','y')
            end
        end
        
        
        % plot events
        if isfield(cfg,'event')
            samplecorr = round(cfg.event.sample) - cfg.trl.trl(cfg.trl.current,1)+1;
            idok = find(samplecorr >= cfg.win.begs & samplecorr <= cfg.win.ends);            
            for iev = idok(:)'
               auxt = double(plotdata.time{cfg.trl.current}(samplecorr(iev)));
                plot(auxt*[1 1],auxYLimmain,'Color',cfg.event.color(iev,:));
                text(auxt,auxYLimmain(2),cfg.event.value{iev},...
                    'FontSize',9,'Color',0.2*[1 1 1],'VerticalAlignment','top','BackgroundColor',[1 1 1])
            end
        end
        

        cfg.redraw.main = false;
        set(cfg.h_main,'YTick',1:1:length(idch_alldata))
        set(cfg.h_main,'YTickLabel',cfg.(cfg.disp.disptype).mainplot.allchlabel(cfg.(cfg.disp.disptype).mainplot.showchid))

        
    end
    %% main plot - plot vertical bars of data selection

    if cfg.redraw.dattypelines
        delete(h_selectionlines)
        aux = [];
        if all(isfield(cfg.dattype,{'t_start_marker','t_stop_marker'})) && cfg.dattype.t_stop_marker > cfg.dattype.t_start_marker
            x = [cfg.dattype.t_start_marker*[1 1], cfg.dattype.t_stop_marker*[1 1]];
            y = [auxYLimmain ,auxYLimmain(2:-1:1)];
            aux.LineWidth = 2; aux.type = 'k--';
        else
            xstart = []; xstop = [];  
            if isfield(cfg.dattype,'t_start_marker')
                xstart =  cfg.dattype.t_start_marker*[1;1];
            end
            if isfield(cfg.dattype,'t_stop_marker')
                xstop = cfg.dattype.t_stop_marker*[1;1];
            end
            aux.LineWidth = 1; aux.type = 'k';
            x = [xstart, xstop]; y = repmat(auxYLimmain(:),1,size(x,2));
        end
        h_selectionlines = plot(cfg.h_main,x,y,aux.type,'LineWidth',aux.LineWidth);
        
        cfg.redraw.dattypelines = false;
    else
        uistack(h_selectionlines,'top');
    end
    
    
    
    %% plot bar in the uicontrol time section
    
    if cfg.redraw.bar
        nsamplestrl = length(plotdata.time{cfg.trl.current});
        sampletype  = ones(1,nsamplestrl,3);
        sampletype(1,idx_current_view,1)= 1; sampletype(1,idx_current_view,2)= 1; sampletype(1,idx_current_view,3)= 0;
        for iaty = [cfg.dattype.background cfg.ica.discarddattypeid(:)']
            if ismember(iaty, cfg.dattype.background)
                auxcolor = cfg.dattype.colors(iaty,:);
            else
               auxcolor = [1 0 0];
            end
            
            
            % convert to samples of the current window.
            %initial sample of the window in sample number of the raw continuous eeg data
            limart_trl = cfg.dattype.(cfg.dattype.labels{iaty}) - cfg.trl.trl(cfg.trl.current,1) +1; %samples relative to TRIAL. = 1 if sample at the begin of the TRIAL
            useart = limart_trl(:,1)>0 & limart_trl(:,1)<= nsamplestrl | ...
                limart_trl(:,2)>0 & limart_trl(:,2)<= nsamplestrl | ...
                limart_trl(:,1)<1 & limart_trl(:,2)>= nsamplestrl;
            limart_trl = limart_trl(useart,:);
            limart_trl(limart_trl<1)=1;
            limart_trl(limart_trl>nsamplestrl)=nsamplestrl;
            for iart = 1:size(limart_trl,1)
                sampletype(1,limart_trl(iart,1):limart_trl(iart,2),1)= auxcolor(1);
                sampletype(1,limart_trl(iart,1):limart_trl(iart,2),2)= auxcolor(2);
                sampletype(1,limart_trl(iart,1):limart_trl(iart,2),3)= auxcolor(3);
            end
        end
        
        subplot(cfg.h_bar); cla(cfg.h_bar); 
        imagesc(sampletype,[0 1])
        %imagesc(cfg.h_bar,sampletype); set(cfg.h_bar,'CLim',[0 1],'XLim',[1 length(sampletype)]); 
        set(cfg.h_bar,'XLim',[1 length(sampletype)],'YLim',[0.5 1.5]); 
        box(cfg.h_bar,'on'); 
        cfg.redraw.bar = false;
    end
    %% plot ica topoplots if wanted
    
    if cfg.figlay.boxtopo.Value ==1 && strcmpi(cfg.disp.disptype,'ica')
        set(cfg.h_topo ,'Layer', 'Top')
        if cfg.redraw.topo
            disp('redraw topo')
            subplot(cfg.h_topo); cla; box on
            set(cfg.h_topo,'CLim',[-1 1])
            ylim(auxYLimmain)
            idok = cfg.raw.haspos(cfg.ica.userawch);
            cpos = cfg.raw.allpos(cfg.ica.userawch(idok),:);
            if size(cpos,1)>1
                % scale pos to -0.5 0.5
                cpos(:,2) = cpos(:,2) - mean(cpos(:,2));
                cpos = 0.4*cpos/max(abs(cpos(:,2)));
                cpos2 = cpos;
                for ic = 1:length(idch_alldata)
                    cpos2(:,2) = cpos(:,2) + ic;
                    vplot =  cfg.ica.a(idok,idch_alldata(ic))';
                    mytopoplot_small(vplot/max(abs(vplot)),cpos2,[],[],40);
                end
                colorbar off;
                box on
            end
            cfg.redraw.topo = 0;
        end            
    else
        axis(cfg.h_topo,'off')
        cla(cfg.h_topo); 
    end
    
    %% plot title
    
    if size(cfg.trl.trl,2)>3
        mystring = sprintf('%s -Cond %i',cfg.disp.disptype,cfg.trl.trl(cfg.trl.current,4));
    else
        mystring = sprintf('%s',cfg.disp.disptype);
    end
    set(cfg.h_title,'String',mystring)
    set(cfg.h_title,'Interpreter','none')
    
    %% save data and continue
    subplot(cfg.h_main)
    guidata(fig,cfg);
    
    
    uiwait;
    if(ishandle(fig)), cfg=guidata(fig); end
    

    
    if cfg.stop;
        break
    end % exit scrolling
end


%% finished - save data if wanted and leave
if ishandle(fig)
    close(fig)
    

    mysave(cfg,'final',data)
    

else
    fprintf('Data not saved! Click <done> to do so!\n')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the SUBFUNCTIONS start taht implement the gui callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next window
function cfg = nextwindow(h, eventdata)
cfg=guidata(h);
trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
if cfg.win.ends >= (trllength_w-1) && cfg.trl.current < cfg.trl.n_trl
    cfg.trl.last = cfg.trl.current;
    cfg.trl.current = cfg.trl.current+1;
    cfg.trl.currentbox.String = num2str(cfg.trl.current);
    trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
    cfg.win.begs=1;
    cfg.win.ends = min([trllength_w cfg.win.n]);
    cfg = cleandatselmarkers(cfg);
elseif cfg.win.ends < trllength_w
    cfg.win.begs = min(cfg.win.begs+cfg.win.n  ,   trllength_w);
    cfg.win.ends = min([cfg.win.begs + cfg.win.n - 1,  trllength_w]);
end
cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
if ~isempty(eventdata),
    guidata(h,cfg);
    uiresume;
end


end

% Previous window
function cfg= prevwindow(h, eventdata)
cfg=guidata(h);
if cfg.win.begs==1
    cfg.trl.last = cfg.trl.current;
    cfg.trl.current = max(cfg.trl.current-1,1);
    cfg.trl.currentbox.String = num2str(cfg.trl.current);
    if cfg.trl.last~=cfg.trl.current
        cfg = cleandatselmarkers(cfg);
        trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
        cfg.win.ends = trllength_w;
        cfg.win.begs = max([1 cfg.win.ends - cfg.win.n+1]);
    end
else
    cfg.win.begs = max(cfg.win.begs - cfg.win.n ,1);
    trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
    cfg.win.ends = min([cfg.win.begs + cfg.win.n - 1,  trllength_w]);
end
cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
if ~isempty(eventdata),
    guidata(h,cfg);
    uiresume;
end
end

% Previous trial
function prevtrial(h, ~)
cfg=guidata(h);
if cfg.trl.current > 1
    cfg.trl.last = cfg.trl.current;
    cfg.trl.current = cfg.trl.current-1;
    cfg.trl.currentbox.String = num2str(cfg.trl.current);
    trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
    cfg.win.begs=1;
    cfg.win.ends = min([trllength_w cfg.win.n]);
    cfg = cleandatselmarkers(cfg);
    cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
end
guidata(h,cfg);
uiresume;
end

% Next trial
function nexttrial(h, ~)
cfg=guidata(h);
if cfg.trl.current < size(cfg.trl.trl,1)
    cfg.trl.last = cfg.trl.current;
    cfg.trl.current = cfg.trl.current+1;
    cfg.trl.currentbox.String = num2str(cfg.trl.current);
    trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
    cfg.win.begs=1;
    cfg.win.ends = min([trllength_w cfg.win.n]);
    cfg = cleandatselmarkers(cfg);
    cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
end
guidata(h,cfg);
uiresume;
end

function trialnumber(h,eventdata)
cfg=guidata(h);
targetnumber = str2double(eventdata.Source.String);
if eventdata.Source.String ~= cfg.trl.currentbox.String %borrar luego
    error('no coinciden')
end
if targetnumber<=0 || targetnumber>size(cfg.trl.trl,1)
    fprintf('Cannot access trial %i - only %i trials\n',targetnumber,size(cfg.trl.trl,1));
    cfg.trl.currentbox.String = num2str(cfg.trl.current);
else
    cfg.trl.current = targetnumber;
    cfg.trl.currentbox.String = num2str(cfg.trl.current);
    trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
    cfg.win.begs=1;
    cfg.win.ends = min([trllength_w cfg.win.n]);
    cfg = cleandatselmarkers(cfg);
end
cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
guidata(h,cfg);
uiresume;
end

function changewindowt(h,eventdata)
cfg=guidata(h);
targett = str2double(eventdata.Source.String);
if targett<=0 || isnan(targett)
    targett = cfg.win.tdefault;
    cfg.win.tbox.String = num2str(targett);
end
cfg.win.t = targett;
cfg.win.n = round(cfg.win.t*cfg.fsample);
trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
cfg.win.ends = min([trllength_w, (cfg.win.begs+cfg.win.n-1)]);
cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
guidata(h,cfg);
uiresume;
end


function changeverts(h,eventdata)
cfg=guidata(h);
action = eventdata.Source.String;
if strcmpi(action,'>')
    targetscale = cfg.(cfg.disp.disptype).amplitude_range*1.2;
elseif strcmpi(action,'<')
    targetscale = cfg.(cfg.disp.disptype).amplitude_range*0.8;
else
    targetscale = str2double(eventdata.Source.String);
end
if targetscale<=0 || isnan(targetscale)
    targetscale = cfg.(cfg.disp.disptype).amplitude_range_default;
end
if mod(targetscale,1)==0
    cfg.win.verts.String = sprintf('%.1e',targetscale);
else
    cfg.win.verts.String = sprintf('%.1e',targetscale);
end
cfg.(cfg.disp.disptype).amplitude_range = targetscale;
cfg.redraw.main = true; cfg.redraw.dattypelines=true;
guidata(h,cfg);
uiresume;
end


function changechanind(h,eventdata)
cfg=guidata(h);
nallch = length(cfg.(cfg.disp.disptype).mainplot.allchid);
cfg.(cfg.disp.disptype).mainplot.showchid = sort(cfg.(cfg.disp.disptype).mainplot.showchid,'ascend');
switch eventdata.Source.Tag
    case '>'
        begch = min([cfg.(cfg.disp.disptype).mainplot.showchid(1)+cfg.(cfg.disp.disptype).nsigblock, ...
            (nallch - cfg.(cfg.disp.disptype).nsigblock +1)]);
        begch = max([begch 1]);
        endch = min([begch + cfg.(cfg.disp.disptype).nsigblock-1, nallch]);
    case '<'
        begch = max([cfg.(cfg.disp.disptype).mainplot.showchid(1)-cfg.(cfg.disp.disptype).nsigblock , 1]);
        endch = min([begch + cfg.(cfg.disp.disptype).nsigblock-1, nallch]);
    case 'begch'
        begch = str2double(eventdata.Source.String);
        if isnan(begch)
            begch = cfg.(cfg.disp.disptype).mainplot.showchid(1);
        else
            begch = max(round(begch),1);
            endch = cfg.(cfg.disp.disptype).mainplot.showchid(end);
        end
    case 'endch'
        endch = str2double(eventdata.Source.String);
        if isnan(endch)
            endch = cfg.(cfg.disp.disptype).mainplot.showchid(end);
        else
            endch = min(round(endch),nallch);
            begch = cfg.(cfg.disp.disptype).mainplot.showchid(1);
        end
end

cfg.(cfg.disp.disptype).mainplot.showchid = endch:-1:begch;
cfg.(cfg.disp.disptype).nsigblock = length(cfg.(cfg.disp.disptype).mainplot.showchid);
cfg.win.mainplotbegch.String = num2str(begch);
cfg.win.mainplotendch.String = num2str(endch);
cfg.redraw.main = true; cfg.redraw.topo = true; cfg.redraw.dattypelines = true;
guidata(h,cfg);

uiresume;
end



function channelselect(h,~)
cfg=guidata(h);
cfg.(cfg.disp.disptype) = mygui_channelselect(cfg.(cfg.disp.disptype));
if strcmpi(cfg.disp.disptype,'raw')
    cfg.win.mainplotbegch.String = num2str(min(cfg.raw.mainplot.showchid));
    cfg.win.mainplotendch.String = num2str(max(cfg.raw.mainplot.showchid));
end
cfg.redraw.main = true; cfg.redraw.topo = true; cfg.redraw.dattypelines = true;
guidata(h,cfg);
uiresume;
end





function dattype_buttonpush(h, eventdata)
cfg=guidata(h);
cfg.dattype.current = find(strcmpi(eventdata.Source.String, cfg.dattype.labels));
% get this data type activated and deactivate the others
cfg.dattype.buttontype{cfg.dattype.current}.Value=1;
cfg.dattype.buttontype{cfg.dattype.current}.FontWeight= 'bold';
for iaty = setdiff(1:length(cfg.dattype.labels),cfg.dattype.current)
    cfg.dattype.buttontype{iaty}.Value = 0;
    cfg.dattype.buttontype{iaty}.FontWeight= 'normal';
end
guidata(h,cfg);
uiresume;
end



% Quit
function stop(h, ~)
cfg=guidata(h);
cfg.stop=1;
guidata(h,cfg);
uiresume;
end

% Allow key - controll
function keypress(h, ~)
disp('here')
cfg=guidata(h);
key= double(get(h, 'CurrentCharacter'));
if  key
    switch key
        case 28     % arrow left
            cfg=prevwindow(h, []);
        case 29     % arrow right
            cfg=nextwindow(h, []);
        case 30
            aux.Source.String = cfg.disp.alllabel;
            aux.Source.Value = mod(cfg.disp.icur,length(cfg.disp.alllabel)) +1;
            aux.EventName = 'keypress';
            cfg=selectdisplaymag(h, aux);
        case 31
            aux.Source.String = cfg.disp.alllabel;
            aux.Source.Value = mod(cfg.disp.icur+2,length(cfg.disp.alllabel)) +1;
            aux.EventName = 'keypress';
            cfg=selectdisplaymag(h, aux);
        case 99 % 'c'
            num = round(h.CurrentAxes.CurrentPoint(1,2));
            if num > 0 && num <= length(cfg.(cfg.disp.disptype).mainplot.showchid)
                num = cfg.(cfg.disp.disptype).mainplot.showchid(num);
                if num == cfg.(cfg.disp.disptype).mainplot.chhighlight
                    cfg.(cfg.disp.disptype).mainplot.chhighlight = [];
                else
                    cfg.(cfg.disp.disptype).mainplot.chhighlight = num;
                end
                cfg.redraw.main = true; cfg.redraw.dattypelines = true;
            else
                fprintf('not selecting any channel\n')
            end
    end
end
set(h, 'CurrentCharacter', '-');
guidata(h, cfg);
uiresume(h);
end

function mouse_click(h,eventdata)

cfg=guidata(h);
value = h.CurrentAxes.CurrentPoint(1,1);
axbarselected = all(cfg.h_bar.Position == h.CurrentAxes.Position);
switch get(h,'Selectiontype');
    case 'normal' % left click;
        if axbarselected
            trllength_w = diff(cfg.trl.trl(cfg.trl.current,1:2)) +1;
            cfg.win.begs = min( max(round(value),1), trllength_w-2);
            cfg.win.ends = min([cfg.win.begs + cfg.win.n - 1,  trllength_w]);
            cfg.redraw.physio = true; cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true;
        else
            cfg.dattype.t_start_marker = max([value, cfg.win.tlim(1)]) ;
            set(cfg.dattype.boxstart,'String',sprintf('%.2f',cfg.dattype.t_start_marker));
            cfg.redraw.dattypelines = true;
        end
    case 'alt'% right-click;
        cfg.dattype.t_stop_marker = min([value, cfg.win.tend]);
        set(cfg.dattype.boxstop,'String',sprintf('%.2f',cfg.dattype.t_stop_marker));
        cfg.redraw.dattypelines = true;
    case 'extend' %middle button - add new data type segment
        button_dattype_change(h,eventdata)
        cfg=guidata(h);
end
guidata(h,cfg);
uiresume;
end

function start_datasel(h,~)
cfg=guidata(h);
cfg.dattype.t_start_marker = str2double(get(h,'String'));
cfg.redraw.dattypelines = true; 
guidata(h,cfg);
uiresume;
end

function stop_datasel(h,~)
cfg=guidata(h);
cfg.dattype.t_stop_marker = str2double(get(h,'String'));
cfg.redraw.dattypelines = true; 
guidata(h,cfg);
uiresume;
end



function button_dattype_change(h,eventdata)
cfg=guidata(h);
if ~isfield(cfg.dattype,'t_stop_marker')
    fprintf('End of selection not defined \n')
elseif ~isfield(cfg.dattype,'t_start_marker')
    fprintf('Begin of selection not defined \n')
elseif cfg.dattype.t_start_marker >= cfg.dattype.t_stop_marker
    fprintf('End of selection should be after the begin \n')
else
    % transform times to samples in the raw continuous eeg data
    tart = [cfg.dattype.t_start_marker cfg.dattype.t_stop_marker];
    limart = round((tart-cfg.win.tlim(1))*cfg.fsample);
    ibegwindow_rawdata = cfg.trl.trl(cfg.trl.current,1) + ...
        cfg.win.begs -1; %initial sample of the window in sample number of the raw continuous eeg data
    newart = limart + ibegwindow_rawdata;
    art = cfg.dattype.(cfg.dattype.labels{cfg.dattype.current});
    
    if strcmpi(eventdata.EventName,'WindowMousePress')
        action = 'add';
    else
        action = eventdata.Source.String;
    end
    
    switch action
        
        case 'del'
            keepart = true(size(art,1),1);
            extraart = [];
            for iart = 1:size(art,1)
                if newart(1)<=art(iart,1) && newart(2)>= art(iart,2)
                    keepart(iart) = false; %old artifact included in the new box
                elseif newart(1)> art(iart,1) && newart(2)<art(iart,2) %new box inside old artifact
                    extraart = [newart(2)+1 art(iart,2)];
                    art(iart,2) = newart(1)-1;
                elseif newart(1)<art(iart,1) && newart(2)>=art(iart,1)
                    art(iart,1) = newart(2)+1;
                elseif newart(1)<=art(iart,2) && newart(2)>art(iart,2)
                    art(iart,2) = newart(1)-1;
                end
            end
            art = art(keepart,:);
            art = [art;extraart];
            [~,goodorder] = sort(art(:,1),'ascend');
            art = art(goodorder,:);
            fprintf('deleting %s - %i artifacts \n',cfg.dattype.labels{cfg.dattype.current},size(art,1));
            
        case 'add'
            art = [art; newart];
            % corrects for overlapping between artifacts
            [~,goodorder] = sort(art(:,1),'ascend');
            art = art(goodorder,:);
            % Checks if one artifact overlaps with the following
            keepart = true(size(art,1),1);
            for iart = 1:size(art,1)-1
                if art(iart,2)+1 >= art(iart+1,1)
                    art(iart+1,1) = art(iart,1);
                    art(iart+1,2) = max([art(iart,2) art(iart+1,2)]);
                    keepart(iart)=false;
                end
            end
            if all(keepart)
                fprintf('adding %s ',cfg.dattype.labels{cfg.dattype.current})
            else
                art = art(keepart,:);
                fprintf('fusing %s ',cfg.dattype.labels{cfg.dattype.current})
            end
            fprintf('- %i segments \n',size(art,1));
    end
    
    cfg.dattype.(cfg.dattype.labels{cfg.dattype.current}) = art;
    cfg = cleandatselmarkers(cfg);
    cfg.redraw.main = true; cfg.redraw.dattypelines = true; cfg.redraw.bar = true; 
    guidata(h,cfg);
end

uiresume;
end




function cfg=selectdisplaymag(h, eventdata)
cfg=guidata(h);
targetdisp = eventdata.Source.String{eventdata.Source.Value};
if ~isfield(cfg.ica,'w')
    fprintf('Cannot display %s - ICA not computed\n',targetdisp);
    cfg.disp.guiselect.Value = cfg.disp.icur;
else
    cfg.disp.icur = eventdata.Source.Value;
    cfg.disp.label = targetdisp;
    cfg.disp.guiselect.Value = eventdata.Source.Value;
    if strcmpi(cfg.disp.label,'ica')
        cfg.disp.disptype = 'ica';
        cfg.disp.guichanselect.String = 'components';
    else
        cfg.disp.disptype = 'raw';
        cfg.disp.guichanselect.String = 'channels';
    end
end

cfg.win.mainplotbegch.String = num2str(min(cfg.(cfg.disp.disptype).mainplot.showchid));
cfg.win.mainplotendch.String = num2str(max(cfg.(cfg.disp.disptype).mainplot.showchid));
cfg.redraw.main = true; cfg.redraw.topo = true; cfg.redraw.dattypelines = true;
if ~strcmpi(eventdata.EventName,'keypress'); guidata(h,cfg); uiresume; end
end



function cfg = cleandatselmarkers(cfg)
if isfield(cfg.dattype,'t_start_marker')
    cfg.dattype = rmfield(cfg.dattype,'t_start_marker');
    cfg.dattype.boxstart.String = '-';
end
if isfield(cfg.dattype,'t_stop_marker')
    cfg.dattype = rmfield(cfg.dattype,'t_stop_marker');
    cfg.dattype.boxstop.String = '-';
end
cfg.redraw.dattypelines = true;
end


function topo_segment(h, ~,plotdata)
cfg=guidata(h);
if ~isfield(cfg.dattype,'t_stop_marker')
    fprintf('End of selection not defined \n')
elseif ~isfield(cfg.dattype,'t_start_marker')
    fprintf('Begin of selection not defined \n')
elseif cfg.dattype.t_start_marker >= cfg.dattype.t_stop_marker
    fprintf('End of selection should be after the begin \n')
else
    tart = [cfg.dattype.t_start_marker cfg.dattype.t_stop_marker];
    limwin = round((tart-cfg.win.tlim(1))*cfg.fsample) + cfg.win.begs;
    idch = intersect(cfg.raw.ctype.(cfg.raw.ctypemainlabel).num,find(cfg.raw.haspos));
    
    
    magdisp = cfg.disp.label; 
    if strcmpi(magdisp,'ica'), magdisp = 'cleanraw'; end
    dataplot = plotdata.trial{cfg.trl.current}(:,limwin(1):limwin(2));
    switch magdisp
        case 'cleanraw'
            idcompbad = setdiff(1:cfg.ica.numofic,cfg.ica.ctype.(cfg.ica.ctypemainlabel).num);
            dataplot(cfg.ica.userawch,:) = dataplot(cfg.ica.userawch,:) - ...
                cfg.ica.a(:,idcompbad)*cfg.ica.w(idcompbad,:)*dataplot(cfg.ica.userawch,:);
            dataplot(setdiff(1:length(cfg.raw.alllabels),cfg.ica.userawch),:) = nan;
        case 'badproj'
            idcompbad = setdiff(1:cfg.ica.numofic,cfg.ica.ctype.(cfg.ica.ctypemainlabel).num);
            dataplot(cfg.ica.userawch,:) = cfg.ica.a(:,idcompbad)*cfg.ica.w(idcompbad,:)*dataplot(cfg.ica.userawch,:);
            dataplot(setdiff(1:length(cfg.raw.alllabels),cfg.ica.userawch),:) = nan;
    end
    dataplot = dataplot(idch,:);
    if cfg.raw.plotavgref
        dataplot = bsxfun(@minus, dataplot, mean(dataplot,1));
    end
    
    
    topoint = mean(abs(dataplot).^2,2);
    cpos = cfg.raw.allpos(idch,:);
    if size(cpos,1)>1
        figure; hold on
        mytopoplot(topoint,cpos,[],[],150)
        set(gca,'CLim',[min(topoint) max(topoint)]);
        colormap('jet');
        title(sprintf('%s - %.1f-%.1fsec',magdisp,tart(1),tart(2)));
    else
        fprintf('No channel positions\n')
    end
    
end
uiresume;
end


% output file
function cfg=changeoutfile(h, ~)
cfg=guidata(h);
answer = inputdlg('Output file name:','Output file name:',1,{cfg.fileout});
if ~isempty(answer)
    cfg.fileout = answer{1};
    cfg.outfiletext.String = cfg.fileout;
    guidata(h,cfg);
end
uiresume;
end


% Add comments
function cfg=addcomments(h, ~)
cfg=guidata(h);
if isempty(cfg.comments)
    answer = inputdlg('Comments:','Comments:',[1 200]);
else
    answer = inputdlg('Comments:','Comments:',[1 200],{cfg.comments});
end
if ~isempty(answer)
    cfg.comments = answer{1};
    guidata(h,cfg);
end
uiresume;
end




function ICA_minitopo(h, ~)
cfg=guidata(h);
cfg.redraw.topo=true;
guidata(h,cfg);
uiresume;
end


function mysave(h, eventdata,data)
if ~ischar(eventdata), cfg=guidata(h); else cfg = h; end
fprintf('Saving information ..')
data.cfg.ica = cfg.ica;
data.cfg.ica.discarddattype = cfg.dattype.labels(cfg.ica.discarddattypeid);
data.cfg.raw = cfg.raw;
if isfield(data.cfg,'chaninfo')
    for itype = 1:length(cfg.raw.ctypelabels)
        data.cfg.chaninfo.chantype(cfg.raw.ctype.(cfg.raw.ctypelabels{itype}).num) = cfg.raw.ctypelabels(itype);
    end
elseif isfield(data,'chaninfo') && isfield(data.chaninfo,'type')
    for itype = 1:length(cfg.raw.ctypelabels)
        data.chaninfo.type(cfg.raw.ctype.(cfg.raw.ctypelabels{itype}).num) = cfg.raw.ctypelabels(itype);
    end
elseif isfield(data,'chaninfo') && isfield(data.chaninfo,'chantype')
    for itype = 1:length(cfg.raw.ctypelabels)
        data.chaninfo.chantype(cfg.raw.ctype.(cfg.raw.ctypelabels{itype}).num) = cfg.raw.ctypelabels(itype);
    end
end




%data.cfg.trl = cfg.trl;
data.cfg.dattype = [];
% data.cfg.dattype.labels = cfg.dattype.labels;
for iaty = 1:length(cfg.dattype.labels)
    data.cfg.dattype.(cfg.dattype.labels{iaty}) = cfg.dattype.(cfg.dattype.labels{iaty});
end
data.cfg.fileout = cfg.fileout;
data.cfg.comments = cfg.comments;
data.cfg.datescroll = datestr(now, 'yyyymmdd');
data.cfg.scrolltool_ver = '1.0';
save(cfg.fileout,'data','-v7.3')
fprintf('Done.\n');
if ~ischar(eventdata), guidata(h,cfg); uiresume; end

end

