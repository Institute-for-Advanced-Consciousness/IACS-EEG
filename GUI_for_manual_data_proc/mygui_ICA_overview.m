
function mygui_ICA_overview(h, ~,data)
cfg=guidata(h);
if ~isfield(cfg.ica,'a') 
    fprintf('ICA has not been computed\n')
elseif sum(cfg.raw.haspos(cfg.ica.userawch))<3
    fprintf('no channel positions\n')
else
    %% CREATE FIGURE
    
    pos      = get(0,'DefaultFigurePosition');
    pos(3:4) = [800 400];
    dlg = figure('Position', pos); %,'WindowStyle','modal');
    %     dlg      = dialog('Name', 'ICA overview', 'Position', pos);
    
    optplot = [];
    optplot.ica = cfg.ica;
    if isfield(cfg.ica,'artifactpower') && isstruct(cfg.ica.artifactpower)
        optplot.artnames = fieldnames(cfg.ica.artifactpower);
        optplot.artnames = optplot.artnames(:)';
    else
        optplot.artnames = [];
    end
    optplot.topoint = cfg.ica.a(:,:); %channels x components
    optplot.cpos = cfg.raw.allpos(cfg.ica.userawch,:);
    idok = cfg.raw.haspos(cfg.ica.userawch);
    optplot.cpos = optplot.cpos(idok,:);
    optplot.topoint = optplot.topoint(idok,:);
    
    optplot.dispselect = uicontrol(dlg,'units','normalized','position',[0.1  0.9  0.3  0.05], ...
        'Style', 'popupmenu', 'String', [{'original'}, optplot.artnames] ,'Value',1,'Callback','uiresume');
    
    
    set(dlg,'WindowButtonDownFcn',@mouse_click);
    param.replot = false;
    iplot = optplot.dispselect.Value;
    
    param.ica = cfg.ica;
    param.raw = cfg.raw;
    param.ordering = 1:cfg.ica.numofic;
    param.currentc = [];
    param.plotc = [];
    
    param.currentlabels = nan(param.ica.numofic,1);
    for ity = 1:length(param.ica.ctypelabels)
        param.currentlabels(param.ica.ctype.(param.ica.ctypelabels{ity}).num) = ity;
    end
    titlecolors = lines(length(param.ica.ctypelabels));
    
    
    %% plot topoplots in initial ordering
    
    ncol = floor(sqrt(cfg.ica.numofic))+4;
    nrow = ceil(cfg.ica.numofic/ncol)+2;
    width = 0.9/ncol;
    height = 0.7/nrow;
    subplotpos = nan(cfg.ica.numofic,4);
    subplotaxes = cell(cfg.ica.numofic,1);
    set(dlg,'Visible','off')
    for ic = 1:cfg.ica.numofic
        icol = mod(ic-1,ncol)+1;
        irow = ceil(ic/ncol);
        subplotpos(ic,:) = [(icol-1)/ncol (nrow-irow-2)/nrow width height];
        subplotaxes{ic} = axes('Position',subplotpos(ic,:)); axis off;
        mytopoplot_small(optplot.topoint(:,ic)/max(abs(optplot.topoint(:,ic))), ...
            optplot.cpos,[],[],40)
        set(subplotaxes{ic},'CLim',[-1 1]);
        colorbar off;
        colormap('jet');
        % label components with proportion of total power 
        title(sprintf('%.3f',cfg.ica.vmax(ic)./sum(cfg.ica.vmax)),'Color',...
            titlecolors(param.currentlabels(ic),:),'FontSize',10);
        axis equal
    end
    set(dlg,'Visible','on')
    param.subplotpos = subplotpos;
    %% Add layout overview
    auxy = (nrow-2)/nrow;
    aux = [(ncol-2)/ncol (nrow-2)/nrow 2*width 1-auxy];
    axes('Position',aux); axis off;
    mytopoplot_small(nan(size(optplot.cpos,1),1), ...
        optplot.cpos,[],[],40); colorbar off;
    plot(optplot.cpos(:,1),optplot.cpos(:,2),'k.','MarkerSize',5);
    aux = cfg.raw.alllabels(cfg.ica.userawch); aux =  aux(idok);
    for ic = 1:size(optplot.cpos,1)
        text(optplot.cpos(ic,1),optplot.cpos(ic,2),aux{ic},'FontSize',7)
    end
    
    
    
    %% PREPARE POWERSPECTRUM A COMPONENT IS SELECTED
    % Create time series
    alldat = nan(length(cfg.ica.userawch),max(cfg.trl.trl(:,2)));
    for itrial = 1:length(data.trial)
        alldat(:,cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2)) = ...
            data.trial{itrial}(cfg.ica.userawch,:);
    end
    % Replace artifacts with nans
    usesample = true(1,max(cfg.trl.trl(:,2)));
    for iaty = cfg.ica.discarddattypeid
        aux =  cfg.dattype.(cfg.dattype.labels{iaty});
        for iart = 1:size(aux,1)
            usesample(aux(iart,1):aux(iart,2)) = false;
        end
    end
    alldat(:,~usesample) = nan;
    [foi,rawpowspctrm] = my_compute_powerspctrm([],cfg.fsample,alldat);
    
    
    
    %%
    
    dlg.UserData = param;
    while(ishandle(dlg))
        param = dlg.UserData;
        % Change ordering
        if param.replot || iplot ~= optplot.dispselect.Value
            iplot = optplot.dispselect.Value;
            if iplot == 1
                ordering = 1:cfg.ica.numofic;
            else
                aux = cfg.ica.artifactpower.(optplot.artnames{iplot-1});
                [ratval,ordering] = sort(aux,'descend');
            end
            for ic = 1:length(ordering)
                set(subplotaxes{ordering(ic)},'Position',subplotpos(ic,:));
                if iplot == 1
                    mystr = sprintf('%s',cfg.ica.alllabels{ordering(ic)});
                else
                    mystr = sprintf('%s (%.1f)',cfg.ica.alllabels{ordering(ic)},ratval(ic));
                end
                title(subplotaxes{ordering(ic)}, mystr,'Color',...
                    titlecolors(param.currentlabels(ordering(ic)),:));
            end
            param.ordering = ordering;
            param.replot = false;
        end
        
        if ~isempty(param.plotc)
            %% INSPECT COMPONENT
            
            
            ic = param.plotc;
            powraw =  rawpowspctrm * abs(cfg.ica.a(:,ic));
            [~,powspctrmwithout] = my_compute_powerspctrm([],cfg.fsample,...
                alldat - cfg.ica.a(:,ic)*cfg.ica.w(ic,:)*alldat);
            powrawcorr = powspctrmwithout * abs(cfg.ica.a(:,ic));
            
            [~,powic] = my_compute_powerspctrm([],cfg.fsample,...
                cfg.ica.w(ic,:)*alldat);
            
            powdiff = 1-powrawcorr./powraw;
            powraw = powraw/max(powraw);
            powic = powic/max(powic);
            
            %% CREATE FIGURE FOR EXTRA PLOT
            
            posnewfig      = get(0,'DefaultFigurePosition');
            posnewfig(3:4) = [800 400];
            % dlg      = dialog('Name', 'Channel selection', 'Position', pos);
            dlg2 = figure('Position', posnewfig);
            %set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears
            
            axestopo  = axes('Position',[0.0 0.2 0.45 0.7],'Parent',dlg2);
            axespow   = axes('Position',[0.5 0.3 0.45 0.6],'Parent',dlg2); axis on
            
            
            uicontrol(dlg2,'units','normalized','position',[0.9  0.03  0.08  0.05], ...
                'Style', 'togglebutton', 'String', 'xlog','Value',0,'callback', {@xlog,axespow});
            
            uicontrol(dlg2,'units','normalized','position',[0.6  0.08  0.08  0.05], ...
                'Style', 'text', 'String', 'fmin');
            optfreq.box.fmin = uicontrol(dlg2,'units','normalized','position',[0.6  0.03  0.08  0.05], ...
                'Style', 'edit', 'tag','fmin','String','2','callback', {@foichange,axespow});
            
            uicontrol(dlg2,'units','normalized','position',[0.7  0.08  0.08  0.05], ...
                'Style', 'text', 'String', 'fmax');
            optfreq.box.fmax = uicontrol(dlg2,'units','normalized','position',[0.7  0.03  0.08  0.05], ...
                'Style', 'edit','tag','fmax', 'String',  '45','callback', {@foichange,axespow});
            
            %% PLOT TOPOPLOT
            topoint = cfg.ica.a(:,ic);
            cpos = cfg.raw.allpos(cfg.ica.userawch,:);
            idok = cfg.raw.haspos(cfg.ica.userawch);
            topoint = topoint(idok); cpos = cpos(idok,:);
            interpmethod = 'v4';
            clim = max(abs(topoint))*[-1 1];
            
            subplot(axestopo); cla; hold on
            mytopoplot(topoint,cpos,interpmethod,[],100)
            set(axestopo,'CLim',clim);
            colorbar off;
            colormap('jet');
            title(sprintf('comp%i',ic))
            %% PLOT POWERSPECTRUM
            subplot(axespow); cla; hold on
            plot(foi,powic,'b','LineWidth',2);
            
            [aux,p1,p2] = plotyy(foi,powraw,foi,powdiff);
            set(aux,{'ycolor'},{'k';[0 0.8 0]});
            set(p1,'color','r'); set(p2,'color',[0 0.8 0]);
            %plot(foi,powraw,'r');
            %yyaxis(axespow,'right');
            %plot(foi,powdiff);
            %plot([min(foi) max(foi)],[0 0],'Color',0.7*[1 1 1])
            %yyaxis(axespow,'left');
            
            
            
            legend('IC power','raw power weighted by a','power decrease (%) weighted by a')
            
            xlabel('frequency(Hz)')
            ylabel('Power spectrum (arbitrary units)');
            title(sprintf('comp%i',ic))
            xlim([2 45])
            %%
            param.plotc = [];
        end
        
        dlg.UserData = param;
        
        
        
        
        uiwait(dlg);
    end
    
    cfg.ica.ctype = param.ica.ctype;
    idbad = setdiff(1:cfg.ica.numofic,cfg.ica.ctype.(cfg.ica.ctypemainlabel).num);
    cfg.ica.mainplot.allchcolor(idbad,:) = repmat(0.6*[1 1 1],length(idbad),1);
    cfg.redraw.main = true; cfg.redraw.dattypelines = true;
    guidata(h,cfg);
    

end

uiresume;
end




function mouse_click(dlg,~)

param = dlg.UserData;

curpos = dlg.CurrentAxes.Position;
[vmin,iplot] = min(sum(abs(param.subplotpos - repmat(curpos,param.ica.numofic,1)),2));
if vmin<0.001
    icomp = param.ordering(iplot);
else
    icomp = [];
end
if isempty(param.currentc) || isempty(icomp)
    param.currentc = icomp;
elseif param.currentc == icomp
    switch get(dlg,'Selectiontype')
        case {'normal','open'} % left click or double clicks
            param.plotc = icomp;
        case 'alt'% right-click;
            [Selection,ok] = listdlg('PromptString',sprintf('Component %i type',icomp),...
                'SelectionMode','single','InitialValue',param.currentlabels(icomp),...
                'ListString',param.ica.ctypelabels);
            if ok && ~isempty(Selection)
                param.currentlabels(icomp) = Selection;
                for ity = 1:length(param.ica.ctypelabels)
                    num = find(param.currentlabels==ity);
                    param.ica.ctype.(param.ica.ctypelabels{ity}).num = num;
                    param.ica.ctype.(param.ica.ctypelabels{ity}).label = param.ica.alllabels(num);
                end
                param.replot = true;
            end
    end
else
    param.currentc = icomp;
end

dlg.UserData = param;
uiresume;
end

function xlog(~, eventdata,axespow)
if eventdata.Source.Value;
    set(axespow,'XScale','log')
else
    set(axespow,'XScale','linear');
end
end


function foichange(~, eventdata,axespow)
targetf = str2double(eventdata.Source.String);
type = eventdata.Source.Tag;

currentlim = get(axespow,'XLim');
switch type
    case 'fmin'
        currentlim(1) = targetf;
    case 'fmax'
        currentlim(2) = targetf;
end
if any(isnan(currentlim)) || currentlim(2)<=currentlim(1)
    fprintf('Cannot use the inserted frequency limits\n')
    
else
    set(axespow,'XLim',currentlim);
end

end