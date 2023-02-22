
function mygui_cplots(h, ~,data)
cfg=guidata(h);

% auxobj = findobj(h,'Enable','on');
% set(auxobj,'Enable','off')


if isempty(cfg.(cfg.disp.disptype).mainplot.chhighlight)
    fprintf('No channel/component has been highlighted\n')
else
    % Select samples without artifacts that are excluded in ICA
    hasartifact = false(1,max(cfg.trl.trl(:)));
    for iaty = cfg.ica.discarddattypeid
        aux =  cfg.dattype.(cfg.dattype.labels{iaty});
        for iart = 1:size(aux,1)
            hasartifact(aux(iart,1):aux(iart,2)) = true;
        end
    end
    idc = cfg.(cfg.disp.disptype).mainplot.allchid(cfg.(cfg.disp.disptype).mainplot.chhighlight);
    clab = cfg.(cfg.disp.disptype).alllabels{idc};
    switch cfg.disp.disptype %return topoint timeseries with nans.
        case 'raw'
            topoint = zeros(length(cfg.raw.alllabels),1);
            topoint(idc) = 1;
            topoint = topoint(cfg.raw.haspos);
            interpmethod = 'nearest';
            cpos = cfg.raw.allpos(cfg.raw.haspos,:);
            clim=[0 1];
            
            alldat = nan(1,max(cfg.trl.trl(:)));
            for itrial = 1:length(data.trial)
                alldat(cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2)) = ...
                    data.trial{itrial}(idc,:);
            end
            
        case 'ica'
            
            topoint = cfg.ica.a(:,idc);
            cpos = cfg.raw.allpos(cfg.ica.userawch,:);
            idok = cfg.raw.haspos(cfg.ica.userawch);
            topoint = topoint(idok); cpos = cpos(idok,:);
            interpmethod = 'v4';
            clim = max(abs(topoint))*[-1 1];
            
            
            alldat = nan(1,max(cfg.trl.trl(:)));
            for itrial = 1:length(data.trial)
                alldat(cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2)) = ...
                    cfg.ica.w(idc,:)*data.trial{itrial}(cfg.ica.userawch,:);
            end
            
            
            
    end
    alldat(hasartifact) = nan;
    
    %% CREATE FIGURE
    
    pos      = get(0,'DefaultFigurePosition');
    pos(3:4) = [800 400];
    %dlg      = dialog('Name', 'Plot channel/component', 'Position', pos);
    % dlg = figure('Position', pos,'WindowStyle','modal','MenuBar','figure','ToolBar','figure');
    dlg = figure('Position', pos);

    
    
    
    %set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears
    
    axestopo  = axes('Position',[0.0 0.2 0.45 0.7]);
    axespow   = axes('Position',[0.5 0.3 0.45 0.6]); axis on
    %     axestopo.Position(
    
    %% PLOT TOPOPLOT
    
    subplot(axestopo); cla; hold on
    if size(cpos,1)>1
    mytopoplot(topoint,cpos,interpmethod,[],100)
    set(axestopo,'CLim',clim);
    colorbar off;
    colormap('jet');
    end

    
    
    %% POWER SPECTRUM DEFAULT PARAMETERS
    
    Fs=cfg.fsample;
    
    % options
    optfreq.deltaf = 0.5; %frequency resolution, Hz
    optfreq.window_shift = 0.5; % fraction of the window to shift
    optfreq.fmin = 2;
    optfreq.fmax = 45;
    optfreq.Fs = Fs;
    optfreq.lengths = length(alldat);
    
    
    optfreq.xlog = false;
    uicontrol(dlg,'units','normalized','position',[0.8  0.03  0.08  0.05], ...
        'Style', 'togglebutton', 'String', 'xlog','Value',optfreq.xlog,'callback', @xlog);
    
    uicontrol(dlg,'units','normalized','position',[0.5  0.08  0.08  0.05], ...
        'Style', 'text', 'String', 'fmin');
    optfreq.box.fmin = uicontrol(dlg,'units','normalized','position',[0.5  0.03  0.08  0.05], ...
        'Style', 'edit', 'tag','fmin','String', num2str(optfreq.fmin),'callback', @foichange);
    
    uicontrol(dlg,'units','normalized','position',[0.6  0.08  0.08  0.05], ...
        'Style', 'text', 'String', 'fmax');
    optfreq.box.fmax = uicontrol(dlg,'units','normalized','position',[0.6  0.03  0.08  0.05], ...
        'Style', 'edit','tag','fmax', 'String',  num2str(optfreq.fmax),'callback', @foichange);
    
    uicontrol(dlg,'units','normalized','position',[0.68  0.08  0.12  0.05], ...
        'Style', 'text', 'String', 'freq smooth (Hz)');
    optfreq.box.deltaf =uicontrol(dlg,'units','normalized','position',[0.7  0.03  0.08  0.05], ...
        'Style', 'edit', 'String', num2str(optfreq.deltaf),'callback', @deltafchange);
    
    optfreq.extraicaplots = false;
    if strcmpi(cfg.disp.disptype,'ica')
        optfreq.box.moreicaplots = uicontrol(dlg,'units','normalized','position',[0.90  0.03  0.08  0.05], ...
            'Style', 'togglebutton', 'String', 'More ICA plots','Value',optfreq.extraicaplots,'callback', 'uiresume');
    else
        optfreq.box.moreicaplots.Value = false;
    end
    
    if strcmpi(cfg.disp.disptype,'raw') && isfield(cfg.ica,'a')
        idc2 = cfg.ica.userawch == idc;
        aux = cfg.ica.ctype.(cfg.ica.ctypemainlabel).num;
        [~,imaxmain] = max(abs(cfg.ica.a(idc2,aux)));
        imaxmain = aux(imaxmain);
        [~,imax] = max(abs(cfg.ica.a(idc2,:)));
        uicontrol(dlg,'units','normalized','position',[0.05  0.03  0.40  0.10], ...
        'Style', 'text', 'String', sprintf('A val max for IC %i (from unselec.comp %i)',imax,imaxmain));
    end
    
    
    %% COMPUTE POWER SPECTRUM
    
    dlg.UserData = optfreq;
    while(ishandle(dlg))
        optfreq = dlg.UserData;
        optfreq.foilim=[optfreq.fmin optfreq.fmax];
        optfreq.ntotwindow = 2.^nextpow2(Fs/optfreq.deltaf);
        
        [foi,powspctrm] = my_compute_powerspctrm(optfreq,Fs,alldat);
        
        subplot(axespow); cla; hold on
        plot(foi,powspctrm,'b','LineWidth',2);
        xlabel('frequency(Hz)')
        ylabel('Power spectrum');
        title(clab)
        if optfreq.xlog
            set(axespow,'XScale','log');
        else
            set(axespow,'XScale','linear');
        end
        xlim([min(foi) max(foi)])
        
        if optfreq.box.moreicaplots.Value
            alldat2 = nan(length(cfg.ica.userawch),max(cfg.trl.trl(:,2)));
            for itrial = 1:length(data.trial)
                alldat2(:,cfg.trl.trl(itrial,1):cfg.trl.trl(itrial,2)) = ...
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
            alldat2(:,~usesample) = nan;
            [~,rawpowspctrm] = my_compute_powerspctrm(optfreq,Fs,alldat2);
            
            powraw =  rawpowspctrm * abs(cfg.ica.a(:,idc));
            [~,powspctrmwithout] = my_compute_powerspctrm(optfreq,Fs,...
                alldat2 - cfg.ica.a(:,idc)*cfg.ica.w(idc,:)*alldat2);
            powrawcorr = powspctrmwithout * abs(cfg.ica.a(:,idc));
            
            
            powdiff = 1-powrawcorr./powraw;
            powraw = max(powspctrm)*powraw/max(powraw);
            
            [aux,p1,p2] = plotyy(foi,powraw,foi,powdiff);
            set(aux,{'ycolor'},{'k';[0 0.8 0]});
            set(p1,'color','r'); set(p2,'color',[0 0.8 0]);
            %plot(foi,powraw,'r');
            %yyaxis(axespow,'right');
            %plot(foi,powdiff);
            %plot([min(foi) max(foi)],[0 0],'Color',0.7*[1 1 1])
            %yyaxis(axespow,'left');
            % legend('IC power','raw power weighted by a','power decrease (%) weighted by a')
            xlim([min(foi) max(foi)])
        end

        uiwait(dlg);

    end
    %%
    
    
end
uiresume;
end



function xlog(dlg, eventdata)
dlg = get(dlg, 'parent');
optfreq = get(dlg, 'userdata');
optfreq.xlog = eventdata.Source.Value;
set(dlg, 'userdata', optfreq);
uiresume
end


function foichange(dlg, eventdata)
dlg = get(dlg, 'parent');
optfreq = get(dlg, 'userdata');
targetf = str2double(eventdata.Source.String);
type = eventdata.Source.Tag;
if isnan(targetf) || targetf <=0 || targetf >= optfreq.Fs/2 %revert to previous value
    fprintf('Cannot use the inserted frequency limits\n')
    optfreq.box.(type).String = num2str(optfreq.(type));
else
    switch type
        case 'fmin'
            if targetf >= optfreq.fmax
                targetf = optfreq.fmax-1;
                optfreq.box.fmin.String = num2str(targetf);
            end
            
        case 'fmax'
            if targetf <= optfreq.fmin
                targetf = optfreq.fmin+1;
                optfreq.box.fmax.String = num2str(targetf);
            end
    end
    optfreq.(type) = targetf;
    
end
set(dlg, 'userdata', optfreq);
uiresume
end


function deltafchange(dlg, eventdata)
dlg = get(dlg, 'parent');
optfreq = get(dlg, 'userdata');
targetdf = str2double(eventdata.Source.String);
if isnan(targetdf) %revert to previous value
    fprintf('Cannot use the inserted frequency limits\n')
    optfreq.box.deltaf.String = num2str(optfreq.deltaf);
else
    optfreq.deltaf = targetdf;
end
set(dlg, 'userdata', optfreq);
uiresume
end