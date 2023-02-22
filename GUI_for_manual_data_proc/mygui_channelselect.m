function chinf = mygui_channelselect(chinf)


%% Set parameters


ctypeothlabel = chinf.ctypelabels( setdiff(1:length(chinf.ctypelabels),chinf.ctypemainid));
ctypegoodlabel = chinf.ctypemainlabel;

if isempty(ctypeothlabel)
   fprintf('Cannot change channel/component selection - there is a single type = %s\n',chinf.ctypemainlabel);
   return
end

%% Create figure


pos      = get(0,'DefaultFigurePosition');
pos(3:4) = [300 600]; %pos(3:4) = [400 800];
dlg      = dialog('Name', 'Channel selection', 'Position', pos);
% dlg = figure('Position', pos);
set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears

%% dimensions


uicontrol(dlg, 'style', 'pushbutton', 'position', [ 55  10    80  20], 'string', 'Cancel',       'callback', 'close');
uicontrol(dlg, 'style', 'pushbutton', 'position', [155  10    80  20], 'string', 'OK',           'callback', 'uiresume');


dim.hborderspace = 20;
dim.lowerregionheight = 40; 
dim.chboxwidth = round((pos(3) - 2*dim.hborderspace )/2.5);

dim.othchheight = round(0.85*(pos(4)-dim.lowerregionheight )/length(ctypeothlabel));
dim.sepheight = round(0.15*(pos(4)-dim.lowerregionheight )/length(ctypeothlabel));

dim.butwidth = round(0.3*dim.chboxwidth);
dim.butheight = 30;

dim.goodchheight = length(ctypeothlabel)*(dim.othchheight+dim.sepheight) -  dim.sepheight;

%% create ui buttons and list

% list for other types channels
for it = 1 : length(ctypeothlabel)
    auxpos = [dim.hborderspace, dim.lowerregionheight+(it-1)*(dim.othchheight+dim.sepheight), ...
        dim.chboxwidth, dim.othchheight];
    chinf.guilist.(ctypeothlabel{it}) = uicontrol(dlg,'units','pixels','position', ...
        auxpos,'Style', 'listbox', 'String',  ...
        chinf.ctype.(ctypeothlabel{it}).label,'min', 0, 'max', 2);
    
    auxpos = [auxpos(1), auxpos(2)+dim.othchheight, dim.chboxwidth, dim.sepheight];
    uicontrol(dlg,'units','pixels','position', ...
        auxpos,'Style', 'text', 'String', (ctypeothlabel{it}));
    
end



% List for type good channel
auxpos = [pos(3)- dim.hborderspace - dim.chboxwidth, dim.lowerregionheight, ...
    dim.chboxwidth, dim.goodchheight];
chinf.guilist.(ctypegoodlabel) = uicontrol(dlg,'units','pixels','position', ...
    auxpos,'Style', 'listbox', 'String', chinf.ctype.(chinf.ctypemainlabel).label,'min', 0, 'max', 2);
    auxpos = [auxpos(1), auxpos(2)+dim.goodchheight, dim.chboxwidth, dim.sepheight];
uicontrol(dlg,'units','pixels','position', ...
    auxpos,'Style', 'text', 'String',(chinf.ctypemainlabel));




% % buttons to switch from one side to the other
for it = 1 : length(ctypeothlabel)
    dim.buttonx = round(( dim.hborderspace + pos(3)- dim.hborderspace - dim.butwidth)/2);
    dim.buttony = dim.lowerregionheight+(it-1)*(dim.othchheight+dim.sepheight) + ...
        round( dim.othchheight/2 )- dim.butheight;
    
    uicontrol(dlg,'units','pixels','position', ...
        [dim.buttonx dim.buttony dim.butwidth dim.butheight],'Style', 'pushbutton', 'String', '<',...
        'Tag',(ctypeothlabel{it}),'Callback', @removegoodch);
    uicontrol(dlg,'units','pixels','position', ...
        [dim.buttonx dim.buttony+dim.butheight dim.butwidth dim.butheight], ...
        'Style', 'pushbutton', 'String', '>','Tag',(ctypeothlabel{it}),'Callback', @addgoodch);
end



%% SAVES DATA
set(dlg, 'userdata', chinf); 
uiwait(dlg);

if ishandle(dlg)
    for itype = 1:length(chinf.ctypelabels)
        auxlabel = chinf.ctypelabels{itype};
        chinf.ctype.(auxlabel).num = find(ismember(chinf.alllabels, ...
            chinf.guilist.(auxlabel).String));
        chinf.ctype.(auxlabel).label = chinf.alllabels(chinf.ctype.(auxlabel).num);
    end
    chinf = rmfield(chinf,'guilist');
    
    if ~isfield(chinf,'numofic') % if eeg
        chinf.mainplot.allchid = chinf.ctype.(chinf.ctypemainlabel).num(:);
        chinf.mainplot.allchcolor = lines(length(chinf.mainplot.allchid));
        if chinf.badchplot
            auxnum = chinf.ctype.bad.num(:);
            chinf.mainplot.allchid = [chinf.mainplot.allchid; auxnum];
            chinf.mainplot.allchcolor = [chinf.mainplot.allchcolor; chinf.badchcolor(ones(1,length(auxnum)),:)];
        end
        [chinf.mainplot.allchid, goodorder] = sort(chinf.mainplot.allchid,'ascend');
        chinf.mainplot.allchcolor = chinf.mainplot.allchcolor(goodorder,:);
        chinf.mainplot.allchlabel = chinf.alllabels(chinf.mainplot.allchid);
        oldbegch = min(chinf.mainplot.showchid);
        begch = min([oldbegch, length(chinf.mainplot.allchid)-chinf.nsigblock+1]);
        begch = max([begch 1]);
        endch = min([begch + chinf.nsigblock-1, length(chinf.mainplot.allchid)]);
        chinf.mainplot.showchid = endch:-1:begch;
    else % if selecting ICA components all are plotted so no need to update chinf.mainplot
        % color gray all component types that are not main type
        idbad = setdiff(1:chinf.numofic,chinf.ctype.(chinf.ctypemainlabel).num);
        chinf.mainplot.allchcolor(idbad,:) = repmat(0.6*[1 1 1],length(idbad),1);
    end
    
    
    
    
    close(dlg);
end
end


function removegoodch(h, eventdata)
h = get(h, 'parent');
chinf = get(h, 'userdata');

% find to which category we should add the channel
addto = eventdata.Source.Tag;

if ~isempty(chinf.guilist.(chinf.ctypemainlabel).String)
% find which labels were selected in the goodchannel box
chlabel =  chinf.ctype.(chinf.ctypemainlabel).label(chinf.guilist.(chinf.ctypemainlabel).Value);

idchkeep = setdiff(1:length(chinf.ctype.(chinf.ctypemainlabel).label), chinf.guilist.(chinf.ctypemainlabel).Value);
chinf.ctype.(chinf.ctypemainlabel).label = chinf.ctype.(chinf.ctypemainlabel).label(idchkeep);
chinf.guilist.(chinf.ctypemainlabel).Value = [];
chinf.guilist.(chinf.ctypemainlabel).String = chinf.ctype.(chinf.ctypemainlabel).label;

chinf.ctype.(addto).label(end+1:end+length(chlabel)) = chlabel;
chinf.guilist.(addto).String = chinf.ctype.(addto).label;

set(h, 'userdata', chinf);
end
end


function addgoodch(h, eventdata)
h = get(h, 'parent');
chinf = get(h, 'userdata');

% find to which category we should removethe 
addfrom = eventdata.Source.Tag;

% find which labels were selected in the goodchannel box
if ~isempty(chinf.guilist.(addfrom).String)
chlabel =  chinf.ctype.(addfrom).label(chinf.guilist.(addfrom).Value);

idchkeep = setdiff(1:length(chinf.guilist.(addfrom).String), chinf.guilist.(addfrom).Value);

chinf.ctype.(addfrom).label = chinf.ctype.(addfrom).label(idchkeep);
chinf.guilist.(addfrom).Value = [];
chinf.guilist.(addfrom).String = chinf.ctype.(addfrom).label;

chinf.ctype.(chinf.ctypemainlabel).label(end+1:end+length(chlabel)) = chlabel;
chinf.guilist.(chinf.ctypemainlabel).String = chinf.ctype.(chinf.ctypemainlabel).label;

set(h, 'userdata', chinf);
end
end