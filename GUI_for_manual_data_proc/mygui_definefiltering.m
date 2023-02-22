function mygui_definefiltering(h,eventdata)
dat=guidata(h);


% Get current filtering parameters

if isfield(dat.raw.visfilt,'hpfreq')
    param.hp.freq = dat.raw.visfilt.hpfreq;
else
    param.hp.freq = [];
end
if isfield(dat.raw.visfilt,'lpfreq')
    param.lp.freq = dat.raw.visfilt.lpfreq;
else
    param.lp.freq = [];
end
param.Fs = dat.fsample;

%% Create popup window

pos      = get(0,'DefaultFigurePosition');
pos(3:4) = [300 400];
% dlg      = dialog('Name', 'Channel selection', 'Position', pos);
dlg = figure('Position', pos);
set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears

%% Create uicontrol elements in the window
uicontrol(dlg, 'style', 'pushbutton', 'position', [ 55  10    80  20], 'string', 'Cancel',       'callback', 'close');
uicontrol(dlg, 'style', 'pushbutton', 'position', [155  10    80  20], 'string', 'OK',           'callback', 'uiresume');

uicontrol(dlg,'units','pixels','position', [20 100 100 30],'Style', 'text', 'String', 'high pass frequency (Hz)');
param.hp.box = uicontrol(dlg,'units','pixels','position',[120 100 50 30], ...
    'Style','edit','Tag','hp','String',num2str(param.hp.freq),'Callback',@freqnumber,'HorizontalAlignment','center');

uicontrol(dlg,'units','pixels','position', [20 200 100 30],'Style', 'text', 'String', 'low pass frequency (Hz)');
param.lp.box = uicontrol(dlg,'units','pixels','position',[120 200 50 30], ...
    'Style','edit','Tag','lp','String',num2str(param.lp.freq),'Callback',@freqnumber,'HorizontalAlignment','center');




%%
dlg.UserData = param; 
uiwait(dlg);
param = dlg.UserData;
if ishandle(dlg)
    dat.raw.visfilt.filter = true;
    dat.raw.visfilt.lpfreq = param.lp.freq;
    dat.raw.visfilt.hpfreq = param.hp.freq;
    close(dlg);
end

guidata(h,dat);
uiresume;
end


function freqnumber(dlg, eventdata)
dlg = get(dlg, 'parent');
param = get(dlg, 'userdata');
targetnumber = str2double(eventdata.Source.String);
tag = eventdata.Source.Tag;
if isempty(eventdata.Source.String)
    param.(tag).freq = [];
elseif isnan(targetnumber) || targetnumber<0 || targetnumber>=param.Fs/2;
    fprintf('cannot use this frequency threshold\n')
    param.(tag).box.String = num2str(param.(tag).freq);   
else
    param.(tag).freq = str2double(param.(tag).box.String);
end
set(dlg, 'userdata', param);
end