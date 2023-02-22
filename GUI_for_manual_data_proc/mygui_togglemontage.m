function mygui_togglemontage(h,~)
dat=guidata(h);

if isempty(dat.raw.plotbanana) || ~dat.raw.plotbanana
    dat.raw.plotbanana = true;
else
    dat.raw.plotbanana = false;
end

guidata(h,dat);
uiresume;
end
