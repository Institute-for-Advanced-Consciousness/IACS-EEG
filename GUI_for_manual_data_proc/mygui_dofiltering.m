function trialdata = mygui_dofiltering(cfgdat,dataft)

if isfield(cfgdat.raw.visfilt,'hpfreq')
    hpfreq = cfgdat.raw.visfilt.hpfreq;
else
    hpfreq =[];
end
if isfield(cfgdat.raw.visfilt,'lpfreq')
    lpfreq = cfgdat.raw.visfilt.lpfreq;
else
    lpfreq =[];
end

% Default filtering parameter
N = 6; %order
Fn = dataft.fsample/2;


trialdata = cell(size(dataft.trial));
for itrial = 1:length(dataft.trial)
    auxtrial = double(dataft.trial{itrial});
    % demean
    meandat = mean(auxtrial,2);
    auxtrial = bsxfun(@minus, auxtrial, meandat);
    if ~isempty(hpfreq)
        [B, A] = butter(N, hpfreq/Fn, 'high');
        auxtrial = filtfilt(B, A, auxtrial')';
    end
    if ~isempty(lpfreq)
        [B, A] = butter(N, lpfreq/Fn);
        auxtrial = filtfilt(B, A, auxtrial')';
    end
    trialdata{itrial}=auxtrial;
end

