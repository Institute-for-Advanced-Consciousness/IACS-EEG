function [foi,powspctrm] = my_compute_powerspctrm(optfreq,Fs,alldat)

% INPUT:
% - optfreq. options for computer powerspctrum, with possible fields:
%       foilim = frequency limits of interest (default [1 150])
%       deltaf = smoothing width (Hz) (default= 0.5)
%       window_shift = window shift, in ratio of window (default = 0.5)
%       window_shift = window length to compute fft (default 2.^nextpow2(Fs/optfreq.deltaf) )
% - Fs = sampling frequency, Hz
% - alldat = input data: Nchannels x Nsamples

if ~isfield(optfreq,'deltaf')
    optfreq.deltaf = 0.5;
end
if ~isfield(optfreq,'window_shift')
    optfreq.window_shift = 0.5;
end
if ~isfield(optfreq,'foilim')
    optfreq.foilim = [1 150];
end
if ~isfield(optfreq,'ntotwindow')
    optfreq.ntotwindow = 2.^nextpow2(Fs/optfreq.deltaf);
end

if optfreq.foilim(1)<5/Fs
    optfreq.foilim(1) = 5/Fs;
end

if optfreq.foilim(2)>0.3*Fs
    optfreq.foilim(2) = 0.3*Fs;
end

% Kernel
T = 1/optfreq.deltaf;       % [s]
n_win   = round(T * Fs);
n_shift = round(n_win * optfreq.window_shift);
KERNEL  = hanning(n_win); KERNEL=KERNEL/sum(KERNEL);


% select windows for fft
possiblewindowsbeg = 1: n_shift : size(alldat,2)-n_win;
possiblewindows = [possiblewindowsbeg', possiblewindowsbeg'+n_win-1];
windowok=true(size(possiblewindows,1),1);
for iwin=1:size(possiblewindows,1)
    if any( isnan(alldat(1,possiblewindows(iwin,1):possiblewindows(iwin,2))) )
        windowok(iwin)=false;
    end
end
windows = possiblewindows(windowok,:);


% Compute fft
idfreqlim = round(optfreq.foilim .* optfreq.ntotwindow / Fs) + 1;
idfreq    = idfreqlim(1):1:idfreqlim(2);
foi = (idfreq-1) ./ (optfreq.ntotwindow / Fs);
powspctrmaux = nan(length(foi),size(alldat,1),size(windows,1)); %frec x channel x window
for iwin = 1:size(windows,1)
    section = double(alldat(:,windows(iwin,1):windows(iwin,2)));
    aux = fft(section'.*repmat(KERNEL,[1,size(section,1)]),optfreq.ntotwindow);
    powspctrmaux(:,:,iwin) = abs(aux(idfreq,:)).^2/ (Fs / optfreq.ntotwindow);
end
powspctrm = trimmean(powspctrmaux,10,'round',3); %average over windows