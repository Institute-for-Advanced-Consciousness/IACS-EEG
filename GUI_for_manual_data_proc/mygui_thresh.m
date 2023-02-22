function[] = mygui_thresh(h, ~,data)
cfg=guidata(h);

hasgamma = sum(ismember(data.label,'GAMMA'));

switch hasgamma
    case true
        ndx = find(ismember(data.label,'GAMMA'));
        logpow = log10(data.trial{1}(ndx,:).^2);
    case false
        cfg.dattype.muscle = [];
        % FIX THIS BELOW SO THERE ARE NO NANS IN THE 
        gamma = data.dat_hp;
        gamma = zscore(gamma')'; % rescale after removing NaNs
        env = nan(size(gamma));
        logpow = nan(size(gamma));
        fnyq = data.fsample/2; % nyquest frequency
        n = fnyq + 1; % make odd number
        N = (n-1)/2;

        % loop through rows b/c hilbert function bad w/2D data
        for i = 1:size(gamma)
            H = hilbert(gamma(i,:));
            amp = sqrt(real(H).^2+imag(H).^2);
            %[amp,~] = envelope(gamma(i,:); % other function that can be used
            smoothEnv = conv(amp,ones(1,n)./n); % smooth the envelope squared (inst. power)
            env(i,:) = smoothEnv(N+1:end-N); % make same size as signal 
            logpow(i,:) = log10(env(i,:).^2); % get log of inst. power
        end
end

[N,X] = hist(mean(logpow,1),1000); 
figure
plot(X,N)
xlabel('log_{10} instanteous gamma power')
ylabel('count')
title('High frequency power distribution ([1Hz - machine LP] - [1 Hz - 32 Hz])')
grid on
display('Pausing ... pressing any key to continue and set threshold.')
pause
thresh = input('Muscle threshold?  ');

isMuscle = (mean(logpow,1) > thresh);
close
% problem with the diff function when very beginning and end are marked

% below appened 0 to beginning and end for taking the derivative 
isMuscle = [zeros(size(isMuscle,1),1) isMuscle zeros(size(isMuscle,1),1)];
d = diff(isMuscle); % find start/stop boundaries
start = find(d==1);
stop = find(d==-1);
gap = [start 0]- [0 stop]; % gap between consecutive blocks
where = gap < data.fsample/2; % close gaps less than half second
where(end) = 0; where(1) = 0; % fix this, else end will always be marked and sizes won't work
start(where(1:end-1)) = []; stop(where(2:end)) = [];
duration = stop-start; % find the number of samples in each block
reject = duration < data.fsample/5; % anything less than a fifth of a second discard 
start(reject) = []; stop(reject) = []; 
cfg.dattype.muscle = [start' stop'];
%cfg.redraw.main = 1; % set redraw as true
guidata(h,cfg);
uiresume
end