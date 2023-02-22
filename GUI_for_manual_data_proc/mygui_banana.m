function[trialdata,newLabels] = mygui_banana(cfgdat,data)

% A = {'FP1','F3','C3','P3','FP2','F4','C4','P4','FP1','F7','T3','T5','FP2',...
%     'F8','T4','T6','FZ','CZ'};
% B = {'F3','C3','P3','O1','F4','C4','P4','O2','F7','T3','T5','O1','F8','T4', ...
%     'T6','O2','CZ','PZ'};

A = {'FP1','F3','C3','P3','FP1','F7','T3','T5','FP2','F4','C4','P4','FP2',...
    'F8','T4','T6','FZ','CZ'};


B = {'F3','C3','P3','O1','F7','T3','T5','O1','F4','C4','P4','O2','F8','T4', ...
    'T6','O2','CZ','PZ'};


% % Combine channel labels until new list
% C = cell(1,length(A));
% for ich = 1:length(A)
%     C{ich} = sprintf('%s-%s',A{ich},B{ich});
% end
    
newLabels = strcat(A,'-',B); % cell array of new labels 
howMany = length(cfgdat.raw.mainplot.allchid);


for itrial = 1:length(data.trial)
    nPhys = size(data.trial{itrial},1) - howMany; % number of physiological channels to keep
    N = length(A) + nPhys;
    auxtrial = zeros(N,length(data.trial{itrial}));
    cfgdat.raw.alllabels
    for ich = 1:length(A)
        I = find(strcmpi(data.label,A{ich}));
        J = find(strcmpi(data.label,B{ich}));
        auxtrial(ich,:) = data.trial{itrial}(I,:) - data.trial{itrial}(J,:); % bipolar reference
    end
    trialdata{itrial}=auxtrial;
end


end


