function data = myload_euaims(cfg)



load(cfg.filename);

trlinf = data.cfg.trl;

fprintf(' ********************* CHANNEL INFO ***********************\n')
nsamplesT=sum(trlinf.trl(:,2)-trlinf.trl(:,1));


sumsmalldiff = nan(length(data.label),length(data.trial));
for itrial=1:size(trlinf.trl,1);
    gadiff = median(median(abs(diff(data.trial{itrial},2)),2 ));
    for ich=1:length(data.label) 
        sumsmalldiff(ich,itrial)= sum( abs(diff(data.trial{itrial}(ich,:)) ) < 0.001*gadiff );
    end
end
propsmalldiff=(sum(sumsmalldiff,2))/nsamplesT;
for ich=1:length(data.label)
    if propsmalldiff(ich) > 0.05 && ~strcmpi(data.label{ich},'FCz')
    fprintf('VERY SMALL FLUCTUATIONS: %s - in %i/%i trials - around %.1f percent time \n', ...
        data.label{ich},sum(sumsmalldiff(ich,:)>0),size(trlinf.trl,1), ...
        100*propsmalldiff(ich));
    end
end

for ich=1:length(data.label) %FCz has to be flat - no artifact
    if ~strcmpi(data.label{ich},'FCz')
        sumzero=nan(1,size(trlinf.trl,1));
        for itrial=1:size(trlinf.trl,1);
            sumzero(itrial)= sum( diff(data.trial{itrial}(ich,:)) ==0  );
        end
        propzero=(sum(sumzero))/nsamplesT;
        if propzero > 0.05
            fprintf('FLAT CHANNEL: %s - diff==0 in %i/%i trials - around %.1f percent time \n', ...
                data.label{ich},sum(sumzero>0),size(trlinf.trl,1), ...
                100*propzero);
        end
    end
end


%%
% options.dirinf='/homebasel/clinbiomarkers/garceslm/projects/euaims/scripts/subject_info/';
% load(fullfile(options.dirinf, 'euaimslog.mat'));
% idsublog=find(strcmpi({euaimslog.code},data.cfg.euaims.id));
% if ~isempty(idsublog)
%    fprintf('Log info:    Capping: %s - %s\n',euaimslog(idsublog).Capping,euaimslog(idsublog).CappingComments)
%    fprintf('Log info: Impedances: %s - %s\n',euaimslog(idsublog).Impedances,euaimslog(idsublog).ImpedanceComments)
% end

for it = 1:length(data.cfg.euaims.comments)
    fprintf('%s\n',data.cfg.euaims.comments{it});
end


% %% IMPEDANCES
% 
% if isfield(data.chanlocs,'impedance')
%    idhighimp = find([data.chanlocs.impedance]>50);
%    fprintf('Channel with >50kOhm impedance: %s \n', data.chanlocs(idhighimp).labels );
% end

%% CHANGE TO NEW FORMAT

if isfield(data,'cfg') && ~isfield(data.cfg,'dattype') && isfield(data.cfg,'art')
    arttypes = fieldnames(data.cfg.art);
    for iaty = 1:length(arttypes)
       data.cfg.dattype.(arttypes{iaty}) = data.cfg.art.(arttypes{iaty}).artifact; 
    end
end