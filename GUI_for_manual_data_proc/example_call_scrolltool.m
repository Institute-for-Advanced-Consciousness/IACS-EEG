%% paths
% eeglab_dir='/homebasel/clinbiomarkers/garceslm/projects/toolboxes/eeglab13_5_4b/';
eeglab_dir='C:/Users/garceslm/Desktop/toolboxes/eeglab13_5_4b/';
addpath(eeglab_dir);
addpath([eeglab_dir 'functions/popfunc/'])
addpath([eeglab_dir 'functions/guifunc'])
addpath([eeglab_dir 'functions/adminfunc'])
addpath([eeglab_dir 'functions/sigprocfunc'])

% addpath('/homebasel/clinbiomarkers/garceslm/projects/toolboxes/fieldtrip-20160112/');
addpath('C:\Users\garceslm\Desktop\toolboxes\fieldtrip-20160215/');
ft_defaults;

addpath('C:\Users\garceslm\Desktop\toolboxes\fieldtrip-20160215\external\fastica');



%%
load('C:\myEEGdata\scripts\preprocessing/testingsubjects.mat');
subjectcode = testingsubjects{12};
% 1 - 6  / 12
% subjectcode = '137282891018'; % smaller
% subjectcode = '327288335476'; % big ocular and muscle
% 5 difficult one
% 6 also interesting

% subject 3 lots of bad components

cfg = [];


cfg.datafile = ['C:\myEEGdata\preproc\' subjectcode ' - RESTING_STATE.set'];
cfg.trialfile =[ 'C:/myEEGdata/RestingPreproc/trlorig/trialdef-' subjectcode '.mat'];
cfg.preprocfun = @mypreproc_euaims;

cfg.fileout = ['C:\myEEGdata\RestingPreproc\preproc/' subjectcode 'Restingpreproc.mat'];



cfg.raw.ctypelabels ={'good','bad','discard','physio'};
cfg.raw.ctypepatterns.discard = {'FT10','FT9','TP9','TP10','P9','P10','PO9','PO10','Iz',...
    'EOGL','EOGR','EOGA','EOGB','ML','MR'};
cfg.raw.ctypepatterns.physio = {'EOGHdiff'};
cfg.raw.badchplot = false; %true;

cfg.art.types = {'ocular','muscle','bigmusclemov','other'};
cfg.ica.discardarttype = {'bigmusclemov','other'};
cfg.ica.numofic = 35;
cfg.ica.ctypelabels = {'good','muscle','ocul','lightmuscle','other'};

cfg.raw.plotavgref = true; %average re-referencing for plotting only

scroll_data(cfg);




