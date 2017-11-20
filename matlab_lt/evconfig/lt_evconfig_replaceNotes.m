function [fname] = lt_evconfig_replaceNotes(config_targ, notenum_targ, config_source, ...
    notenum_source)

%% ============ transfer template params from one config to another
% LT 11/17/17

% Works only for evtafv4. 

% the template (and therefore number of c0lumns) must be identical in
% source and target configs.

% does not overwrite, but saves as name_v2

% config_source = '/bluejay5/lucas/birds/pu46wh04/config111817_onecol.evconfig2';
% notenum_source = 1; % 1, 2, ...
% 
% config_targ = '/bluejay5/lucas/birds/pu46wh04/config111817.evconfig2';
% notenum_targ = [1 3]; % array (1, 2. ...)

%% ================ run

NDsource = ReadEvTAFv4ConfigFile(config_source);
[NDtarg, OP] = ReadEvTAFv4ConfigFile(config_targ);

%% for each target note, replace params with that of the source

for i=1:length(notenum_targ)
   
    ind = notenum_targ(i);
    
    
    % ------------ 1) confirm that is same template file
    assert(strcmp(NDtarg(ind).TemplFile, NDsource(notenum_source).TemplFile), 'problem - diff templatews...');
        
    % ------------- 2) confirm that is same num columns
    assert(length(NDtarg(ind).CntRng) == length(NDsource(notenum_source).CntRng), 'proplem - diff num columns. ..');
    
    
    % ============================== replace target note with source note
    % --- target
    NDtarg(ind) = NDsource(notenum_source);
    disp( ' ')
    disp(['============ replaceing targ config notenuk: ' num2str(ind) ' with source config notenum '  num2str(notenum_source)']);
end


%% =============== save

indtmp = strfind(config_targ, '.evconfig2');
fname = [config_targ(1:indtmp-1) '_v2.evconfig2'];

disp(' ');
disp(['=============== Saving config file named : ' fname]);
WriteEvTAFv4ConfigFile(fname,NDtarg,OP);
