function [savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES_ori, SummaryStruct, ...
    prms, ListOfTimeWindows, ListOfFrBinsizes, savenotes)

prms.ClassGeneral.Nmin = 6; %
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

%%

tstamp = lt_get_timestamp(0);
% savedir = [savedir '/Results_' prms.Extract.strtype '_' tstamp '_' savenotes];
savedir = [savedir '/Results_' prms.Extract.strtype '_AlgnSyl' num2str(prms.alignWhichSyl) 'Onset' ...
    num2str(prms.alignOnset) '_' tstamp '_' savenotes];
mkdir(savedir);

%%
cd(savedir)

save('SummaryStruct', 'SummaryStruct');

%% lt 8/12/17 - will run classifier for each combination of time window and bin size

counter = 1;
for i=1:size(ListOfTimeWindows,1);
    
    prms.ClassGeneral.frtimewindow =ListOfTimeWindows(i,:); % on and off, relative to syl onset
    
    for ii=1:length(ListOfFrBinsizes)
        prms.ClassGeneral.frbinsize = ListOfFrBinsizes(ii); % in s.
       
        % === run
        CLASSES = lt_neural_v2_CTXT_ClassGeneral(CLASSES_ori, SummaryStruct, prms);
        
        % ======== save output
        classes_fname = ['classes' num2str(counter) '.mat'];
        prms_fname = ['params' num2str(counter) '.mat'];
        
        save(classes_fname, 'CLASSES');
        save(prms_fname, 'prms');

        % - increment
        counter = counter+1;
        
    end
end



