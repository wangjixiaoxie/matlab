function DATAllSwitches = lt_neural_BatchSmthCOMPILE(basedir, ListOfDatstructs);
%% lt 11/17/17 - given multiple switches, compiles into one structure


% basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/';
% ListOfDatstructs = {...
%     'DATSTRUCT_BatchSm_17Nov2017_1235.mat',...
%     'DATSTRUCT_BatchSm_17Nov2017_1227.mat'}; % needs to be in basedir.
%     Should enter these in chrono order.


%%

DATAllSwitches = struct;
for i=1:length(ListOfDatstructs)
    
    datstruct = load([basedir ListOfDatstructs{i}]);
    dirfields = fieldnames(datstruct.DATSTRUCT);
    
    DATAllSwitches.switch(i).datstructname = ListOfDatstructs{i};
    DATAllSwitches.switch(i).basedir = basedir;
    DATAllSwitches.switch(i).params = datstruct.DATSTRUCT.(dirfields{1}).params;
    
    % ---- figure out order of dirfields (chronoligcla)
    tmins = []; % first and last song in each field
    tmaxs = [];
    for dd = 1:length(dirfields)
        t = datstruct.DATSTRUCT.(dirfields{dd}).motifnum(1).t;
        
        tmins = [tmins min(t)];
        tmaxs = [tmaxs max(t)];
    end
    [~, inds1] = sort(tmins);
    [~, inds2] = sort(tmaxs);
    assert(all(inds1 == inds2), 'mixed ordering ...');
    dirfields = dirfields(inds1); % put in correct order!
    
    
    % =============== for each channel and motif, plot post and pre
    nummotifs = length(datstruct.DATSTRUCT.(dirfields{1}).motifnum);
    allchans = find(~cellfun('isempty', datstruct.DATSTRUCT.(dirfields{1}).motifnum(1).DatAll));
    
    for mm = 1:nummotifs
        
        %         for cc = allchans
        %
        %             for dd =1:length(dirfields)
        %             Datsm = datstruct.DATSTRUCT.(dirfields{dd}).motifnum(mm).DatAll{cc};
        %             datstruct.DATSTRUCT.(dirfields{dd}).motifnum(mm)
        %
        %             % ========= out
        %             DATAllSwitches.switch(i).motif(mm).chan(cc).batchinorder(dd).condition = dirfields{dd};
        %             DATAllSwitches.switch(i).motif(mm).chan(cc).batchinorder(dd).Datsm = Datsm;
        %
        %
        %             DATAllSwitches.switch(i).motif(mm).chan(cc).batchinorder(dd) = ...
        %                 datstruct.DATSTRUCT.(dirfields{dd}).motifnum(mm);
        %             end
        %         end
        motif = datstruct.DATSTRUCT.(dirfields{1}).motifnum(mm).motifname;
        
        DATAllSwitches.switch(i).motif(mm).motif = motif;
        DATAllSwitches.switch(i).motif(mm).batchinorder = [];
        
        for dd =1:length(dirfields)
            
            datstruct.DATSTRUCT.(dirfields{dd}).motifnum(mm).condition = dirfields{dd};
            
            
            % ========= out
            DATAllSwitches.switch(i).motif(mm).batchinorder = ...
                [DATAllSwitches.switch(i).motif(mm).batchinorder ...
                datstruct.DATSTRUCT.(dirfields{dd}).motifnum(mm)];
        end
        
    end
end


% ====================== SAVE
fname = [basedir 'DATSTRUCT_BatchSmCombined.mat'];
save(fname, 'DATAllSwitches', '-v7.3');
