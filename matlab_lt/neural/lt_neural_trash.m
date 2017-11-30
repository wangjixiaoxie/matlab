clear all; close all;
load CLASSESv2_xaa_Algn2Ons1_27Oct2017_1114_XLMAN25ms

%%

savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M/xaa_Algn2Ons1_27Oct2017_1114_XLMAN25ms/SHUFFDECODE';



%%

numbirds = length(CLASSES.birds);

for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    
    for ii=1:numneurons
        
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for j=1:numbranches
            
            % ---------------
            savefname = [savedir '/bird' num2str(i) '_neur' num2str(ii) ...
                '_branch' num2str(j) '_tbin' num2str(1) '.mat'];
            
            % ---------------- prepare file to save
            if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(j), 'SHUFFDECODE')
                disp('skipped');
                continue
            elseif isempty(CLASSES.birds(i).neurons(ii).branchnum(j).SHUFFDECODE)
                               disp('skipped - empty');
                continue
            else
                disp('doing!');
            end
            dattmp = CLASSES.birds(i).neurons(ii).branchnum(j).SHUFFDECODE.timebin(1);
            
            decodestruct = struct;
            decodestruct.window_relonset = dattmp.window_relonset;
            decodestruct.ConfMatAll_DAT = dattmp.ConfMatAll_DAT;
            decodestruct.ConfMatAll_NEG = dattmp.ConfMatAll_NEG;
            decodestruct.Pdat = dattmp.Pdat;
            
            save(savefname, 'decodestruct');
            
        end
    end
    
end

%%

savename_par = [savedir '/Params.mat'];
TimeWindows = CLASSES.SHUFFDECODEpar.TimeWindows_relOnsetOffset;
save(savename_par, 'TimeWindows')