function lt_neural_v2_LEARNING_NeurPairCorr(MOTIFSTATS_Compiled, SwitchStruct)




%% go thru all switches. if have greater than N simult neurons (with concurrent trials) then calculate pairwise corr for each syl

Numbirds = length(SwitchStruct.bird);

for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
%         SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
%         SameSyls = MotifStats.params.SameTypeSyls;
%         DiffSyls = MotifStats.params.DiffTypeSyls;
%         
%         targsyls = MotifStats.params.TargSyls;
        nummotifs = length(motiflist);
        
%         WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
%             MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            

            goodneurons = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons;
            
            if isempty(goodneurons)
                disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                continue
            end
            
%             swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
%             swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
%             swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
%             
%             plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
            
            % - things about target
%             numtargs = length(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs)/2;
%             targssamesyl = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl;
            
%             if length(unique([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2:2:end}])) ==1
%                 TargsSameDir = 1;
%             else
%                 TargsSameDir = 0;
%             end
%             
%             if TargsSameDir ==0
%                 targlearndir = nan;
%             else
%                 targlearndir = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2};
%             end
            
            
            % --- skip if multiple targs and in different directions
%             if skipMultiDir ==1
%                 if numtargs>1 & TargsSameDir==0
%                     disp(['SKIPPING ' birdname '-' exptname '-sw' num2str(iii) ' [multidir]']);
%                     continue
%                 end
%             end
            

            % ============== see if there is any overlap between these
            % neurons
            motifnum = 1;
            
            NeuronNumall = [];
            RendIDall = [];
            for j = goodneurons
               
%                 tmp = [num2str(1e4*[MotifStats.neurons(j).motif(motifnum).SegmentsExtract.song_datenum]') ...
%                 num2str(ceil(1e4*[MotifStats.neurons(j).motif(motifnum).SegmentsExtract.]'))];
% 
%                 tmp = [1e4*[MotifStats.neurons(j).motif(motifnum).SegmentsExtract.song_datenum]' ...
%                 ceil(1e4*[MotifStats.neurons(j).motif(motifnum).SegmentsExtract.Dur_syl]') ...
%                 ceil([MotifStats.neurons(j).motif(motifnum).SegmentsExtract.Dur_gappost]')];
% 
                RendIDall = [RendIDall; tmp];
                NeuronNumall = [NeuronNumall; j*ones(size(tmp,1), 1)];
                
                if any(ceil(1e6*[MotifStats.neurons(j).motif(motifnum).SegmentsExtract.Dur_syl]')==85467)
                    keyboard
                end
                
            end
            
            RendIDall = mat2cell(RendIDall, ones(size(RendIDall,1),1), size(RendIDall,2)); % convert to cell, each one a rend
            RendIDunique = unique(RendIDall);
            
            % --- for each Rend, determine if all neurons have it.
            tabulate(RendIDall)
            NeuronNumall(strcmp(RendIDall, '7368278980.208389467'))
            
            
%             SwitchStruct.bird(1).exptnum(1).switchlist(1).neuron(1).DATA.motif(1).tvals
%             
%             {MOTIFSTATS_Compiled.birds(i).exptnum(1).MOTIFSTATS.neurons(1).motif(1).SegmentsExtract.song_filename};
%             {MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.neurons(2).motif(1).SegmentsExtract.song_filename}
            
        end
    end
end
           