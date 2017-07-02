function SummaryStruct = lt_neural_v2_PostInfo(SummaryStruct)



%% lt 5/3/17 - adds info onto SummaryStruct


NumBirds = length(SummaryStruct.birds);
failures = 0;
for j=1:NumBirds
    
    NumNeurons = length(SummaryStruct.birds(j).neurons);
    
    for jj=1:NumNeurons
        
        birdname = SummaryStruct.birds(j).birdname;
        exptname = SummaryStruct.birds(j).neurons(jj).exptID;
        
        
        % === WHAT MOTIFS FOR THIS BIRD/EXPT
        MotifsActual = {};
        if strcmp(birdname, 'wh6pk36')
            MotifsActual = {'nlcchb', 'jklcchb', 'ga', 'mksd' ,'vb'};
        elseif strcmp(birdname, 'bk7')
            MotifsActual = {'nnhh', 'gh', 'vbbb', 'gb', 'jkk', 'kl', 'gv', 'yoo', 'rs'}; % COULD IMPROVE
        elseif strcmp(birdname, 'bu77wh13')
%             MotifsActual = {'kspj', 'ab', 'bh', 'jb', 'nkrs'}; % COULD IMPROVE
            MotifsActual = {'ab', 'jb', 'ijbhnkrs', 'oks', 'pjabbhgks'};
        elseif strcmp(birdname, 'br92br54')
            MotifsActual = {'nkh', 'ddd', 'dh', 'agc', 'cc', 'cd', 'hd', 'sd', 'md'}; % COULD IMPROVEe
        elseif strcmp(birdname, 'or74bk35');
            MotifsActual = {'abg', 'anb', 'anbg', 'hbg',}; % IN PROGRESS
        else
            disp('NO MOTIFS INPUTED!!');
            failures = failures+1;
        end
        
        
        % ==== BASED ON THOSE MOTIFS FIGURE OUT ALL REG EXP STRINGS
        motif_regexpr_str = {}; % will put motifs for extraction here
        for i=1:length(MotifsActual)
            motif_actual = MotifsActual{i};
            
            % for each vocalization in the motif, extract one segment
            numvocals = length(motif_actual);
            
            for ii=1:numvocals
                if ii==1
                    segmentmotif = ['(' motif_actual(ii) ')' motif_actual(ii+1:end)];
                    
                elseif ii == numvocals
                    segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')'];
                    
                else
                    
                    segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')' motif_actual(ii+1:end)];
                    
                end
                
                %        disp(segmentmotif);
                
                motif_regexpr_str = [motif_regexpr_str segmentmotif];
                
            end
        end
        
        % display output
        if (0)
            disp(['--------------------------------' birdname '-' exptname])
            disp(MotifsActual);
            disp(motif_regexpr_str);
        end
        
        % ========== FIGURE OUT LIST OF SINGLE SYLS
        SingleSylAll = {};
        for i=1:length(motif_regexpr_str)
           motifsyl = motif_regexpr_str{i};
           ind = strfind(motifsyl, '(');
            singlesyl = motifsyl(ind+1);
            SingleSylAll = [SingleSylAll singlesyl];
        end
        SingleSylAll = unique(SingleSylAll);
        
        % ==== SAVE OUTPUT INTO STRUCT
        SummaryStruct.birds(j).neurons(jj).POSTINFO.MotifsActual = MotifsActual;
        SummaryStruct.birds(j).neurons(jj).POSTINFO.MotifsActual_regexpStr = motif_regexpr_str;
         SummaryStruct.birds(j).neurons(jj).POSTINFO.SingleSyls = SingleSylAll;
       
    end
end

if failures>0
    lt_figure;
    lt_plot_text(0, 0.5, 'FAILURE!, no motif specificed');
end
