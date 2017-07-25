function SummaryStruct = lt_neural_v2_PostInfo(SummaryStruct)

%% lt modified on 6/26/17 to be accurate for elarning experiments
% PRoblem - WN over motifs means can't have entire motif name for elaring
% experiments. solution - enter motifs separately for each learning expt.
% NOTE: if want to analyze encoding would look at pre-learning for WN
% training, and use the bird's general motifs (which would always be used
% for experiments without WN)
% NOTE: still not great becuase does not account for variable transitions
% betwen motifs.

% NOTE; if motif has no (), then will go through all syls. if has (), then
% will only take that syl


%% lt 5/3/17 - adds info onto SummaryStruct


NumBirds = length(SummaryStruct.birds);
failures = 0;
for j=1:NumBirds
    
    NumNeurons = length(SummaryStruct.birds(j).neurons);
    
    for jj=1:NumNeurons
        
        birdname = SummaryStruct.birds(j).birdname;
        exptname = SummaryStruct.birds(j).neurons(jj).exptID;
        
        %% === what motifs, modified 6/26, see above.
        
        islearning = SummaryStruct.birds(j).neurons(jj).INFO_islearning;
        MotifsActual = {};
        
        if islearning ==1
            
            if strcmp(birdname, 'wh6pk36')
                if strcmp(exptname, 'LMANlearn2')
                    MotifsActual = {'mnl', 'mksd', 'jkl', 'vbga', 'cchb'}; % OK - NOTE; c1 has diff preceding...
                end
                
                
            elseif strcmp(birdname, 'bk7')
                if strcmp(exptname, 'LearnLMAN1') | strcmp(exptname, 'LearnLMAN2')
                    MotifsActual = {'nnhh', 'hrs', 'gn', 'hpp', 'jkk', 'kkl', ...
                        'gvb', 'gh', 'gbb', 'yoo', 'opt'}; % DIFFICULT... (fits more the later days)
                end
                
                
            elseif strcmp(birdname, 'bu77wh13')
                if strcmp(exptname, 'LMANlearn1')
                    MotifsActual = {'pjab', 'bhgks', 'ijbhnk', 'oks', 'rs'};
                end
                
                
            elseif strcmp(birdname, 'br92br54')
                if strcmp(exptname, 'LMANlearn2') | strcmp(exptname, 'LMANlearn4')
                    MotifsActual = {'ddd', 'agcc', 'nkh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                elseif strcmp(exptname, 'LMANlearn3')
                    MotifsActual = {'ddd', 'agc', 'nkh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                elseif strcmp(exptname, 'LMANlearn5')
                    MotifsActual = {'(d)dd', 'd(d)d', 'agcc', 'nkh', 'dh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                elseif strcmp(exptname, 'LMANlearn6')
                    MotifsActual = {'h(d)', 'c(d)', 'dd(d)', 'agcc', ...
                        'nkh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                end
                
                
            elseif strcmp(birdname, 'or74bk35');
                if strcmp(exptname, 'LMANneural2')
                    MotifsActual = {'(a)aa', 'a(a)a', '(a)ba', 'a(b)a', 'abgj', ...
                        'anb', 'agj', 'am', 'hbg'}; %
                    
                elseif strcmp(exptname, 'LMANlearn3')
                    MotifsActual = {'(a)aa', 'a(a)a', '(a)ba', 'a(b)a', 'abg', ...
                        '(a)nba', 'a(n)ba', 'an(b)a', 'anbg', ...
                        'ag', 'am', 'hbg'}; %
                end
                
            else
                disp('NO MOTIFS INPUTED!!');
                failures = failures+1;
            end
            
        else
            % THEN USE GENERAL MOTIFS
            if strcmp(birdname, 'wh6pk36')
                MotifsActual = {'mnl', 'mksd', 'jkl', 'vbga', 'cchbga'}; % OK - NOTE; c1 has diff preceding...
                
            elseif strcmp(birdname, 'bk7')
                MotifsActual = {'nnhh', 'hrs', 'gn', 'hpp', 'jkk', 'kkl', ...
                    'gvb', 'ghh', 'gbb', 'yoo', 'opt'}; % DIFFICULT... (fits more the later days)
                
                
            elseif strcmp(birdname, 'bu77wh13')
                MotifsActual = {'pjabbh', 'ijbhnk', 'oks', 'gks', 'rs'};
                %                 MotifsActual = {'ab', 'jb', 'ijbhnkrs', 'oks',
                %                 'pjabbhgks'}; this
                
            elseif strcmp(birdname, 'br92br54')
                %                 MotifsActual = {'nkh', 'ddd', 'dh', 'agc', 'cc', 'cd', 'hd', 'sd', 'md'}; % COULD IMPROVEe
                MotifsActual = {'ddd', 'agcc', 'nkh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                
                
            elseif strcmp(birdname, 'or74bk35');
                %                 MotifsActual = {'abg', 'anb', 'anbg', 'hbg',}; % IN PROGRESS
                MotifsActual = {'(a)aa', 'a(a)a', '(a)ba', 'a(b)a', 'abgj', ...
                    '(a)nba', 'a(n)ba', 'an(b)a', 'anbg', ...
                    'agj', 'am', 'hbg'}; %
                
            else
                disp('NO MOTIFS INPUTED!!');
                failures = failures+1;
            end
        end
        
        assert(~isempty(MotifsActual), 'no motifs given...');
        
        
        %% === WHAT MOTIFS FOR THIS BIRD/EXPT
        if (0)
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
        end
        %%
        % ==== BASED ON THOSE MOTIFS FIGURE OUT ALL REG EXP STRINGS
        motif_regexpr_str = {}; % will put motifs for extraction here
        for i=1:length(MotifsActual)
            motif_actual = MotifsActual{i};
            
            if ~isempty(strfind(motif_actual, '('));
                % then just want this syl
                
                motif_regexpr_str = [motif_regexpr_str motif_actual];
                
                
            else
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
        end
        
%         disp(motif_regexpr_str);
%         if any(strcmp(motif_regexpr_str, 'h(()d)'))
%             keyboard
%         end
        
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
        SummaryStruct.birds(j).neurons(jj).POSTINFO.SingleSyls_unique = SingleSylAll;
        
    end
end

if failures>0
    lt_figure;
    lt_plot_text(0, 0.5, 'FAILURE!, no motif specificed');
end
