function SummaryStruct = lt_neural_v2_PostInfo(SummaryStruct, IgnoreWhetherLearn, ...
    MotifsToCollect)
%% lt 1/17/18 - added ability to manually enter which motifs want to extract
% note: will only use this if bird is entered. if not, then will still do
% automatic motif extraction

if ~exist('MotifsToCollect', 'var')
    MotifsToCollect = [];
end

% e.g. MotifsToCollect = {'pu69wh78', {'(j)jbhhg', 'aabh'}};
% bird cell pair. cell contains strings. if a given string has parantheses then is a specific regexp. if no
% parantheses then will automaticalyl extract all strings for that "motif"

%% lt 11/29/17 - added IgnoreWhetherLearn
% IgnoreWhetherLearn = 1; % then uses same motifs for all experiments for
% a given bird, regardless of whether is learning. default is 0, i.e. takes
%  learning specific motifs

if ~exist('IgnoreWhetherLearn', 'var')
    IgnoreWhetherLearn=0;
end


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
        
        if IgnoreWhetherLearn==1
            islearning=0;
        end
        
        if islearning ==1
            %% ##### USE LEARNING SPECIFIC MOTIFS (accounts for disruption by WN)
            if strcmp(birdname, 'wh6pk36')
                if strcmp(exptname, 'LMANlearn2')
                    %                     MotifsActual = {'mnl', 'mksd', 'jkl', 'vbga', 'cchb'}; % OK - NOTE; c1 has diff preceding...
                    MotifsActual = {'mksd', 'jkl', 'vbga', ...
                        'mnl', 'chb', 'kl(c)', 'nl(c)'}; % NOTE: CHEKCEKD, GOOD ...
                end
                
                
            elseif strcmp(birdname, 'bk7')
                if strcmp(exptname, 'LearnLMAN1') | strcmp(exptname, 'LearnLMAN2')
                    MotifsActual = {'nhh', 'j(k)', 'jk(k)', 'kl', 'y(o)', 'yo(o)', 'o(p)', 'op(p)', ...
                        'b(y)', ...
                        'h(j)', 'p(j)', ...
                        'g(h)', ...
                        'hh(r)', 'p(r)', 'l(r)', 'c(r)', ...
                        'r(s)', ...
                        'hr(g)', 'b(g)', 'l(g)', 's(g)', 'p(g)', 'h(g)', 'v(g)', 'c(g)',...
                        'g(v)', 'u(v)', ...
                        'gv(b)', 'gvb(b)', 'g(b)', 'gb(b)', 'c(b)', ...
                        'o(t)', 's(t)', ...
                        'f(c)', 'g(c)'}; % NOTE: good
                end
                
                
            elseif strcmp(birdname, 'bu77wh13')
                if strcmp(exptname, 'LMANlearn1')
                    %                     MotifsActual = {'pjab', 'bhgks', 'ijbhnk', 'oks',
                    %                     'rs'}; % OLD VERSION -
                    MotifsActual = {'pjab', 'ijbhnk', 'rs', ...
                        'ok', 'gk', 's(k)', ...
                        'k(s)'};
                end
                
                
            elseif strcmp(birdname, 'br92br54')
                if strcmp(exptname, 'LMANlearn2') | strcmp(exptname, 'LMANlearn4')
                    MotifsActual = {'gcc', ...
                        'm(s)', 's(s)', 'd(s)', 'c(s)', ...
                        'k(h)', 'd(h)',...
                        'n(k)', 'm(k)', ...
                        'o(n)', ...
                        'c(d)','d(d)d', 'm(d)', ...
                        'd(m)', 's(m)', 'm(m)', 'c(m)', 'o(m)', ...
                        'm(a)', 'n(a)', ...
                        's(o)', 'y(o)', 'yo(o)', ...
                        's(y)'}; % NEW, looks good.
                elseif strcmp(exptname, 'LMANlearn3')
                    %                     MotifsActual = {'ddd', 'agc', 'nkh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                    MotifsActual = {'gc', ...
                        'm(s)', 's(s)', 'd(s)', 'c(s)', ...
                        'k(h)', 'd(h)',...
                        'n(k)', 'm(k)', 'h(k)', ...
                        'o(n)', 's(n)', ...
                        'c(d)','d(d)d', 'h(d)', 'm(d)', ...
                        'd(m)', 's(m)', 'm(m)', 'c(m)', 'h(m)', 'o(m)', ...
                        'm(a)', 'n(a)', ...
                        's(o)', 'y(o)', 'yo(o)', ...
                        's(y)', 'h(y)'};
                elseif strcmp(exptname, 'LMANlearn5')
                    %                     MotifsActual = {'(d)dd', 'd(d)d', 'agcc', 'nkh', 'dh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                    MotifsActual = {'gcc', ...
                        'm(s)', 's(s)', 'd(s)', 'c(s)', ...
                        'k(h)', 'd(h)',...
                        'n(k)', 'm(k)', ...
                        'o(n)', 's(n)', ...
                        'c(d)','d(d)d', 'm(d)', ...
                        'd(m)', 's(m)', 'm(m)', 'c(m)', 'o(m)', ...
                        'm(a)', 'n(a)', ...
                        's(o)', 'y(o)', 'yo(o)', ...
                        's(y)'}; % NEW, looks good.
                elseif strcmp(exptname, 'LMANlearn6')
%                     MotifsActual = {'h(d)', 'c(d)', 'dd(d)', 'agcc', ...
%                         'nkh'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                    MotifsActual = {'gcc', ...
                        'm(s)', 's(s)', 'd(s)', 'c(s)', ...
                        'k(h)', 'd(h)',...
                        'n(k)', 'm(k)', 'h(k)', ...
                        'o(n)', 's(n)', ...
                        'c(d)','dd(d)', 'h(d)', 'm(d)', ...
                        'd(m)', 's(m)', 'm(m)', 'c(m)', 'h(m)', 'o(m)', ...
                        'm(a)', 'n(a)', ...
                        's(o)', 'y(o)', 'yo(o)', ...
                        's(y)', 'h(y)'}; % NEW, looks good.
                elseif strcmp(exptname, 'LMANlearn7')
                    MotifsActual = {'gc', ...
                        'm(s)', 's(s)', 'd(s)', ...
                        'k(h)', 'd(h)',...
                        'n(k)', 'm(k)', 'h(k)', ...
                        'o(n)', 's(n)', ...
                        'd(d)d', 'h(d)', 'm(d)', ...
                        'd(m)', 's(m)', 'm(m)', 'h(m)', 'o(m)', ...
                        'm(a)', 'n(a)', ...
                        's(o)', 'y(o)', 'yo(o)', ...
                        's(y)', 'h(y)'}; % NEW, looks good.
                end
                
            elseif strcmp(birdname, 'or74bk35');
                if strcmp(exptname, 'LMANneural2')
                    MotifsActual = {...
                        'ja(b)', 'ba(b)', 'an(b)', 'h(b)', ...
                        'ja(n)', 'ba(n)', ...
                        'm(a)', 'a(a)', 'j(a)', ...
                        'g(j)', 'a(j)', ...
                        'a(h)', 'ah(h)'}; % GOOD
                    
                elseif strcmp(exptname, 'LMANlearn3')
                    MotifsActual = {'(a)aa', 'a(a)a', '(a)ba', 'a(b)a', 'abg', ...
                        '(a)nba', 'a(n)ba', 'an(b)a', 'anbg', ...
                        'ag', 'am', 'hbg'}; %
                end
                
            elseif strcmp(birdname, 'pu69wh78')
                if strcmp(exptname, 'RALMANlearn1')
                    MotifsActual = {...
                        'aab', ...
                        'jjbhh', 'jjbhh(g)'};
                elseif strcmp(exptname, 'RALMANlearn2')
                    MotifsActual = {...
                        'aab', ...
                        'jjbhh', 'jjbhh(g)'};
                elseif strcmp(exptname, 'RALMANOvernightLearn1')
                    MotifsActual = {...
                        'aabh', ...
                        'jjbhh', 'jjbhh(g)'};
                elseif strcmp(exptname, 'RAlearn1')
                    MotifsActual = {...
                        'aab', ...
                        'jjb', ...
                        'h(g)'};
                end
                
            elseif strcmp(birdname, 'wh44wh39')
                if strcmp(exptname, 'RALMANlearn1')
                    MotifsActual = {'nhh', 'dkccb', ...
                        'dkccb(b)', ...
                        '(j)n', '(m)d'};
                end
                
            else
                disp('NO MOTIFS INPUTED!!');
                failures = failures+1;
            end
            
        else
            %% ###### THEN USE GENERAL MOTIFS (with no WN disruption)
            if strcmp(birdname, 'wh6pk36')
                MotifsActual = {'mksd', 'jkl', 'vbga', ...
                    'mnl', 'chbga', 'kl(c)', 'nl(c)'}; % NOTE: CHEKCEKD, GOOD ...
                
            elseif strcmp(birdname, 'bk7')
                %                 MotifsActual = {'nnhh', 'h(h)r', 'rs', 'gn', 'hpp', '(h)jk', 'jkk', 'kkl', ...
                %                         'gvb', 'gh', 'gbb', 'yoo', 'opt'}; % OLD VERSION, GOT TOO MESSY
                MotifsActual = {'nhh', 'j(k)', 'jk(k)', 'kl', 'y(o)', 'yo(o)', 'we', ...
                    'o(p)', 'op(p)', 'h(p)', ...
                    'w(e)', 'j(e)', ...
                    'b(y)', ...
                    'h(j)', 'p(j)', 'r(j)', ...
                    'g(h)', 'gh(h)' ...
                    'hh(r)', 'p(r)', 'l(r)', 'c(r)', ...u
                    'r(s)', 'h(s)', ...
                    'r(g)', 'b(g)', 'l(g)', 's(g)', 'p(g)', 'h(g)', 'v(g)', 'c(g)',...
                    'g(v)', 'u(v)', 't(v)', ...
                    'gv(b)', 'gvb(b)', 'g(b)', 'gb(b)', 'm(b)', 'c(b)', ...
                    'o(t)', 's(t)', ...
                    'v(m)', ...
                    'e(u)', ...
                    'f(c)', 'g(c)'};
                
                
            elseif strcmp(birdname, 'bu77wh13')
                %                 MotifsActual = {'pjabbh', 'ijbhnk', 'oks', 'gks', 'rs'};
                MotifsActual = {'pjabbh', 'ijbhnk', 'rs', ...
                    'ok', 'gk', 's(k)', ...
                    'k(s)'};
                
            elseif strcmp(birdname, 'br92br54')
                %                 MotifsActual = {'nkh', 'ddd', 'dh', 'agc', 'cc', 'cd', 'hd', 'sd', 'md'}; % COULD IMPROVEe
                %                 MotifsActual = {'ddd', 'agcc', 'nkh', 'ddd(h)'}; % IMPROVE - NEED TO ADD THINGS (E.G. []d)
                MotifsActual = {'gcc', ...
                    'm(s)', 's(s)', 'd(s)', 'c(s)', ...
                    'k(h)', 'd(h)',...
                    'n(k)', 'm(k)', 'h(k)', ...
                    'o(n)', 's(n)', ...
                    'c(d)','d(d)d', 'h(d)', 'm(d)', ...
                    'd(m)', 's(m)', 'm(m)', 'c(m)', 'h(m)', 'o(m)', ...
                    'm(a)', 'n(a)', ...
                    's(o)', 'y(o)', 'yo(o)', ...
                    's(y)', 'h(y)'}; % NEW, looks good.
                
            elseif strcmp(birdname, 'or74bk35');
                %                 MotifsActual = {'(a)aa', 'a(a)a', '(a)ba', 'a(b)a', 'abgj', ...
                %                     '(a)nba', 'a(n)ba', 'an(b)a', 'anbg', ...
                %                     'agj', 'am', 'hbg'}; %
                MotifsActual = {...
                    'nb(g)', 'ab(g)', ...
                    'ja(b)', 'ba(b)', 'an(b)', 'h(b)', ...
                    'ja(n)', 'ba(n)', ...
                    'm(a)', 'a(a)', 'j(a)', 'b(a)', ...
                    'g(j)', 'a(j)', ...
                    'a(h)', 'ah(h)'}; % GOOD
            elseif strcmp(birdname, 'pu69wh78')
                MotifsActual = {...
                    'aabh', 'aabh(h)', 'aabhh(g)', ...
                    'jjbhh', 'jjbhh(g)'};
            elseif strcmp(birdname, 'wh44wh39')
                MotifsActual = {'nhh', 'dkccbb', ...
                    '(j)n', '(m)d'};
            else
                disp('NO MOTIFS INPUTED!!');
                failures = failures+1;
            end
        end
        
        
        
        
        %% === WHAT MOTIFS FOR THIS BIRD/EXPT
        %         if (0)
        %             MotifsActual = {};
        %             if strcmp(birdname, 'wh6pk36')
        %                 MotifsActual = {'nlcchb', 'jklcchb', 'ga', 'mksd' ,'vb'};
        %             elseif strcmp(birdname, 'bk7')
        %                 MotifsActual = {'nnhh', 'gh', 'vbbb', 'gb', 'jkk', 'kl', 'gv', 'yoo', 'rs'}; % COULD IMPROVE
        %             elseif strcmp(birdname, 'bu77wh13')
        %                 %             MotifsActual = {'kspj', 'ab', 'bh', 'jb', 'nkrs'}; % COULD IMPROVE
        %                 MotifsActual = {'ab', 'jb', 'ijbhnkrs', 'oks', 'pjabbhgks'};
        %             elseif strcmp(birdname, 'br92br54')
        %                 MotifsActual = {'nkh', 'ddd', 'dh', 'agc', 'cc', 'cd', 'hd', 'sd', 'md'}; % COULD IMPROVEe
        %             elseif strcmp(birdname, 'or74bk35');
        %                 MotifsActual = {'abg', 'anb', 'anbg', 'hbg',}; % IN PROGRESS
        %             else
        %                 disp('NO MOTIFS INPUTED!!');
        %                 failures = failures+1;
        %             end
        %         end
        

        
        %% ALTERNATIVE - MANUALLY ENTER MOTIFS
        if ~isempty(MotifsToCollect)
            % then overwrite stuff from above, use hand entered stuff
            indtmp = find(strcmp(MotifsToCollect, birdname));
            if ~isempty(indtmp)
               MotifsActual = MotifsToCollect{indtmp+1};
                
            end
        end
        
%%
%         assert(~isempty(MotifsActual), 'no motifs given...');
        
        if isempty(MotifsActual)
            
            return
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
        

        %%
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
