function [Params,RawDatStruct_SeqFilt]=lt_seq_dep_pitch_SeqFilter_singleday(Params,RawDatStruct)
%% LT 1/31/15 - modified from lt_compile_seq_dep_pitch_data_SEQFILTER
% takes single day structure, filters it, and spits out single day
% structure. Does not save or do anything permanent

%% LT 11/13/14 - Use on structure made with function lt_compile_seq_dep_pitch_data, to filter syls(e.g. b) to sequence specific syl (e.g. a[b]);

% INPUTS
% RawDatStruct - structure that is output from
% lt_compile_seq_dep_pitch_data. This contains data for the single syls
% (e.g. b) of pitch, contour, sequence, etc.
% SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should
% combine. in this case the two sequences are ac[c]b and ab[b].
% SylTargList={'c','b'}; % these must already have raw data compiled above
% SeqPostList={'b',''};
%
% OUTPUTS:
% RawDatStruct - overrides input structure by adding fields
% related tot the new sequences

%% EXTRACT PARAMS

SeqPreList=Params.SeqFilter.SeqPreList;
SylTargList=Params.SeqFilter.SylTargList;
SeqPostList=Params.SeqFilter.SeqPostList;
curr_dir=pwd;

%% FILTER syllables based on sequence context (e.g. only b from ab);
% check whether seq filter is supposed to be done (i.e. this code does both
% sequence and repeat filtering)
% if length(SylTargList)

for j=1:length(SylTargList);
    SeqPre=SeqPreList{j};
    SylTarg=SylTargList{j};
    SeqPost=SeqPostList{j};
    
    % Automatic params
    TargStr=[SeqPre SylTarg SeqPost]; % entire string
    TartPos=length(SeqPre)+1; % target syl position in string
    
    % field names for structure - capitalize the targ syl only
    SeqFieldName=[SeqPre upper(SylTarg) SeqPost];
    
    % First make sure this targ syl is present this day.
    
    
    % DATA
    NumRends=size(RawDatStruct.data.(SylTarg),1);
    IndsToKeep=[];
    for ii=1:NumRends; % all rends of the targ syl
        sylpos=RawDatStruct.data.(SylTarg){ii,8}; % ind in song of candidate
        IndsToCheck=(sylpos-length(SeqPre)):(sylpos+length(SeqPost)); % what within song inds to look at
        
        if IndsToCheck(1)>0 && IndsToCheck(end)<=length(RawDatStruct.data.(SylTarg){ii,7}); % sequence is expected to be within the 1st and last indices of the actual data.
            SeqToCheck=RawDatStruct.data.(SylTarg){ii,7}(IndsToCheck); % is this the desired string?
            
            if strcmp(TargStr,SeqToCheck)==1; % if is desired string, then keep this data pt
                IndsToKeep=[IndsToKeep ii];
            end
        end
    end
    
    if ~isempty(IndsToKeep);
        RawDatStruct.data.(SeqFieldName)=...
            RawDatStruct.data.(SylTarg)(IndsToKeep,:);
    end
    
    % DECIDED NOT to filter data that has not had outliers removed. - to
    % save file space
    
    % OUTLIERS REMOVED
    %     NumRends=size(RawDatStruct.data_WithOutlier.(SylTarg),1);
    %     IndsToKeep=[];
    %     for ii=1:NumRends;
    %         sylpos=RawDatStruct.data_WithOutlier.(SylTarg){ii,8}; % ind in song of candidate
    %         IndsToCheck=(sylpos-length(SeqPre)):(sylpos+length(SeqPost)); % what within song inds to look at
    %
    %         if IndsToCheck(1)>0 && IndsToCheck(end)<=length(RawDatStruct.data.(SylTarg){ii,7});
    %             SeqToCheck=RawDatStruct.data_WithOutlier.(SylTarg){ii,7}(IndsToCheck); % is this the desired string?
    %
    %             if strcmp(TargStr,SeqToCheck)==1; % if is desired string, then keep this datapt
    %                 IndsToKeep=[IndsToKeep ii];
    %             end
    %         end
    %     end
    %
    %     RawDatStruct.data_WithOutlier.(SeqFieldName)=...
    %         RawDatStruct.data_WithOutlier.(SylTarg)(IndsToKeep,:);
end



%% FILTER BASED ON REPEATS

% WARNING, does not have
% ability to filter based on post syl.
% METHOD 1 - similar to sequence filter. takes each rendition, e.g. acbbB,
% as one datapoint. cannot filter based on post syl.  does not put all
% renditions of a single repeat instance together.
if isfield(Params.SeqFilter,'Repeats');
    disp('WARNING - compiling repeats only works if all repeat rends have 20 or fewer');
    
    % how many repeat classes?
    NumRepClasses=length(Params.SeqFilter.Repeats);
    
    for i=1:NumRepClasses;
        Repeats = Params.SeqFilter.Repeats{i};
        
        % find upper case - target syl, and/or pre syls
        ind=find(lower(Repeats)~=Repeats);
        
        % error if mroe than one ind:
        if length(ind)>1
            disp(['error, repeat ' Repeats ' has more than one capitalized letter - Rerun correctly']);
        else
            
            % === parse repeat sequence into target, presyls, and postsyls
            targsyl=lower(Repeats(ind));
            if ind>1;
                presyl=Repeats(1:ind-1);
            else
                presyl = '';
            end
            
            if ind<length(Repeats);
                postsyl=Repeats(ind+1:end);
            else
                postsyl='';
            end
            
            
            
            % === FILTER OUT ALL REPEATS LENGTHS (iteratively)
            % method - for each rendition of the target, loops through long
            % to short repeats. once matches, gets data then moves on to
            % next rendition, that way a syllable is never used >1 time.
            % e.g. if acBj is repeat, then for a rendition it could be
            % abccccBj, where B is the rendition detected.
            NumRends=size(RawDatStruct.data.(targsyl),1);
            IndsToKeep=[];
            RepeatsToTry=linspace(20,1,20); % start with 20. assumes that no song will ever have repeat larger than 20.
            
            RendsNotMatchedToRepeat=[]; % will make a tally.
            
            for ii=1:NumRends; % all rends of the targ syl
                sylpos=RawDatStruct.data.(targsyl){ii,8}; % ind in song of candidate
                WasThisRendMatched=0;
                
                for iii=1:length(RepeatsToTry);
                    currentrepeat=RepeatsToTry(iii);
                    
                    % what sequence preceding and including the target
                    % shoudl look like
                    desiredtargmotif=[presyl repmat(targsyl,1,currentrepeat)];
                    
                    % is this a match?
                    IndsToCheck=(sylpos-length(presyl)-currentrepeat+1:sylpos);
                    
                    if IndsToCheck(1)<1; % then is problem (starts before actual song first syl)
                        continue % try next lower number of repeats
                    elseif IndsToCheck(end)>length(RawDatStruct.data.(targsyl){ii,7}); % then is problem, end of desired sequence would go past end of song.
                        continue
                    elseif ~strcmp(RawDatStruct.data.(targsyl){ii,7}(IndsToCheck), desiredtargmotif); % fail, then no match
                        continue
                    else % this is match,
                        WasThisRendMatched=1;    % Note down that the repeat type of this rendition was determined.
                        SeqFieldName=[presyl repmat(targsyl,1,currentrepeat-1) upper(targsyl)]; % e.g. converts jbbb to jbbB
                        
                        % save to structure
                        if isfield(RawDatStruct.data,SeqFieldName);
                            RawDatStruct.data.(SeqFieldName)=[RawDatStruct.data.(SeqFieldName); RawDatStruct.data.(targsyl)(ii,:)];
                        else  % make new field if first encounter of repeat
                            RawDatStruct.data.(SeqFieldName)=RawDatStruct.data.(targsyl)(ii,:);
                        end
                        
                        break % because this rendition has been classified.
                        
                    end
                end
                
                % -- if not matched to repeat, then note that down
                if WasThisRendMatched==0;
                    RendsNotMatchedToRepeat=[RendsNotMatchedToRepeat ii];
                end
                
            end
        end
    end
    
    RawDatStruct.repeats.RendsNotMatchedToRepeat=RendsNotMatchedToRepeat;
    
end


if (0) % STILL IN PROGRESS - stopped at line 272 - decided to write a more general general expressions code that works for all both sequence and repeats.
    % =================  LOOK FOR REPEATS
    % METHOD 2 - creates array, each row is one rendition of the motif.
    % Goes through all renditions of the single syl. For each rend checks
    % whether it is the final rend in a repeat (because the next syl is NOT
    % that syl, or the next syl is the defined post syl, if it is defined).
    % Then sticks that motif in, and saves FF for all those rends. Also saves
    % the rend number so can go back and find raw data if required.
    
    if isfield(Params.SeqFilter,'Repeats');
        disp('WARNING - compiling repeats only works if all repeat rends have 20 or fewer');
        
        % how many repeat classes?
        NumRepClasses=length(Params.SeqFilter.Repeats);
        
        for i=1:NumRepClasses;
            Repeats = Params.SeqFilter.Repeats{i};
            
            % find upper case - target syl, and/or pre syls
            ind=find(lower(Repeats)~=Repeats);
            
            % error if mroe than one ind:
            if length(ind)>1
                disp(['error, repeat ' Repeats ' has more than one capitalized letter - Rerun correctly']);
            else
                
                % === parse repeat sequence into target, presyls, and postsyls
                targsyl=lower(Repeats(ind));
                if ind>1;
                    presyl=Repeats(1:ind-1);
                else
                    presyl = '';
                end
                
                if ind<length(Repeats);
                    postsyl=Repeats(ind+1:end);
                    NumPostSyls=length(postsyl);
                else
                    postsyl=['[^' targsyl ']']; % for regular expression.
                    NumPostSyls=1;
                end
                
                
                
                % === FIND FINAL REND IN EACH REPEAT RENDITION.
                % method - for each rendition of the target, loops through long
                % to short repeats. once matches, that has to be the last syl in rendition.
                % e.g. if acBj is repeat, then for a rendition it could be
                % abccccBj, where B is the rendition detected.
                NumRends=size(RawDatStruct.data.(targsyl),1);
                RepeatsToTry=linspace(20,1,20); % start with 20. assumes that no song will ever have repeat larger than 20.
                
                RendsNotMatchedToRepeat=[]; % will make a tally.
                
                for ii=1:NumRends; % all rends of the targ syl
                    sylpos=RawDatStruct.data.(targsyl){ii,8}; % ind in song of candidate
                    WasThisRendMatched=0;
                    
                    for iii=1:length(RepeatsToTry);
                        currentrepeat=RepeatsToTry(iii);
                        
                        % === what sequence preceding and including the target
                        % shoudl look like
                        desiredtargmotif=[presyl repmat(targsyl,1,currentrepeat)];
                        % add on either post syl, or a not(regexp)
                        desiredtargmotif=[desiredtargmotif postsyl];
                        
                        % is this a match?
                        IndsToCheck=(sylpos-length(presyl)-currentrepeat+1:sylpos+NumPostSyls);
                        
                        if IndsToCheck(1)<1; % then is problem (starts before actual song first syl)
                            continue % try next lower number of repeats
                        elseif IndsToCheck(end)>length(RawDatStruct.data.(targsyl){ii,7}); % then is problem, end of desired sequence would go past end of song.
                            continue
                        else
                            string_match=regexp(RawDatStruct.data.(targsyl){ii,7}(IndsToCheck), desiredtargmotif, 'once'); % fail, then no match
                            if isempty(string_match); % then failed to match, try next repeat number
                                continue
                                
                            else % THIS IS MATCH!!!! (i.e. this rend is the last repeated syllable in a repeat rendition.
                                WasThisRendMatched=1;    % Note down that the repeat type of this rendition was determined.
                                SeqFieldName=[presyl repmat(targsyl,1,currentrepeat-1) upper(targsyl)]; % e.g. converts jbbb to jbbB
                                
                                % save to structure
                                if isfield(RawDatStruct.data,SeqFieldName);
                                    RawDatStruct.data.(SeqFieldName)=[RawDatStruct.data.(SeqFieldName); RawDatStruct.data.(targsyl)(ii,:)];
                                else  % make new field if first encounter of repeat
                                    RawDatStruct.data.(SeqFieldName)=RawDatStruct.data.(targsyl)(ii,:);
                                end
                                
                                break % because this rendition has been classified.
                            end
                        end
                    end
                    
                    % -- if not matched to repeat, then note that down
                    if WasThisRendMatched==0;
                        RendsNotMatchedToRepeat=[RendsNotMatchedToRepeat ii];
                    end
                    
                end
            end
        end
        
        RawDatStruct.repeats.RendsNotMatchedToRepeat=RendsNotMatchedToRepeat;
        
        
    end
end


%% SUMMARY STATISTICS - e.g. mean, SD, etc.


syl_fields=fields(RawDatStruct.data); % what seqs are there?
for i=1:length(syl_fields);
    
    X=cell2mat(RawDatStruct.data.(syl_fields{i})(:,1));
    RawDatStruct.summary_stats.(syl_fields{i}).meanFF=mean(X);
    RawDatStruct.summary_stats.(syl_fields{i}).medianFF=median(X);
    RawDatStruct.summary_stats.(syl_fields{i}).sdFF=std(X);
    RawDatStruct.summary_stats.(syl_fields{i}).n=length(X);
    RawDatStruct.summary_stats.(syl_fields{i}).semFF=std(X)/sqrt(length(X)-1);
    
    % get mean pitch contour
    X=cell2mat(RawDatStruct.data.(syl_fields{i})(:,2)); % rows: rend, cols: timebin
    RawDatStruct.summary_stats.(syl_fields{i}).pitchcountour_mean=mean(X,1);
    
    %     WILL IGNORE GETTING MEAN SPECTROGRAM - large file size.
    % get mean spectrogram
    %     X=[];
    %     for ii=1:size(RawDatStruct.data.(syl_fields{i})(:,4),1); % num rends
    %         X(:,:,ii)=RawDatStruct.data.(syl_fields{i}){ii,4}; % extract individual specs
    %     end
    %     RawDatStruct.summary_stats.(syl_fields{i}).spec_mean=mean(X,3);
    
    
    % IGNORE DATA THAT HAS OUTLIERS.
    %     % NO OUTLIERS
    %     X=cell2mat(RawDatStruct.data_WithOutlier.(syl_fields{i})(:,1));
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).meanFF=mean(X);
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).medianFF=median(X);
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).sdFF=std(X);
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).n=length(X);
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).semFF=std(X)/sqrt(length(X)-1);
    %
    %     % get mean pitch contour
    %     X=cell2mat(RawDatStruct.data_WithOutlier.(syl_fields{i})(:,2)); % rows: rend, cols: timebin
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).pitchcountour_mean=mean(X,1);
    %
    %     % get mean spectrogram
    %     X=[];
    %     for ii=1:size(RawDatStruct.data_WithOutlier.(syl_fields{i})(:,4),1); % num rends
    %         X(:,:,ii)=RawDatStruct.data_WithOutlier.(syl_fields{i}){ii,4}; % extract individual specs
    %     end
    %     RawDatStruct.summary_stats_WithOutlier.(syl_fields{i}).spec_mean=mean(X,3);
end


%% OUTPUT
RawDatStruct_SeqFilt=RawDatStruct;

