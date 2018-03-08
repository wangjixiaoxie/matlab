function [startinds, tokenExtents, endinds, matchlabs] = lt_neural_QUICK_regexp(AllLabels, regexpr_str)
%% note:
% if want to use square brackets, currently only works if square brackets
% are the first term, and no token within those brackets 
% e.g. '[a-z]b(b)' works, but (b)[a-z] does not ...

%% lt 3/6/18 - given string, pulls out motifs
% there has to be one token syllable (e.g. 'ab(c)', where c is token)
% 

    indstmp = regexp(regexpr_str, '[\(\)]'); % find left and right parantheses
    
    % confirm that there is exactly one token syl.
    assert(length(indstmp)==2,'sdafasd');
    assert(indstmp(2)-indstmp(1) ==2, 'sdfsdaf');
    
    % -
    tmptmp = [regexpr_str(1:indstmp(1)-1) ...
        regexpr_str(indstmp(1)+1) regexpr_str(indstmp(2)+1:end)];
    
    %     regexpr_str2 = [regexpr_str(1:indstmp(1)-1) '(?=' ...
    %         regexpr_str(indstmp(1)+1) regexpr_str(indstmp(2)+1:end) ')']; % for lookahead assertion
    
% ========== are there are square brackets?
brackinds = [regexp(regexpr_str, '[') regexp(regexpr_str, ']')];
if length(brackinds)==2
    assert(brackinds(1) ==1, 'currently following only works if square brackets is at start of string');
    assert(indstmp(1)>brackinds(2), 'currenrly onlly works if token is to the right of the closing sqare bracket');
    
    % -- then ..
    
    regexpr_str2 = [regexpr_str(brackinds(1):brackinds(2)) '(?=' tmptmp(brackinds(2)+1:end) ')'];
%     
%     regexpr_str2 = [tmptmp(1) '(?=' tmptmp(2:end) ')']; % for lookahead assertion
%     
    [startinds, ~, ~]=regexp(AllLabels, regexpr_str2, 'start', 'end', ...
        'match');
    
    strlength = length(regexpr_str)-2 - (brackinds(2) - brackinds(1));
    
    tokenExtents = startinds + indstmp(1)-1 - (brackinds(2) - brackinds(1)); % i.e. where token was
    endinds = startinds + strlength -1;
    
    % - get match syls
    if size(startinds,2)==1
        startinds = startinds';
    end
    
    indmat = [];
    for j=1:strlength
        indmat = [indmat startinds'+j-1];
    end
    
    %     if size(AllLabels,1)==1
    %     matchlabs = mat2cell(AllLabels(indmat)', ones(size(indmat,1),1))';
    %     else
    if size(indmat,2)==1
        % then is one col (i.e. only one syl in motif) so need to make sure
        % output is col
        matchlabs = mat2cell(AllLabels(indmat)', ones(size(indmat,1),1))';
    else
        matchlabs = mat2cell(AllLabels(indmat), ones(size(indmat,1),1))';
    end
    
elseif length(brackinds)>2
    % have not coded how to deal with this ...
    
    disp('NOT SURE ...');
    return
elseif isempty(brackinds)    
    % ================== THEN THIS IS ORIGNAL VERSION 
    regexpr_str2 = [tmptmp(1) '(?=' tmptmp(2:end) ')']; % for lookahead assertion
    
    [startinds, ~, ~]=regexp(AllLabels, regexpr_str2, 'start', 'end', ...
        'match');
    
    strlength = length(regexpr_str)-2;
    
    tokenExtents = startinds + indstmp(1)-1; % i.e. where token was
    endinds = startinds + strlength -1;
    
    % - get match syls
    if size(startinds,2)==1
        startinds = startinds';
    end
    
    indmat = [];
    for j=1:strlength
        indmat = [indmat startinds'+j-1];
    end
    
    %     if size(AllLabels,1)==1
    %     matchlabs = mat2cell(AllLabels(indmat)', ones(size(indmat,1),1))';
    %     else
    if size(indmat,2)==1
        % then is one col (i.e. only one syl in motif) so need to make sure
        % output is col
        matchlabs = mat2cell(AllLabels(indmat)', ones(size(indmat,1),1))';
    else
        matchlabs = mat2cell(AllLabels(indmat), ones(size(indmat,1),1))';
    end
    
    % ---- compare methods
    if (0)
        if length(startinds) == length(startinds1)
            assert(isempty(setxor(endinds, endinds1)), 'asfasd');
            assert(isempty(setxor(startinds1, startinds)), 'asfasd');
            assert(isempty(setxor(matchlabs1, matchlabs)), 'asdfsd');
            assert(all(tokenExtents1 == tokenExtents), 'asdfsd');
            
        else
            assert(isempty(setdiff(startinds1, startinds)), 'asfasdf'); % i.e. old version is proper subset of new version
        end
    end
else
    disp('NOT SURE ...');
    return
    
end