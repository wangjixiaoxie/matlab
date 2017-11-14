function [tokenExtents, startinds, endinds, matchlabs] = lt_batchsong_regexp(labels, regexpr_str)
%% lt 11/14/17, finds occurance of syl motif

% regexpstr = 'ab(h)'; % follows regexp rules. but have to have one syl as a token (in parantheses).
% will output position in labels of that token, as well as the first syl in
% motif. This is better than matlab's regexp becuase it works for repeats.
% e.g. if a bird sings abbbb and you search for 'b(b)', then it will pick
% out positions 3, 4, 5. Matlab's wil only pick out 3, 5 (since it does not
% overlap). This code uses matlab's but with lookahead assertion.

% labels = string o labels

%% RUn
% ========== find token
indstmp = regexp(regexpr_str, '[\(\)]'); % find left and right parantheses

% =========== confirm that there is exactly one token syl.
assert(length(indstmp)==2,'sdafasd');
assert(indstmp(2)-indstmp(1) ==2, 'sdfsdaf');

% ====== rewrite rexeprt string for lookahead assertion]
tmptmp = [regexpr_str(1:indstmp(1)-1) ...
    regexpr_str(indstmp(1)+1) regexpr_str(indstmp(2)+1:end)];
regexpr_str2 = [tmptmp(1) '(?=' tmptmp(2:end) ')']; % for lookahead assertion

% ====== do regexp.
[startinds, ~, ~]=regexp(labels, regexpr_str2, 'start', 'end', ...
    'match');

% ====== find positions of tokens and ends of motifs.
strlength = length(regexpr_str)-2;

tokenExtents = startinds + indstmp(1)-1; % i.e. where token was
endinds = startinds + strlength -1;

% ========== get match syls
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
    matchlabs = mat2cell(labels(indmat)', ones(size(indmat,1),1))';
else
    matchlabs = mat2cell(labels(indmat), ones(size(indmat,1),1))';
end
