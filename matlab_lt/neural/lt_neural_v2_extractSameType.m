function [SameTypeSyls, DiffTypeSyls, motif_regexpr_str, SingleSyls] = ...
    lt_neural_v2_extractSameType(motif_regexpr_str, TargSyls)

%% lt 6/29/17 -
% modified so takes in motif_regexpr_str instead of MotifsActual

%%
% TargSyls = {'ab(c)', 'x(y)'};
% MotifsActual = {'abc', 'xyz', 'ddc'};

% NOTE: if targets are different syls then extracts sametypes that can
% match to either.

%% 1) extract all possible syls (regexp)
% motif_regexpr_str = {}; % will put motifs for extraction here
% for i=1:length(MotifsActual)
%     motif_actual = MotifsActual{i};
%     
%     % for each vocalization in the motif, extract one segment
%     numvocals = length(motif_actual);
%     
%     for ii=1:numvocals
%         if ii==1
%             segmentmotif = ['(' motif_actual(ii) ')' motif_actual(ii+1:end)];
%             
%         elseif ii == numvocals
%             segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')'];
%             
%         else
%             
%             segmentmotif = [motif_actual(1:ii-1) '(' motif_actual(ii) ')' motif_actual(ii+1:end)];
%             
%         end
%         
%         %        disp(segmentmotif);
%         
%         motif_regexpr_str = [motif_regexpr_str segmentmotif];
%         
%     end
% end
% 
%% do

% get lower case for targ
TargSyls_lower = {};
for i=1:length(TargSyls)
targsyl = TargSyls{i};

ind = strfind(targsyl, '(');
lowersyl = targsyl(ind+1);

TargSyls_lower = [TargSyls_lower lowersyl];
end

numsyls = length(motif_regexpr_str);
SameTypeSyls = {};
DiffTypeSyls = {};
SingleSyls = {};
% === for each syl, compare it to the targ syls
for i=1:numsyls   
    syl = motif_regexpr_str{i};
    syllower = syl(strfind(syl, '(')+1);
    
    % --- save syl lower
    SingleSyls = [SingleSyls syllower];
    
    if any(strcmp(syllower, TargSyls_lower)) & ~any(strcmp(syl, TargSyls))
        % i.e. is not target, and is same type
        SameTypeSyls = [SameTypeSyls syl];
        
    elseif ~any(strcmp(syllower, TargSyls_lower))
        DiffTypeSyls = [DiffTypeSyls syl];
        
    end
end

disp(' -------------------- ');
disp(['Targ: '  TargSyls]);
disp(['Same: ' SameTypeSyls]);
disp(['Diff: ' DiffTypeSyls]);


assert((length(TargSyls) + length(SameTypeSyls) + length(DiffTypeSyls)) ...
    ==length(motif_regexpr_str), 'asasdf');



