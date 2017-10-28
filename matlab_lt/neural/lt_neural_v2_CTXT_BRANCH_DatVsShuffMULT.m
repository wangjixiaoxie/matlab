function lt_neural_v2_CTXT_BRANCH_DatVsShuffMULT(allanalyfnames)

%% lt 10/27/17 - takes multiple analyses and plots
% INPUT fnames
% allanalyfnames = {...
%     'xaa_Algn2Ons1_26Oct2017_1257_testLMAN2birds', ...
%     'xaa_Algn2Ons1_26Oct2017_1257_testLMAN2birds', ...
%     };


%%

DATSTRUCT = struct;

for i=1:length(allanalyfnames)
    
    afname = allanalyfnames{i};
    dstruct = lt_neural_v2_CTXT_BRANCH_DatVsShuffPLOT(afname);
    
    % ====== save structs
    DATSTRUCT.analynum(i).dat = dstruct;
    DATSTRUCT.analynum(i).fname = afname;
    
%     if ~exist('DATSTRUCT', 'var')
%         DATSTRUCT = dstruct;
%     else
%         DATSTRUCT = [DATSTRUCT; dstruct];
%     end
    
end

%% compare two analyses

lt_figure; hold on;
numanalyses = length(DATSTRUCT.analynum);

Yvals ={};
Fnames = {};
for i=1:numanalyses
    
    yvals = DATSTRUCT.analynum(i).dat.AllDecode_z;
    fname = DATSTRUCT.analynum(i).fname;
    
    Yvals = [Yvals yvals];
    Fnames = [Fnames fname];
    
end

lt_plot_MultDist(Yvals, 1:length(Yvals), 1, 'k', 0, 0);
% set(gca, 'XTickLabel', Fnames);
% rotateXLabels(gca, 45)
ylabel('decode (z)');
xlabel('analysis num');
lt_plot_zeroline


%% ============== 
