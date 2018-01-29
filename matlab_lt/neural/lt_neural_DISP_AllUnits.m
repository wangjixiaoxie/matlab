function lt_neural_DISP_AllUnits(SummaryStruct)

numbirds = length(SummaryStruct);

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(birdname);
    
    exptnames = unique({SummaryStruct.birds(i).neurons.exptID});
    
    for ii=1:length(exptnames)
        
        neurlist = find(strcmp({SummaryStruct.birds(i).neurons.exptID}, exptnames{ii}));
        
        for iii=neurlist
            
            loc = SummaryStruct.birds(i).neurons(iii).NOTE_Location;
            plot(ii, iii, 'or');
            lt_plot_text(ii, iii, [num2str(iii) '-' loc], 'k');
        end
        
    end
    grid on;
    set(gca, 'XTickLabel', exptnames);
    rotateXLabels(gca, 90);
    
end

