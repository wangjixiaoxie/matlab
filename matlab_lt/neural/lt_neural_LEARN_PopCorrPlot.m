
NumBirds = length(SwitchStruct.bird);

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds
    
   numexpts = length(SwitchStruct.bird(i).exptnum);
   
   assert(length(MOTIFSTATS_pop.birds(i).exptnum) == numexpts, 'asdfsd');
   
   for ii=1:numexpts
       
      numsw = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
      
      for iii=1:numsw
         
          
          %% =========== sanity check - plot all neuron sets for this switch
          
          Datswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
          
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(['b' num2str(i) '-e' num2str(ii) '-sw' num2str(iii)]);
ylabel('neuron set (bk = in this sw)');
xlabel('time of song');

          % -------- go thru all sets. if time of set overlaps this switch,
          % then plot
          
          numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum);
          
          for ss = 1:numsets
             
              songfiles = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_songfiles{ss};
              songfiles = datenum(songfiles, 'yymmdd_HHMMSS');
              neurinset = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
              
              if any(songfiles<Datswitch.switchdnum & songfiles>Datswitch.switchdnum_previous) ...
                  & any(songfiles>Datswitch.switchdnum & songfiles<Datswitch.switchdnum_next)
                 % Then this set has data for this switch
                 
                 plot(songfiles, ss, 'ok');
              else
                  plot(songfiles, ss, 'o', 'Color', [0.7 0.7 0.7]);
              end
                 lt_plot_text(min(songfiles), ss+0.2, ['n' num2str(neurinset)], 'b')
              
          end
          
          ylim([0 numsets+1]);
          line([Datswitch.switchdnum Datswitch.switchdnum], ylim, 'Color', 'r');
          
          
          %% ========== PLOT LEARNING
          
          
           %% ========== COLLECT LEARNING
           
          
          
          
      end
       
       
   end
   
    
end