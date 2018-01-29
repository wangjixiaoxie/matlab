function SegExtract =lt_neural_QUICK_GetRawNeural(SegExtract, SummaryStruct, birdnum, neurnum)
%% lt 1/24/18 - extracts raw neural (unfiltered) aligned to data for each trial in segextract

%%



%% go and load neural dat for this seg

[SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, birdnum, neurnum);
cd(SummaryStruct.birds(birdnum).neurons(neurnum).dirname);
tmp = load('data.mat');

%% extract neural dat for period in all trials of SegExtract

for i=1:length(SegExtract)
   ons = SegExtract(i).global_ontime_motifInclFlank;
   off = SegExtract(i).global_offtime_motifInclFlank;
   
   fs = SegExtract(i).fs;
   
   % -- sanity check
   indtmp = SegExtract(i).global_tokenind_DatAlignedToOnsetOfThis;
   
   assert(abs(ons - SongDat.AllOnsets(indtmp))<0.5, 'probably problem...');
   
   
   % -- find corresponding samples
   ons_samp = round(fs*ons);
   off_samp = round(fs*off);
   
   datseg = tmp.data(ons_samp:off_samp);
   datseg = single(datseg);

   SegExtract(i).neural_rawdat = datseg;
  
end

