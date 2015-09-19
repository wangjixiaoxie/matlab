function [fcbd,pitchcbd,fcbud,pitchcbud,fcbudd,pitchcbudd]=jc_do
[fcbd,pitchcbd]=jc_PitchData611(1024,1020,1,3000,4000,'batchD','c','c','c');
[fcbud,pitchcbud]=jc_PitchData611(1024,1020,1,3000,4000,'batchUD','c','c','c');
[fcbudd,pitchcbudd]=jc_PitchData611(1024,1020,1,3000,4000,'batchUDd','c','c','c');
