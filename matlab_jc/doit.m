function [fbd,pitchbd,fbud,pitchbud,fbudd,pitchbudd]=doit
[fbd,pitchbd]=jc_PitchData611(1024,1020,1,1900,2800,'batchD','d','d','d');
[fbud,pitchbud]=jc_PitchData611(1024,1020,1,1900,2800,'batchUD','d','d','d');
[fbudd,pitchbudd]=jc_PitchData611(1024,1020,1,1900,2800,'batchUDd','d','d','d');
a=1;
