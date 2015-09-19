function get_rand_nt_pb(batch,NT,PRET,PSTT,NNotes,TSilence,CS);
% get_rand_nt_pb(batch,NT,PRET,PSTT,NNotes,TSilence,CS);
% all times are in ms



if (~exist('TSilence'))
	TSilence=10;
elseif (length(TSilence)==0)
	TSilence=10;
end
scil=zeros([ceil(TSilence*44.1),1]);

if (~exist('NT'))
	NT='a';
elseif (length(NT)==0)
	NT='a';
end

if (~exist('NNotes'))
	NNotes=20;
elseif (length(NNotes)==0)
	NNotes=20;
end

if (~exist('CS'))
	CS='obs0';
elseif (length(CS)==0)
	CS='obs0';
end

ff=load_batchf(batch);
inds=randperm(length(ff));

sndf=[];
ntfnd=0;ii=0;
while (ntfnd<NNotes)
	if (ii>=length(inds))
		break;
	end
	pp=[];
	while (length(pp)==0)
		ii=ii+1;
		fn=ff(inds(ii)).name;
		load([fn,'.not.mat']);
		pp=findstr(labels,NT);
	end
	disp(fn);
	[dat,fs]=evsoundin('',fn,CS);

	% pick one of the notes randomly
	tmp=randperm(length(pp));
	pp=pp(tmp(1));
	datinds=[floor((onsets(pp)-PRET)*1e-3*fs):ceil((offsets(pp)+PSTT)*1e-3*fs)];
	dattmp=dat(datinds);
	mask=mk_mask(dattmp,fs,5e-3,length(dattmp)/fs-5e-3,1e-3);
	dattmp=dattmp.*mask;
	dattmp=resample(dattmp,44100,fix(fs));
	dattmp=dattmp./max(abs(dattmp))./1.01;
	sndf=[sndf;scil;dattmp];

	ntfnd=ntfnd+1;
end


wavwrite(sndf,44100,16,'TEMP.wav');
return;
