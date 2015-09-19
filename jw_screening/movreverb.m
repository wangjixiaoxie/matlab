function remov=movreverb(mov,reverb_length)
% remov=movreverb(mov,reverb_length);
% mov is a matlab movie, reverb_length is number of frames for reverb to carry across

remov=mov;

a=1/(reverb_length+1);
for i=1:reverb_length
	steps(i)=(a*i)^(.5);	
end

for f=length(remov):-1:reverb_length+1		%have to go backwards for no-feedforward-style reverb.
	remov(f).cdata=1.5*a*(remov(f).cdata);
	for r=1:reverb_length			%r=1 is most distant frame in delay line so steps(r) is lowest steps value
		remov(f).cdata=remov(f).cdata+1.5*a*steps(r)*(remov(f-reverb_length+r-1).cdata);
	end
end

%now all frames but the first few should have reverb.

