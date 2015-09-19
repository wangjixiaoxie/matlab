function [out]=findh(list)
% 
%list is a vector that is a list of integers. findh picks out which of these numbers forms a set of harmonics: numbers that are integer multiples of a single number.

for i=1:length(list)
	test(i,:)=list./list(i);
end

result=0.*test;

for i=1:numel(test)
%       if vpa(test(i),3)==vpa(fix(test(i)),3);
	if abs(test(i)-round(test(i)))<.026*test(i)
		result(i)=1;
	end
end

% now there is the business of L-shapes through results to find pitches
% and scoring those pitches by how many harmonics are present... and
% how loud they were in the initial signal?

test
result

return
