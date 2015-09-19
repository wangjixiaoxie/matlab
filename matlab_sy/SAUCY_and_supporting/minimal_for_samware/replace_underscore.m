% replaced '_' with '-'
function out=foo(fname)
id=find(fname=='_');
fname(id)='-';
out=fname;
