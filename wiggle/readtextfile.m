function list = readtextfile(listname);
%reads a text file into a 2-D character array
list = [];
fid = fopen(listname,'r');
ln = 1;
while feof(fid) == 0
    list = char(list,fgetl(fid));
    ln = ln+1;
    if rem(ln,1000) == 0
        disp(['Reading line # ', num2str(ln)])
    end
end
list(1,:)=[];
fclose(fid);