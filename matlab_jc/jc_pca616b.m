function m=jc_pca616b(input)


%put the data into a matrix for pca
b=size(input);
for i=1:401
    m(i)=mean(input(i,:));
end
