function  row_vect=makerow(vect);

%makes sure that argument is a row vector;


[nrows, ncols]=size(vect);


if min([nrows, ncols]) > 1
  disp('make_row: argument is not a vector');
else
  if ncols > nrows
    row_vect=vect;
  else
    row_vect=vect';
  end
end
