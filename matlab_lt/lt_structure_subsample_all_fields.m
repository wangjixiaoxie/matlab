function Struct_out=lt_structure_subsample_all_fields(Struct_in, inds)
%% LT 9/23/15 - given a structure with multiple fields, subsample indices of all the fields (each field can be vector or cell array
% NOTE: max of inds must be <= field with minimum numel

% e.g.
% == INPUT:
% Struct_in = 
% 
%     x: [1x100 double]
%     y: [1x100 double]
%     z: {[1]  'asdf'  '3'  [4]  [4]  'asdfasd'}
%     
% inds=1:5;
% 
% 
% == Output:
% Struct_out=
%     x: [1 2 3 4 5]
%     y: [1 1.4949494949495 1.98989898989899 2.48484848484848 2.97979797979798]
%     z: {[1]  'asdf'  '3'  [4]  [4]}



%% RUN

% === make function handle to take subset of inds
func=@(x)x(inds);

% === get subset of all fields
Struct_out=structfun(func, Struct_in,'UniformOutput',0);




    
