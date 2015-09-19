%% LT 4/19/15 - takes multiple batches - sorts each by label vs. unlabeled.
% enter each batch as a string, a different variable.

function lt_sort_batch_by_labeled_MULT(varargin)

NumBatches=nargin;

for i=1:NumBatches;

    disp(varargin{i});
    lt_sort_batch_by_labeled(varargin{i});
    
end

disp('DONE!');


