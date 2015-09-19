cleandir4('batch',1e5,500,5,5);
    % tfs is the threshold for segmentation
mk_rmdata('batch630',1)
!csh yes | ./rmdata
    % do it for .tmp,.rec as well
    
    
%%%%%%  Linux Commands   %%%%%%
% cp -r ampon_1031 /cardinal1/...   % -r copies entire directory
% smbmount //plover/jcharles plovmount/
% smbumount plovmount/