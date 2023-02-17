%-------------------------------------------------------------------------------
%
% Compute the pairwise distance (squared) matrix between each sample. We assume
% the data is in the form (number of sample) x (number of dimension/features)
%
%-------------------------------------------------------------------------------
function D = comp_pairwise_distmat(X)



% Compute pairwise distance matrix
% tic
sum_X = sum(X .^ 2, 2);
D     = bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (X * X')));
% toc

%-------------------------------------------------------------------------------
% Checking
%-------------------------------------------------------------------------------
% 1000 point data
% Above code: Elapsed time is 0.007676 seconds.
% Below code: Elapsed time is 0.608363 seconds.
% Factor improvement: 79 
%-------------------------------------------------------------------------------
% 2000 point data
% Above code: Elapsed time is 0.032209 seconds.
% Below code: Elapsed time is 2.296322 seconds.
% Factor improvement: 71
%-------------------------------------------------------------------------------
% 4000 point data
% Above code: Elapsed time is 0.109497 seconds.
% Below code: Elapsed time is 9.464961 seconds.
% Factor improvement: 86
%-------------------------------------------------------------------------------
% tic
% D0 = zeros(size(X,1),1);
% for n = 1:size(X,1)
%     for m = 1:size(X,1)
%         D0(n,m) = norm( X(n,:) - X(m,:) )^2;
%     end
% end
% toc
% norm( D(:) - D0(:))