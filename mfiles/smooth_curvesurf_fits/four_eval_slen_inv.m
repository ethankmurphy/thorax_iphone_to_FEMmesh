%-------------------------------------------------------------------------------
%
% Fit the data to an ellipse
%
%-------------------------------------------------------------------------------
function [ds]  = four_eval_slen_inv(ts,b,N,s0)

%-------------------------------------------------------------------------------
% For each angle calculate the arc length
ss = 0*ts;
for n = 1:length(ts)
    ss(n) = quad(@(t) four_eval_slen(t,b,N), 0, ts(n));
end

ds = abs(ss - s0);

