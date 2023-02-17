%-------------------------------------------------------------------------------    
%
% Convert the angles to a distance by assuming a cylinder with the mean
% radius of all the radii
%
%-------------------------------------------------------------------------------
function [ts,transobj,r_four] = fourb_trans(rxys,ts,transobj,fwd_flg)


%-------------------------------------------------------------------------------
N = 4;

%-------------------------------------------------------------------------------
% Perform the forward or inverse transform. 
% Forward: goes from angles to an distance variable
% Inverse: goes from the distance variable back to angles
%-------------------------------------------------------------------------------
if fwd_flg == 1
    %---------------------------------------------------------------------------
    % Forward & Initialize: goes from angles to an distance variable. This
    % case also defines the mapping via the radii. 
    [r_four, b]  = four_fit(rxys,ts,4);
    for n = 1:length(ts)
        ts(n) = quad(@(t) four_eval_slen(t,b,N), 0, ts(n));
    end
    transobj.b   = b;
elseif fwd_flg == 2
    %---------------------------------------------------------------------------
    % Forward: goes from angles to an distance variable. This case uses the
    % mean radii from a previous initialized run
    b = transobj.b;
    %     tsin = ts;
    for n = 1:length(ts)
        ts(n) = quad(@(t) four_eval_slen(t,b,N), 0, ts(n));
    end
    % [tsin' ts']
    
elseif fwd_flg == -1
    %---------------------------------------------------------------------------
    % Inverse: goes from the distance variable back to angles
    b    = transobj.b;
    for n = 1:length(ts)
        ts(n) = fminbnd(@(t) four_eval_slen_inv(t,b,N,ts(n)), -2*pi, 2*pi+2*pi);
    end
    r_four = eval_four_fit(ts,4,b);
end