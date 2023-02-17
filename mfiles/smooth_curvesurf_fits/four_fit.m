%--------------------------------------------------------------------------
%
% Fit the data to an ellipse
%
%--------------------------------------------------------------------------
function [r_four, b]  = four_fit(rs,thts,N)

%--------------------------------------------------------------------------
% Construct the data matrix
A = zeros(size(rs,1),2*N+1);
for m = 1:size(rs,1)
    A(m, 1 ) = 1;
    for k = 1:2        
        for n = 1:N
            if k == 1
                A(m,1+n)   = cos(n*thts(m));
            else
                A(m,1+n+N) = sin(n*thts(m));
            end            
        end
    end
end

%-------------------------------------------------------------------------------
% Calculate the best fit
% ATA  = A'*A;
% lamb = 1e-3;
% b    = pinv(ATA + lamb*eye(size(ATA,1)) )*A'*rs;
b    = inv(A'*A)*A'*rs;

%-------------------------------------------------------------------------------
% Construct the fit data
r_four   = A*b;
% four_dat = [ ...
%     r_four.*cos(thts)+xy0(1) ...
%     r_four.*sin(thts)+xy0(2)];


