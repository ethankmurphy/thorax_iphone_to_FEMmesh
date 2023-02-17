%-------------------------------------------------------------------------------
%
% Fit the data to an ellipse
%
%-------------------------------------------------------------------------------
function r_four = eval_four_fit(thts,N,b)

%-------------------------------------------------------------------------------
% Construct the data matrix
A = zeros(length(thts),2*N+1);
for m = 1:length(thts)
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
% Construct the fit data
r_four   = A*b;
% four_dat = [ ...
%     r_four.*cos(thts)+xy0(1) ...
%     r_four.*sin(thts)+xy0(2)];


