function LogXY = ky_spd_log(X, Y)

% KY_SPD_LOG computes the logarithmic map of Y at X.
%
%   * USAGE
%       LogAB = KY_SPD_LOG(X, Y)
%
%   * INPUT
%       X      an (n-by-n) SPD matrix
%       Y      an (n-by-n) SPD matrix
%
%   * OUTPUT
%       LogXY  an (n-by-n) symmetric matrix
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [06/2022] initial implementation.

% logarithmic map
symm  = @(x) .5*(x+x');
LogXY = symm(real(X*real(logm(X\Y))));

end