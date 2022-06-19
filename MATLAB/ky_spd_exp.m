function Y = ky_spd_exp(X, eta, t)

% KY_SPD_LOG computes the exponential map at X in the direction of eta.
%
%   * USAGE
%       Y  = KY_SPD_EXP(X, eta, t)
%
%   * INPUT
%       X      an (n-by-n) SPD matrix
%       eta    an (n-by-n) symmetric matrix
%       t      the amount of marching (optional; default=1.0)
%
%   * OUTPUT
%       Y      an (n-by-n) symmetric matrix
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [06/2022] initial implementation.

% exponential mapping

if (nargin < 3)
    t = 1.0;
end
symm = @(x) .5*(x+x');
Y    = symm(real(X*real(expm(X\(t*eta)))));

end