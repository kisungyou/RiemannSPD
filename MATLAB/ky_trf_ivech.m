function mat = ky_trf_ivech(vec)

% KY_TRF_IVECH takes elements from a given vector and stacks it onto the main
% diagonal and upper triangular part of a square matrix. Then, it fills the
% lower triangular part in a symmetric manner. 
%
%   * USAGE
%       mat = KY_TRF_IVECH(vec)
%
%   * INPUT
%       vec    a column vector of length n*(n+1)/2
%
%   * OUTPUT
%       mat    an (n-by-n) symmetric matrix
%
%   * AUTHOR   Kisung You (kyou@nd.edu)
%   * HISTORY
%       0.1. [06/2022] initial implementation.
%
%   See also KY_TRF_VECH


%% preprocessing
%   1. must be a vector
if (~isvector(vec))
    error('* ky_trf_ivech : an input must be a vector');
end

%% main computation
k = length(vec);
n = round(sqrt(2*k + 0.25)-0.5);
mat = zeros(n,n);

mat(triu(ones(n))==1)=vec;
mat = mat + (mat');
mat = mat - diag(diag(mat))/2;
end