function vec = ky_trf_vech(mat)

% KY_TRF_VECH extracts the upper triangular part as well as the main diagonal
% part of a matrix and gathers the elements as a column vector. The
% function performs a job in correspondance with the wikipedia page.
%
%   * USAGE
%       vec = KY_TRF_VECH(mat)
%
%   * INPUT
%       mat    an (n-by-n) matrix
%
%   * OUTPUT
%       vec    a column vector of length n*(n+1)/2
%
%   * AUTHOR   Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [06/2022] initial implementation.
%
%   See also KY_TRF_IVECH

%% preprocessing
%   1. should be a matrix
if (~ismatrix(mat))
    error("* ky_trf_vech : an input must be a matrix.");
end
%   2. square matrix
if (size(mat,1)~=size(mat,2))
    error("* ky_trf_vech : an input must be a square matrix.");
end

%% main part
mask = triu(true(size(mat)));
vec  = mat(mask);
end