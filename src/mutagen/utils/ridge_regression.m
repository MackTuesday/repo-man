
function coeffs = ridge_regression(x, y, lambda)

    covariances = x' * x;
    lambda_i = lambda * eye(size(covariances));
    coeffs = (covariances + lambda_i) \ (x' * y);

    % To speed things up, do
    % u = cholesky_decomposition(x);
    % Now u is upper triangular and easy to invert,
    % or you can use elimination directly:
    % First substitute u'z = x'y, solve for z, then
    % u beta = z, solve for beta

endfunction