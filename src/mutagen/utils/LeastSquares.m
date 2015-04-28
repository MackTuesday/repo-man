
function coeffs = LeastSquares(x, y)

    column_dot_products = x' * x;
    coeffs = inv(column_dot_products) * x' * y;
    
    % To speed things up, do
    % u = cholesky_decomposition(x);
    % Now u is upper triangular and easy to invert,
    % or you can use elimination directly:
    % First substitute u'z = x'y, solve for z, then
    % u beta = z, solve for beta

endfunction