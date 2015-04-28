
function coeffs = least_squares(x, y)

    dot_products = x' * x;
    % This preconditioner doesn't seem to do much good...
    preconditioner = diag(1 ./ norm(dot_products, 2, "cols"));
    coeffs = (preconditioner * dot_products) \ (preconditioner * x' * y);

    % To speed things up, do
    % u = cholesky_decomposition(x);
    % Now u is upper triangular and easy to invert,
    % or you can use elimination directly:
    % First substitute u'z = x'y, solve for z, then
    % u beta = z, solve for beta

endfunction