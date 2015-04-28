
function [mu, residual] = LeastAngle(in, column_vectors)

    % This implementation is incomplete!
    
    in = in(:);
    
    % Augment active set
    inner_products = inv(column_vectors' * column_vectors);
    normalizer = 1 / sqrt(ones(1, size(active_set_indices,1)) * inner_products * ...
                          ones(size(active_set_indices,1), 1));
    weights = normalizer * inner_products * ones(size(active_set_indices,1));
    equiangular = active_set * weights;
    dot_products_with_equiangular = column_vectors' * equiangular;
    step_candidates_left  = (-little_c_hat + big_c_hat) ./ ...
                            (-dot_products_with_equiangular + normalizer);
    step_candidates_right = ( little_c_hat + big_c_hat) ./ ...
                            ( dot_products_with_equiangular + normalizer);
    pos_step_cands_left =  step_candidates_left  .* (step_candidates_left  > 0);
    pos_step_cands_right = step_candidates_right .* (step_candidates_right > 0);
    [step_size, new_vector_index] = min(min(pos_step_cands_left, pos_step_cands_right));

endfunction