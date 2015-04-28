%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_grains(grains, brightness, contrast)
    
    plot_matrix = zeros(800,1032);
    ranges = [1024 799 64 256];

    for i = 1:size(grains,2)
        mins(i) = min(grains(:,i));
        maxs(i) = max(grains(:,i));
        ords(:,i) = ...
            (grains(:,i) - mins(i)) / (maxs(i) - mins(i)) * ranges(i) + 1;
    endfor
    for i = 1:size(grains,1)
        plot_matrix(801-round(ords(i,2)), [round(ords(i,1)):round(ords(i,1)+ords(i,4))-1]) = ...
        plot_matrix(801-round(ords(i,2)), [round(ords(i,1)):round(ords(i,1)+ords(i,4))-1]) + ...
        abs(ords(i,3));
    endfor
    
    plot_matrix = abs(plot_matrix);
    plot_matrix = plot_matrix - min(plot_matrix(:));
    plot_matrix = plot_matrix / max(plot_matrix(:));
    image(contrast * (log(plot_matrix) + brightness));
    
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
