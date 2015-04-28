function [x,y,rsqrd] = fit_circle(P1, P2, P3)

    A = [ P1 P2 P3; 1 1 1 ]';
    b = [ -sum(P1.^2) -sum(P2.^2) -sum(P3.^2) ]';
    
    s = A \ b;

    x = -s(1)/2;
    y = -s(2)/2;
    
    rsqrd = -s(3) + x*x + y*y;
    
endfunction