
%% Tarin Ziyaee / tarin.ziyaee@gmail.com
%  Recursive Extrema Detector (RED)

%% Function Summary
%  This function will find all local extrema of a signal that satisfy a
%  given recursive threshold. (i.e. Local mins and maxes).

%% Inputs
% x           : The signal vector
% T           : Recursive threshold (between 0 and 1)
% varargin{1} : 'normal': Will excecute in normal mode (default value if left empty)
%               'debug' : Will excecute in debug mode for vizualization

%% Outputs
% dmaxes      : Indices of all local maxima
% dmins       : Indices of all local minima

function [dmaxes, dmins] = RED(x,T, varargin)

    %Make sure program doesn't recursively explode...
    if T <= 0 error('T must be positive'); end
        
    %Initialize data and scanner
    x = x(:).'; x = x./(max(x)-min(x)); N.x = length(x);
    k.max = 1;  k.min = 1;
    detected.max = zeros(2,N.x); detected.min = zeros(2,N.x);
    count.max = 0; count.min = 0;
    

    %Initialize Scanner State Machine
    state.maxScanner = 1;
    state.minScanner = 1;
    state.cold = 1;

    %Check if normal or debug mode
    V = length(varargin);
    if V == 0
        varargin{1} = 'normal';
    end
    
    switch varargin{1}
        case 'normal' 
            
            %Begin Recursive Detection
            for kk = 2:N.x
               if state.maxScanner
                   dy = x(kk) - x(k.max);
                   if (dy < -T)
                       count.max = count.max + 1;
                       detected.max(:,count.max) = [k.max ; dy];
                       state.maxScanner = 0; 
                       state.minScanner = 1;
                       k.min = kk;
                       if state.cold
                           detected.max(1) = 0;
                           count.max = 0;
                           state.cold = 0;
                       end                      
                   elseif (x(kk) > x(k.max)) 
                       k.max = kk;                         
                   end       
               end
               if state.minScanner
                   dy = x(kk) - x(k.min);
                   if (dy > T)
                       count.min = count.min + 1;
                       detected.min(:,count.min) = [k.min ; dy];
                       state.minScanner = 0;
                       state.maxScanner = 1;
                       k.max = kk;
                       if state.cold
                           detected.min(1) = 0;
                           count.min = 0;
                           state.cold = 0;
                       end                      
                   elseif (x(kk) < x(k.min))
                       k.min = kk;                     
                   end
               end
            end
            dmaxes = detected.max(:,1:count.max);
            dmins = detected.min(:,1:count.min);
          
        case 'debug'
            
            %Plot for visualization/debug
            figure(1); clf(1); plot(x); grid on; hold on;
            
            %Begin Recursive Detection
            for kk = 2:N.x
               if state.maxScanner
                   dy = x(kk) - x(k.max);
                   if (dy < -T)
                       count.max = count.max + 1;
                       detected.max(:,count.max) = [k.max ; dy];
                       state.maxScanner = 0; 
                       state.minScanner = 1;
                       k.min = kk;
                       if state.cold
                           detected.max(1) = 0;
                           count.max = 0;
                           state.cold = 0;
                       end                      
                   elseif (x(kk) > x(k.max)) 
                       k.max = kk;
                       plot(k.max, x(k.max), '.g');               
                   end       
               end
               if state.minScanner
                   dy = x(kk) - x(k.min);
                   if (dy > T)
                       count.min = count.min + 1;
                       detected.min(:,count.min) = [k.min ; dy];
                       state.minScanner = 0;
                       state.maxScanner = 1;
                       k.max = kk;
                       if state.cold
                           detected.min(1) = 0;
                           count.min = 0;
                           state.cold = 0;
                       end                      
                   elseif (x(kk) < x(k.min))
                       k.min = kk;
                       plot(k.min, x(k.min), '.r');
                   end
               end
            end
            dmaxes = detected.max(:,1:count.max);
            dmins = detected.min(:,1:count.min);
            for ii = 1:length(dmaxes)
                 plot(dmaxes(1,:), x(dmaxes(1,:)), '.g', 'MarkerSize', 20);
            end
            for ii = 1:length(dmins)
                plot(dmins(1,:), x(dmins(1,:)), '.r','MarkerSize', 20);
            end
%             axis([0 N.x min(x)-sign(min(x))*0.1*min(x) max(x)+sign(max(x))*0.1*max(x)]);
            ylim([-1 1]);
    otherwise 
        error('Last argument not valid.');
    end
end

