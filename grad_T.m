% Parameters
% fn : handler or function
% x0 : coordinates
% Result
% g := grad vector
function [g, v] = grad_T(fn , x0, x_min)

    dim = length(x0); % dim
    g = zeros(1, dim);
    xm = cell(1, dim);

    for i = 1:dim

        xm{i} = x_min(i);
    
    end
    
    for i = 1:dim % For each dimension, calculate the derivative
    
        % Dual Number Vector
    
        x = cell(1, dim);
    
        for j = 1:dim
    
            if j == i
    
                x{j} = DualN(x0(j), 1);
    
            else
    
                x{j} = x0(j);
    
            end
    
        end
        
        f = tunnel_c(fn, x, xm); % Dual number, f.f. value of the function, f.d. value of the derivative
        g(i) = f.d;
    
    end
    
    v = f.f;
    
end
