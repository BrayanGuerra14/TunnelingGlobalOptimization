% Parámetros
% fn : handler o función
% x0 : coordenadas
% Resultado
% g := vector gradiente
function [g, v] = grad_T(fn , x0, x_min)

    dim = length(x0); % Número de dimensiones
    g = zeros(1, dim);
    xm = cell(1, dim);

    for i = 1:dim

        xm{i} = x_min(i);
    
    end
    
    for i = 1:dim % Para cada dimensión calcular la derivada
    
        % Vector de números duales
    
        x = cell(1, dim);
    
        for j = 1:dim
    
            if j == i
    
                x{j} = DualN(x0(j), 1);
    
            else
    
                x{j} = x0(j);
    
            end
    
        end
        
        f = tunnel_c(fn, x, xm); % Número dual, f.f valor de la función, f.d valor de la derivada
        g(i) = f.d;
    
    end
    
    v = f.f;
    
end