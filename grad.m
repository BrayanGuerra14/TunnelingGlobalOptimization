% Parameters
% fn : handler or function
% x0 : coordinates
% Result
% g := grad vector
function [g, v] = grad(fn , x0)

dim = length(x0); % dim

g = zeros(1, dim);

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

    f = fn(x); % Dual number, f.f. value of the function, f.d. value of the derivative
    g(i) = f.d;

end

v = f.f;

end
