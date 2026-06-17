function [x_min, f_min, pto] = SteepestDescent(fn, x_ini, tolSD, maxIterSD, delta, maxIterLS, tolX)

       pto = [];

       dim = length(x_ini);
       x = x_ini;
       [G, v] = grad(fn, x);

       iterSD = 0;
       x_ant = x + 10*tolX;
       pasoAceptado = true;

       while (norm(G) > tolSD) && (iterSD < maxIterSD) && (norm(x - x_ant) >= tolX)

           x_ant = x;
           v_ant = v;
           dir = G/norm(G);

           alfaValor = 0;
           alfa = DualE(alfaValor, 1, 0);
           xt = cell(dim, 1);

           for i = 1:dim
               xt{i} = x(i) - alfa*dir(i);
           end

           ft = fn(xt);
           alfaIzq = alfaValor;
           derIzq = ft.d;
           iterN = 0;
           cambio = false;
           paso = delta;
           maxPaso = max(50*delta, 0.05*(1 + norm(x)));
           factorCrecimiento = 1.25;

           while ~cambio && (iterN < maxIterLS)

               alfaDer = alfaValor + paso;
               alfa = DualE(alfaDer, 1, 0);
               xt = cell(dim, 1);

               for i = 1:dim
                   xt{i} = x(i) - alfa*dir(i);
               end

               ft = fn(xt);
               derDer = ft.d;

               if sign(derIzq) ~= sign(derDer)
                   cambio = true;
               else
                   alfaValor = alfaDer;
                   alfaIzq = alfaDer;
                   derIzq = derDer;
                   paso = min(factorCrecimiento*paso, maxPaso);
               end

               iterN = iterN + 1;
               pto = [pto; xt{1}.f xt{2}.f ft.f];

           end

           if cambio

               tolAlfa = 1e-8;
               maxIterBis = 50;

               for iterBis = 1:maxIterBis

                   alfaValor = (alfaIzq + alfaDer)/2;
                   alfa = DualE(alfaValor, 1, 0);
                   xt = cell(dim, 1);

                   for i = 1:dim
                       xt{i} = x(i) - alfa*dir(i);
                   end

                   ft = fn(xt);
                   derMed = ft.d;
                   pto = [pto; xt{1}.f xt{2}.f ft.f];

                   if (abs(derMed) < tolSD) || (abs(alfaDer - alfaIzq) < tolAlfa)
                       break
                   end

                   if sign(derIzq) ~= sign(derMed)
                       alfaDer = alfaValor;
                       derDer = derMed;
                   else
                       alfaIzq = alfaValor;
                       derIzq = derMed;
                   end

               end

           end

           if ft.f < v_ant

               for i = 1:dim
                   x(i) = xt{i}.f;
               end

           else

               x = x_ant;
               pasoAceptado = false;
               break

           end

           [G, v] = grad(fn, x);
           iterSD = iterSD + 1;

       end

       if norm(G) <= tolSD
           fprintf('SteepestDescent stopped due to a small gradient: ||G|| = %.10e\n', norm(G))
       elseif iterSD >= maxIterSD
           fprintf('SteepestDescent stopped due to maximum number of iterations: iterSD = %d\n', iterSD)
       elseif norm(x - x_ant) < tolX
           fprintf('SteepestDescent stopped by a small step: ||x_n - x_{n-1}|| = %.10e\n', norm(x - x_ant))
       elseif ~pasoAceptado
           fprintf('SteepestDescent stopped because it could not find an acceptable downward path.\n')
       else
           fprintf('SteepestDescent stopped due to an unidentified condition.\n')
       end

       x_min = x;
       f_min = v;

end
