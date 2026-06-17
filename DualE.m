classdef DualE

    properties

        f % Valor función
        d % Valor primera derivada
        dd % Valor segunda derivada

    end

    methods

        %% Constructor

        function y = DualE(a, b, c)

            if nargin == 0 % Número dual extendido cero
                y.f = 0;
                y.d = 0;
                y.dd = 0;
            elseif nargin == 1 % Número dual extendido constante
                y.f = a;
                y.d = 0;
                y.dd = 0;
            elseif nargin == 2 % Número dual extendido relacionado con la variable independiente
                
                if abs(b - 1) < 2*eps('double')
                    y.f = a;
                    y.d = 1;
                    y.dd = 0;
                else
                    error('DualE: parámetros incorrectos')
                end

            else
                y.f = a;
                y.d = b;
                y.dd = c;
            end

        end

        %% SUMA

        function y = plus(u, v)

            if isa(u, 'DualE') && isa(v, 'DualE')
                F = u.f + v.f;
                D = u.d + v.d;
                DD = u.dd + v.dd;
            elseif isa(u, 'DualE') && isa(v, 'double')
                F = u.f + v;
                D = u.d;
                DD = u.dd;
            elseif isa(u, 'double') && isa(v, 'DualE')
                F = u + v.f;
                D = v.d;
                DD = v.dd;
            else
                error('NDual:plus:tipos no definidos')
            end

            y = DualE(F, D, DD);

        end

        %% RESTA

        function y = minus(u, v)

            if isa(u, 'DualE') && isa(v, 'DualE')
                F = u.f - v.f;
                D = u.d - v.d;
                DD = u.dd - v.dd;
            elseif isa(u, 'DualE') && isa(v, 'double')
                F = u.f - v;
                D = u.d;
                DD = u.dd;
            elseif isa(u, 'double') && isa(v, 'DualE')
                F = u - v.f;
                D = -v.d;
                DD = -v.dd;
            else
                error('NDual:minus:tipos no definidos')
            end

            y = DualE(F, D, DD);

        end

        %% MULTIPLICACIÓN

        function y = mtimes(u, v)

            if isa(u, 'DualE') && isa(v, 'DualE')
                F = u.f*v.f;
                D = u.f*v.d+u.d*v.f;
                DD = u.f*v.dd + u.d*v.d + u.d*v.d + u.dd*v.f;
            elseif isa(u, 'DualE') && isa(v, 'double')
                F = u.f*v;
                D = u.d*v;
                DD = u.dd*v;
            elseif isa(u, 'double') && isa(v, 'DualE')
                F = u*v.f;
                D = u*v.d;
                DD = u*v.dd;
            else
                error('DualE:mtimes:tipos no definidos')
            end

            y = DualE(F, D, DD);

        end

        %% División

        function y = mrdivide(u, v)

            if isa(u, 'DualE') && isa(v, 'DualE')
                F = u.f / v.f;
                D = (u.d * v.f - u.f * v.d)/ v.f^2;
                DD = ((u.dd * v.f + u.d * v.d - u.d * v.d - u.f * v.dd) * v.f^2 - (u.d * v.f - u.f * v.d) * 2 * v.f * v.d) / v.f^4;
            elseif isa(u, 'DualE') && isa(v, 'double')
                F = u.f / v;
                D = u.d / v;
                DD = u.dd / v;
            elseif isa(u, 'double') && isa(v, 'DualE')
                F = u / v.f;
                D = - u * v.d / v.f^2;
                DD = (- u * v.dd * v.f^2 + u * v.d * 2 * v.f * v.d) / v.f^4;
            else
                error('DualE:rdivide:Tipos no definidos')
            end

            y = DualE(F, D, DD);
            
        end

        %% Cambio de signo

        function y = uminus(u)

            y = (-1)*u;

        end

        %% Potenciación

        function y = mpower(u, v) %u^v

            if isa(u, 'DualE') && isa(v, 'DualE')
                F = NaN;
                D = NaN;
                DD = NaN;
            elseif isa(u, 'DualE') && isa(v, 'double')
                F = u.f^v;
                D = v*u.f^(v-1)*u.d; % Regla de la cadena
                DD = v*(((v-1)*u.f^(v-2)*u.d)*u.d + u.f^(v-1)*u.dd); % Regla de la cadena
            elseif isa(u, 'double') && isa(v, 'DualE')
                F = NaN;
                D = NaN;
                DD = NaN;
            else
                error('DualE:mpower:tipos no definidos')
            end

            y = DualE(F, D, DD);
        
        end

        %% Radicación

         function y = sqrt(u)

        	F = sqrt(u.f);
        	D = u.d*(1/(2*sqrt(u.f)));
        	DD = (sqrt(u.f)*u.dd - u.d*u.d/(2*sqrt(u.f))) / (2*u.f);

        	y = DualE(F, D, DD);

         end

         %% Función exponencial

         function y = exp(u)

             F = exp(u.f);
             D = u.d*exp(u.f);
             DD = u.dd*exp(u.f) + u.d^2*exp(u.f);

             y = DualE(F, D, DD);

         end

        %% Funciones trigonométricas

        function y = sin(u)

            F = sin(u.f);
            D = cos(u.f)*u.d; % Regla de la cadena
            DD = (-sin(u.f)*u.d)*u.d + cos(u.f)*u.dd; % Regla de la cadena
            
            y = DualE(F, D, DD);
            
        end

        function y = cos(u)

            F = cos(u.f);
            D = -sin(u.f)*u.d; % Regla de la cadena
            DD = -(u.dd*sin(u.f) + u.d^2*cos(u.f)); % Regla de la cadena
            
            y = DualE(F, D, DD);
            
        end

    end

end







