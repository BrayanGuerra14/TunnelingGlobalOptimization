classdef DualN

    properties

        f % Valor de la función
        d % Valor de la derivada

    end

    methods

        %% Constructor

        function y = DualN(a, b)

            if nargin == 0 % Número dual 0
                y.f = 0;
                y.d = 0;
            elseif nargin == 1 % Número dual constante
                y.f = a;
                y.d = 0;
            else
                y.f = a;
                y.d = b;
            end

        end

        %% CastDualN a Vector (to vector, 2 vector)

        function y = DN2V(u)
            
            y = zeros(1, 2);
            y(1) = u.f;
            y(2) = u.d;

        end

        %% Suma

        function y = plus(u, v)

            if isa(u, 'DualN') && isa(v, 'DualN')
                F = u.f + v.f;
                D = u.d + v.d;
            elseif isa(u, 'DualN') && isa(v, 'double')
                F = u.f + v;
                D = u.d;
            elseif isa(u, 'double') && isa(v, 'DualN')
                F = u + v.f;
                D = v.d;
            else 
                error('NDual:plus:Tipos no definidos')
            end

            y = DualN(F, D);

        end

        %% Resta

        function y = minus(u, v)

            if isa(u, 'DualN') && isa(v, 'DualN')
                F = u.f - v.f;
                D = u.d - v.d;
            elseif isa(u, 'DualN') && isa(v, 'double')
                F = u.f - v;
                D = u.d;
            elseif isa(u, 'double') && isa(v, 'DualN')
                F = u - v.f;
                D = -v.d;
            else 
                error('NDual:minus:Tipos no definidos')
            end

            y = DualN(F, D);

        end

        %% Multiplicación

        function y = mtimes(u, v)

            if isa(u, 'DualN') && isa(v, 'DualN')
                F = u.f * v.f;
                D = u.f * v.d + u.d * v.f;
            elseif isa(u, 'DualN') && isa(v, 'double')
                F = u.f * v;
                D = u.d * v;
            elseif isa(u, 'double') && isa(v, 'DualN')
                F = u * v.f;
                D = u * v.d;
            else 
                error('NDual:mtimes:Tipos no definidos')
            end

            y = DualN(F, D);

        end

        %% División

        function y = mrdivide(u, v)

            if isa(u, 'DualN') && isa(v, 'DualN')
                F = u.f / v.f;
                D = (u.d * v.f - u.f * v.d)/ v.f^2;
            elseif isa(u, 'DualN') && isa(v, 'double')
                F = u.f / v;
                D = u.d / v;
            elseif isa(u, 'double') && isa(v, 'DualN')
                F = u / v.f;
                D = - u * v.d / v.f^2;
            else
                error('DualN:mrdivide:Tipos no definidos')
            end

            y = DualN(F, D);

        end
        
        %% Cambio de signo (operador unario que cambia de signo)
        
        function y = uminus(u)

            y = (-1)*u;

        end

        %% Potenciación

        function y = mpower(u, v)

            if isa(u, 'DualN') && isa(v, 'DualN')
                F = NaN;
                D = NaN;
            elseif isa(u, 'DualN') && isa(v, 'double')
                F = u.f^v;
                D = v*u.f^(v-1)*u.d; % Regla de la cadena
            elseif isa(u, 'double') && isa(v, 'DualN')
                F = NaN;
                D = NaN;
            else
                error('NDual:mpower:Tipos no definidos')
            end

            y = DualN(F, D);

        end

        %% Radicación

        function y = sqrt(u)

        	F = sqrt(u.f);
        	D = u.d*(1/(2*sqrt(u.f)));
        	y = DualN(F, D);

        end

        %% Función exponencial

        function y = exp(u)

            F = exp(u.f);
            D = u.d*exp(u.f);

            y = DualN(F, D);

        end

        %% Funciones trigonométricas

        function y = sin(u)

            F = sin(u.f);
            D = cos(u.f)*u.d; % Derivada con regla de la cadena
            
            y = DualN(F, D);

        end

        function y = cos(u)

            F = cos(u.f);
            D = -sin(u.f)*u.d; % Regla de la cadena

            y = DualN(F, D);
            
        end
        
    end

end


%% Revisar Peñuñuri y Neidinger