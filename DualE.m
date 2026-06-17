classdef DualE

    properties

        f % Function value
        d % First derivative
        dd % Second derivative

    end

    methods

        %% Constructor

        function y = DualE(a, b, c)

            if nargin == 0 % Zero dual number 
                y.f = 0;
                y.d = 0;
                y.dd = 0;
            elseif nargin == 1 % Constant dual number
                y.f = a;
                y.d = 0;
                y.dd = 0;
            elseif nargin == 2 % Extended dual number associated with the independent variable
                
                if abs(b - 1) < 2*eps('double')
                    y.f = a;
                    y.d = 1;
                    y.dd = 0;
                else
                    error('DualE: wrong parameters')
                end

            else
                y.f = a;
                y.d = b;
                y.dd = c;
            end

        end

        %% SUM

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
                error('NDual:plus:undefined types')
            end

            y = DualE(F, D, DD);

        end

        %% SUBTRACT

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
                error('NDual:minus:undefined types')
            end

            y = DualE(F, D, DD);

        end

        %% MULTIPLICATION

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
                error('DualE:mtimes:undefined types')
            end

            y = DualE(F, D, DD);

        end

        %% Rational function

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
                error('DualE:rdivide:undefined types')
            end

            y = DualE(F, D, DD);
            
        end

        %% Minus sing

        function y = uminus(u)

            y = (-1)*u;

        end

        %% Power function

        function y = mpower(u, v) %u^v

            if isa(u, 'DualE') && isa(v, 'DualE')
                F = NaN;
                D = NaN;
                DD = NaN;
            elseif isa(u, 'DualE') && isa(v, 'double')
                F = u.f^v;
                D = v*u.f^(v-1)*u.d; % Chain rule
                DD = v*(((v-1)*u.f^(v-2)*u.d)*u.d + u.f^(v-1)*u.dd); % Chain rule
            elseif isa(u, 'double') && isa(v, 'DualE')
                F = NaN;
                D = NaN;
                DD = NaN;
            else
                error('DualE:mpower:undefined types')
            end

            y = DualE(F, D, DD);
        
        end

        %% Square root

         function y = sqrt(u)

        	F = sqrt(u.f);
        	D = u.d*(1/(2*sqrt(u.f)));
        	DD = (sqrt(u.f)*u.dd - u.d*u.d/(2*sqrt(u.f))) / (2*u.f);

        	y = DualE(F, D, DD);

         end

         %% Exponential function

         function y = exp(u)

             F = exp(u.f);
             D = u.d*exp(u.f);
             DD = u.dd*exp(u.f) + u.d^2*exp(u.f);

             y = DualE(F, D, DD);

         end

        %% Trigonometric Functions

        function y = sin(u)

            F = sin(u.f);
            D = cos(u.f)*u.d; % Chain rule
            DD = (-sin(u.f)*u.d)*u.d + cos(u.f)*u.dd; % Chain rule
            
            y = DualE(F, D, DD);
            
        end

        function y = cos(u)

            F = cos(u.f);
            D = -sin(u.f)*u.d; % Chain rule
            DD = -(u.dd*sin(u.f) + u.d^2*cos(u.f)); % Chain rule
            
            y = DualE(F, D, DD);
            
        end

    end

end







