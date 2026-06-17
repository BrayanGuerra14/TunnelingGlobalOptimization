clc
clear all

rng('shuffle')

%% Data


% Objective function

% fn = @ackley_n; % (0, 0)
% fn = @three_hump_n; % (0, 0)
% fn = @levy_n13; % (1, 1)
fn = @drop_wave; % (0, 0)
% fn = @esfera; 

% Initial conditions

% X0 = [-1.5, 1];
X0 = 5*rand(1, 2);
X0_i = X0;
X0_C = {X0_i(1), X0_i(2)};

%%

TRA = {}; % Initialize the TRA variable, which stores the search trajectories, in the function

LOG = [];

maxTotalIter = 45; % Maximum number of iterations for the global algorithm (local minimization and tunneling)

iterTotal = 0; % The variable iterTotal is initialized to 0         

tolX = 1e-5; % Tolerance for the difference between steps in Steepest Descent (|x_n - x_n+1| < tolX)
tolSD = 1e-5; % Tolerance for the f-gradient norm in Steepest Descent (|grad f| < tolSD)
maxIterSD = 50; % Maximum number of iterations for Steepest Descent
delta = 0.0001; % delta for Line Search in Steepest Descent
maxIterLS = 50; % Maximum number of iterations for Line Search in Steepest Descent

deltaTunnel = 1e-3; % Step to check for a sign change in the tunneling function
maxIterTunnelLS = 500; % Maximum number of attempts per tunneling attempt
randomSearchAttempts = 200; % Recovery attempts when tunneling fails to find a level
searchBounds = [-8, 8]; % Domain used for random restarts

while iterTotal < maxTotalIter
    
    % Find the local minimum

    [x_min, f_min, tra] = SteepestDescent(fn, X0, tolSD, maxIterSD, delta, maxIterLS, tolX);
    
    TRA{end+1} = tra;

    LOG = [LOG; iterTotal + 1, x_min(1), x_min(2), f_min];

    % Tunneling Stage. Find a point on a contour line

    NC = true;
    iterTunnel = 0;
    maxTunnelIter = 60;
    x_niv = X0;
    f_niv = NaN;

    while NC && (iterTunnel < maxTunnelIter)

        eps_1 = 2*(tolSD)^(1/5)*(1 + norm(x_min)); % Epsilon for random updating
        dir = rand(1, 2)*2 - 1;
        dir = dir/norm(dir);
        XT = x_min + eps_1*dir;

        [x_niv, f_niv, tra2, NC] = TunnelRootSearch(fn, XT, x_min, maxIterTunnelLS, deltaTunnel);

        iterTunnel = iterTunnel + 1;

    end

    % Update point for next global iteration

    if NC
        best_x = x_min;
        best_f = f_min;

        for intento = 1:randomSearchAttempts
            X_rand = searchBounds(1) + (searchBounds(2) - searchBounds(1))*rand(1, 2);
            [x_rand_min, f_rand_min, tra_rand] = SteepestDescent(fn, X_rand, tolSD, maxIterSD, delta, maxIterLS, tolX);

            if f_rand_min < best_f
                best_x = x_rand_min;
                best_f = f_rand_min;
                tra2 = tra_rand;
            end
        end

        if best_f < f_min
            x_niv = best_x;
            f_niv = best_f;
            NC = false;
            TRA{end+1} = tra2;
            disp('Tunneling did not find a level; a random reset with a lower value was used.')
        else
            warning('No valid tunneling point was found, and no restart with a lower value was found.')
            break
        end
    end

    X0 = x_niv;
    iterTotal = iterTotal + 1;

end


%% GRAPHS

fprintf('Minimum point found: (%.10f, %.10f)\n', x_min(1), x_min(2))
fprintf('Minimum value found: %.10e\n', f_min)

close all

XX = linspace(-4, 4, 350);
YY = linspace(-4, 4, 350);
ZZ1 = zeros(350, 350);

ptoM = {X0(1), X0(2)};

for i = 1:350
    for j = 1:350
        pto = {XX(i), YY(j)};
        ZZ1(j, i) = fn(pto);
    end
end

% Graph of the objective function, search paths, 
% and the start and end points of the simulation

figure(1)
clf
hold on

% Surface of the objective function
surf(XX, YY, ZZ1)
shading("interp")
contour3(XX, YY, ZZ1, 15, 'LineColor', [0.1 0.1 0.1], 'LineWidth', 0.1)

% Trayectorias de búsqueda
h1 = gobjects(length(TRA), 1);
for i = 1:length(TRA)
    h1(i) = plot3(TRA{i}(:, 1), TRA{i}(:, 2), TRA{i}(:, 3), ...
        'LineWidth', 1.2, 'Color', '#ff5900');
end

% Starting point (pink)
h2 = scatter3(X0_i(1), X0_i(2), fn(X0_C), 60, [1, 0.1059, 0.5922], "filled");

% End point (red)
h3 = scatter3(x_min(1), x_min(2), f_min, 60, [1 0 0.0745], "filled");

% Legend
legend([h1(1), h2, h3], ...
    {'Search paths', 'Initial point', 'End point'}, ...
    'Location', 'best', 'FontSize', 10)

xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$f(x,y)$', 'Interpreter', 'latex')

% Graph showing the function value versus the number of iterations. 
% This allows us to observe how the global optimization 
% method behaves as the number of iterations increases.

figure(2)
clf
hold on

plot(LOG(:, 1), log(max(LOG(:, 4), realmin)), "LineWidth", 1.5, "LineStyle", ":", "Color", '#bd00ff')

xlabel('Iteración')
ylabel('$\ln{(\max(f,\mathrm{realmin}))}$', 'Interpreter', 'latex')
grid on



%% OBJECTIVE FUNCTIONS

function z = three_hump_n(XX)

    x = XX{1};
    y = XX{2};

    z = 2*x^2 - 1.05*x^4 + x^6/6 + x*y + y^2;

end

function z = levy_n13(XX)

    x1 = XX{1};
    x2 = XX{2};

    term1 = sin(3*pi*x1)^2;
    term2 = (x1 - 1)^2 * (1 + sin(3*pi*x2)^2);
    term3 = (x2 - 1)^2 * (1 + sin(2*pi*x2)^2); 

    z = term1 + term2 + term3;

end

function z = esfera(XX)

    x = XX{1};
    y = XX{2};

    z = x^2 + y^2;

end

function z = ackley_n(XX)

x = XX{1};
y = XX{2};

z1 = 20*exp(-0.2*sqrt(0.5*(x^2 + y^2)));
z2 = exp(0.5*(cos(2*pi*x)+cos(2*pi*y)));

z = -z1 - z2 + exp(1) + 20;

end

function z = drop_wave(XX)

x = XX{1};
y = XX{2};

z1 = 1 + cos(12*sqrt(x^2 + y^2));
z2 = 0.5*(x^2 + y^2) + 2;

z = -z1/z2;

end