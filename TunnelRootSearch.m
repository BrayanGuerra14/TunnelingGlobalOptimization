function [x_raiz, f_raiz, pto, NC] = TunnelRootSearch(fn, x_ini, x_min, maxIterLS, delta)

    if nargin < 5
        delta = 0.00001;
    end

    pto = [];
    NC = true;

    dim = length(x_ini);
    x = x_ini;
    xm = cell(1, dim);

    for i = 1:dim
        xm{i} = x_min(i);
    end

    [G, ~] = grad_T(fn, x, x_min);
    normaG = norm(G);
    tolGradTunnel = 1e-12;

    if normaG < tolGradTunnel
        x_raiz = x;
        xt0 = cell(dim, 1);

        for i = 1:dim
            xt0{i} = x(i);
        end

        ft0 = tunnel_c(fn, xt0, xm);
        f_raiz = ft0.f;
        return
    end

    dirBase = G/normaG;
    signosDireccion = [-1, 1];
    maxDelta = max(50*delta, 0.05*(1 + norm(x_min)));
    factorCrecimiento = 1.25;

    x_raiz = x;
    f_raiz = NaN;

    for k = 1:length(signosDireccion)

        signoDireccion = signosDireccion(k);
        alfaValor = 0;
        paso = delta;
        alfa = DualE(alfaValor, 1, 0);
        xt = cell(dim, 1);

        for i = 1:dim
            xt{i} = x(i) + signoDireccion*alfa*dirBase(i);
        end

        ft = tunnel_c(fn, xt, xm);
        signoAnterior = sign(ft.f);
        iterN = 0;
        cambio = false;

        while ~cambio && (iterN < maxIterLS)

            alfaValor = alfaValor + paso;
            alfa = DualE(alfaValor, 1, 0);
            xt = cell(dim, 1);

            for i = 1:dim
                xt{i} = x(i) + signoDireccion*alfa*dirBase(i);
            end

            ft = tunnel_c(fn, xt, xm);
            signoActual = sign(ft.f);

            if signoAnterior ~= signoActual
                cambio = true;
            else
                signoAnterior = signoActual;
                paso = min(factorCrecimiento*paso, maxDelta);
            end

            iterN = iterN + 1;
            pto = [pto; xt{1}.f xt{2}.f ft.f];

        end

        if cambio

            for i = 1:dim
                x_raiz(i) = xt{i}.f;
            end

            f_raiz = ft.f;
            NC = false;
            return

        end

    end

end
