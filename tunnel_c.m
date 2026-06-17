function z = tunnel_c(fn, XX, XM)

n=1;

num = fn(XX)-fn(XM);

den = sqrt( ( XX{1} - XM{1} )^2 +  ( XX{2} - XM{2} )^2);

z = num/den^n;

end