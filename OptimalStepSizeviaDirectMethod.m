% Finding the optimal Step Size for forward, backward and center
% finite differentiation through direct method via Absolute Error


syms x;
F = x - 1 - sin(x)/2; f = @(x) x - 1 - sin(x)/2;
ff = diff(F); fff = diff(F,2);
E = pi/2*10^-16; aa = subs(fff,x,pi/2);
M = double(aa); h = nthroot((3*E/M),3);
abs_ff = vpa(subs(ff,x,pi/2)); abs_fff = vpa(subs(fff,x,pi/2));
x_i = pi/2; x_iplus1 = pi/2 + h; x_iminus1 = pi/2 - h;
f_xi = f(x_i); f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);


%forward_first
forw1 = (f_xplusi - f_xi) / h;
rel_err_for1 = abs(abs_ff - forw1)*100;


%backward_first
back1 = (f_xi - f_xminusi) / h;
rel_err_back1 = abs(abs_ff - back1)*100;


%center_first
cent1 = (f_xplusi - f_xminusi) / (2*h);
rel_err_cent1 = abs(abs_ff - cent1)*100;


%second derivative
x_iplus2 = x_iplus1 + h; x_iminus2 = x_iminus1 - h;
f_xplusi2 = f(x_iplus2); f_xminusi2 = f(x_iminus2);


%forward second
forw2 = (f_xplusi2 - 2*f_xplusi + f_xi) / h^2;
rel_err_for2 = abs(abs_fff - forw2)*100;


%backward second
back2 = (f_xi - 2*f_xminusi + f_xminusi2) / h^2;
rel_err_back2 = abs(abs_fff - back2)*100;


%center second
cent2 = (f_xplusi - 2*f_xi + f_xminusi) / h^2;
rel_err_cent2 = abs(abs_fff - cent2)*100;