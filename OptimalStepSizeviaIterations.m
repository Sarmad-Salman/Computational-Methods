% Finding the optimal Step Size for forward, backward and center
% finite differentiation through iterations via Absolute Error


syms x; h(1) = 1; warning('off');
F = x - 1 - sin(x)/2; f = @(x) x - 1 - sin(x)/2;
ff = diff(F); fff = diff(F,2); aa = subs(fff,x,pi/2); %Finding Derivatives
abs_ff = vpa(subs(ff,x,pi/2)); abs_fff = vpa(subs(fff,x,pi/2));


%First Derivative Approximations
x_i = pi/2; x_iplus1 = pi/2 + h(1); x_iminus1 = pi/2 - h(1);
f_xi = f(x_i); f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);


%Forward first iteration
forw1 = (f_xplusi - f_xi) / h(1);
rel_err_for1(1) = abs(abs_ff - forw1)*100;


%Backward first iteration
back1 = (f_xi - f_xminusi) / h(1);
rel_err_back1(1) = abs(abs_ff - back1)*100;


%Center first iteration
cent1 = (f_xplusi - f_xminusi) / (2*h(1));
rel_err_cent1(1) = abs(abs_ff - cent1)*100;


%Second derivative Approximations
x_iplus2 = x_iplus1 + h(1); x_iminus2 = x_iminus1 - h(1);
f_xplusi2 = f(x_iplus2); f_xminusi2 = f(x_iminus2);


%Forward first iteration
forw2 = (f_xplusi2 - 2*f_xplusi + f_xi) / h(1)^2;
rel_err_for2(1) = abs(abs_fff - forw2)*100;


%Backward first iteration
back2 = (f_xi - 2*f_xminusi + f_xminusi2) / h(1)^2;
rel_err_back2(1) = abs(abs_fff - back2)*100;


%Center first iteraton
cent2 = (f_xplusi - 2*f_xi + f_xminusi) / h(1)^2;
rel_err_cent2(1) = abs(abs_fff - cent2)*100;


%Iterations for 1st Order Forward Approximation
for i = 2:inf
    h(i) = h(i-1)/10;
    x_iplus1 = pi/2 + h(i); x_iminus1 = pi/2 - h(i);
    f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);
    forw1 = (f_xplusi - f_xi) / h(i);
    rel_err_for1(i) = abs(abs_ff - forw1)*100;
    % When Absolute error starts getting bigger
     if (rel_err_for1(i) >= rel_err_for1(i-1))
         h_for1 = h(i) * 10;
         disp(['1st Order Forward Approximation Step Size = ' num2str(h_for1)]);
         break;
     end
end
subplot(321);loglog(h,rel_err_for1); title('1st Order Forward Approximation');
xlabel('Log Step Size'); ylabel('Log Error');
grid on; clear h; h(1)=1;
%Iterations for 1st Order Backward Approximation
for i = 2:inf
    h(i) = h(i-1)/10;
    x_iplus1 = pi/2 + h(i); x_iminus1 = pi/2 - h(i);
    f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);
    back1 = (f_xi - f_xminusi) / h(i);
    rel_err_back1(i) = abs(abs_ff - back1)*100;
    % When Absolute error starts getting bigger
    if (rel_err_back1(i) >= rel_err_back1(i-1))
        h_back1 = h(i) * 10;
        disp(['1st Order Backward Approximation Step Size = ' num2str(h_back1)]);
        break;
    end
end
subplot(322);loglog(h,rel_err_back1); title('1st Order Backward Approximation');
xlabel('Log Step Size'); ylabel('Log Error');
grid on; clear h ; h(1)=1;


%Iterations for 1st Order Centered Approximation
for i = 2:inf
    h(i) = h(i-1)/10;
    x_iplus1 = pi/2 + h(i); x_iminus1 = pi/2 - h(i);
    f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);
    cent1 = (f_xplusi - f_xminusi) / (2*h(i));
    rel_err_cent1(i) = abs(abs_ff - cent1)*100;
    % When Absolute error starts getting bigger
    if (rel_err_cent1(i) >= rel_err_cent1(i-1))
        h_cent1 = h(i) * 10;
        disp(['1st Order Center Approximation Step Size = ' num2str(h_cent1)]);
        break;
    end
end
subplot(323);loglog(h,rel_err_cent1); title('1st Order Centered Approximation');
xlabel('Log Step Size'); ylabel('Log Error');
grid on; clear h; h(1)=1;


%Iterations for 2nd Order Forward Approximation
for i = 2:inf
    h(i) = h(i-1)/10;
    x_iplus1 = pi/2 + h(i); x_iminus1 = pi/2 - h(i);
    f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);
    x_iplus2 = x_iplus1 + h(i); x_iminus2 = x_iminus1 - h(i);
    f_xplusi2 = f(x_iplus2); f_xminusi2 = f(x_iminus2);
    forw2 = (f_xplusi2 - 2*f_xplusi + f_xi) / h(i)^2;
    rel_err_for2(i) = abs(abs_fff - forw2)*100;
    % When Absolute error starts getting bigger
    if (rel_err_for2(i) >= rel_err_for2(i-1))
        h_for2 = h(i) * 10;
        disp(['2nd Order Forward Approximation Step Size = ' num2str(h_for2)]);
        break;
    end
end
subplot(324);loglog(h,rel_err_for2); title('2nd Order Forward Approximation');
xlabel('Log Step Size'); ylabel('Log Error');
grid on; clear h; h(1)=1;


%Iterations for 2nd Order Backward Approximation
for i = 2:inf
    h(i) = h(i-1)/10;
    x_iplus1 = pi/2 + h(i); x_iminus1 = pi/2 - h(i);
    f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);
    x_iplus2 = x_iplus1 + h(i); x_iminus2 = x_iminus1 - h(i);
    f_xplusi2 = f(x_iplus2); f_xminusi2 = f(x_iminus2);
    back2 = (f_xi - 2*f_xminusi + f_xminusi2) / h(i)^2;
    rel_err_back2(i) = abs(abs_fff - back2)*100;
    % When Absolute error starts getting bigger
    if (rel_err_back2(i) >= rel_err_back2(i-1))
        h_back2 = h(i) * 10;
        disp(['2nd Order Backward Approximation Step Size = ' num2str(h_back2)]);
        break;
    end
end
subplot(325);loglog(h,rel_err_back2); title('2nd Order Backward Approximation');
xlabel('Log Step Size'); ylabel('Log Error');
grid on; clear h; h(1)=1;


%Iterations for 2nd Order Centered Approximation
for i = 2:inf
    h(i) = h(i-1)/10;
    x_iplus1 = pi/2 + h(i); x_iminus1 = pi/2 - h(i);
    f_xplusi = f(x_iplus1); f_xminusi = f(x_iminus1);
    x_iplus2 = x_iplus1 + h(i); x_iminus2 = x_iminus1 - h(i);
    f_xplusi2 = f(x_iplus2); f_xminusi2 = f(x_iminus2);
    cent2 = (f_xplusi - 2*f_xi + f_xminusi) / h(i)^2;
    rel_err_cent2(i) = abs(abs_fff - cent2)*100;
    % When Absolute error starts getting bigger
    if (rel_err_cent2(i) >= rel_err_cent2(i-1))
        error=rel_err_cent2(i);
        h_cent2 = h(i) * 10;
        disp(['2nd Order Center Approximation Step Size = ' num2str(h_cent2)]);
        break;
    end
end
subplot(326);loglog(h,rel_err_cent2); title('2nd Order Centered Approximation');
xlabel('Log Step Size'); ylabel('Log Error'); grid on; clear h; h(1)=1;