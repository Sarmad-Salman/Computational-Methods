%Open Methods
%______Fixed point Iteration_______________

format short g; syms H;  %Initial Conditions
Q = 5; S = 0.0002; B = 20; n = 0.03; iter = 10;
g = @(H) ((Q*n)^(3/5)*(B+2*H)^(2/5)/(B*S^(3/10))); %g(x) from f(x)
xi = input('Enter value of Xi = '); %Initial Value
xr = xi; tr = 0.7023;
SD = input('Enter number of Significant Digits = ');
es = (0.5*10^(2-SD));   %N significant digits error
t = 1:iter; tempxi = xi;
fprintf('       \t\t\t\t\t\t Fixed Point Iteration\n');
fprintf('\t\tIterations \t\tXi \t\t\tXr \t\t\tEa \t\t\tEt\n');
table=[]; %Tabular Form
for i=1:iter
    xro = xr;   %Updating Values
    xr = g(xro);    %Fixed point Iteration
    if(xr~=0)
        ea = abs((xr-xro)/xr)*100;
        Rel_err(i) = ea;    %%Array of Relative Errors
        et = abs ((tr - xr)/tr) * 100;
        Abs_err(i) = et;    %Array of True Errors
    end
    table=[table; i xi xr ea et];   %Tabular Data
    if(ea<es)
        break;  %Sufficient error achieved
    end
end
disp(double(table));    %Display Table
figure(1); plot(t(1:length(Rel_err)),Rel_err); grid on;
hold on; plot(t(1:length(Abs_err)), Abs_err); figure(2);
plot(t(1:length(Rel_err)),Rel_err); grid on; hold on;
plot(t(1:length(Abs_err)), Abs_err); xlabel('Iterations');
ylabel('True Percent Relative Error'); title('Fixed Point Iteration');
legend('Fixed Point Iteration Relative','Fixed Point Iteration True');figure(1);
%______Newton Raphson Method______________________________
clear Rel_err Abs_err; xi = tempxi;
fprintf('       \t\t\t\t\t\t Newton Raphson Method\n');
fprintf('\t\tIterations \t\tXi \t\t\tXr \t\t\tEa \t\t\tEt\n');
table=[]; xr = xi;
f = @(H)(((sqrt(S).*((B.*H).^(5/3))) ./ (n.*(B+2*H).^(2/3))) - Q);
der_f = matlabFunction(diff(f(H),H));   %Derivative of function
for i=1:iter
    xro = xr;   %Updating Values
    xr = xro - (f(xro)/der_f(xro)); %Newton Raphson
    if(xr~=0)
        ea = abs((xr-xro)/xr)*100;
        Rel_err(i) = ea;    %Array of Relative Errors
        et = abs ((tr - xr)/tr) * 100;
        Abs_err(i) = et;    %Array of True Errors
    end
    table=[table; i xi xr ea et];   %Tabular Data
    if(ea<es)
        break;  %Sufficient error achieved
    end
end
disp(double(table));    %Display Table
plot(t(1:length(Rel_err)),Rel_err);xlabel('Iterations');
plot(t(1:length(Abs_err)),Abs_err);title('Open Methods');
ylabel('True Percent Relative Error'); figure (3);
plot(t(1:length(Rel_err)),Rel_err);xlabel('Iterations'); hold on;
plot(t(1:length(Abs_err)),Abs_err);ylabel('True Percent Relative Error')
title('Newton Raphson Method'); grid on; legend('Newton Raphson Relative',...
    'Newton Raphson True'); figure(1);
%_________Secant Method_______________________
clear Rel_err Abs_err; xi = tempxi;
xi_1 = input('Enter value of a previous point = '); %Previous step
fprintf('       \t\t\t\t\t\t Secant Method\n');
fprintf('\tIterations \t\tXi-1 \t\t\tXi \t\t\tXi+1 \t\tEa \t\t\tEt\n');
table=[];
f = @(H)(((sqrt(S).*((B.*H).^(5/3))) ./ (n.*(B+2*H).^(2/3))) - Q);
for i=1:iter
    xr = xi - ((f(xi)*(xi_1 - xi)) / (f(xi_1) - f(xi))); %Secant
    if(xr ~= 0)
        ea = abs ((xr - xi) / xr) * 100;
        Rel_err(i) = ea;    %Array of Relative Errors
        et = abs ((tr - xr) / tr) * 100;
        Abs_err(i) = et;    %Array of True Errors
    end
    table=[table; i xi_1 xi xr ea et];  %Tabular Data
    if(ea<es)
        break;  %Sufficient error achieved
    end
    xi_1 = xi;      %Updating Values
    xi = xr;
end
disp(double(table));    %Display Table
plot(t(1:length(Rel_err)),Rel_err); plot(t(1:length(Abs_err)),Abs_err);
legend('Fixed Point Iteration Relative','Fixed Point Iteration True',...
    'Newton Raphson Relative','Newton Raphson True',...
    'Secant Relative','Secant True'); figure(4);
plot(t(1:length(Rel_err)),Rel_err);hold on; plot(t(1:length(Abs_err)),Abs_err);
xlabel('Iterations'); ylabel('True Percent Relative Error');
title('Secant Method'); grid on; legend('Secant Relative','Secant True');
figure(1);