%Finding value of H by Graphical Method

format short g; syms H;
Q = 5; S = 0.0002; B = 20; n = 0.03; %Initial Conditions
Eq1 = ((sqrt(S).*((B.*H).^(5/3))) ./ (n.*(B+2*H).^(2/3))) - Q == 0;
H_val = double(solve(Eq1, H)); disp(H_val); %Finding value of H
H = 0:0.01:10;
Eq1 = ((sqrt(S).*((B.*H).^(5/3))) ./ (n.*(B+2*H).^(2/3))) - Q; figure(1);
plot (H, Eq1); grid on; title('Graphical Approach'); %Plot of H and f[H]
xlabel('H = depth[m]'); ylabel('f(H)');  figure(2);
%______________Bisection Method_________________________
iter = 10; xro = 0; tr = H_val(1); t = 1:iter;
xu = input('Enter value of Xu = '); tempxu = xu; %Upper Bracket
xl = input('Enter value of Xl = '); tempxl = xl; %Lower Bracket
SD = input('Enter number of Significant Digits = '); tempSD = SD;
es = (0.5*10^(2-SD)); %For error upto 'n' significant figures
f = @(H)(((sqrt(S).*((B.*H).^(5/3))) ./ (n.*(B+2*H).^(2/3))) - Q);
fprintf('      \t\t\t\t\t\t Bisection Method\n');
fprintf('\t\tIterations \t\tXl \t\t\tXu \t\t\tXr \t\t\tEa \t\t\tEt\n');
table=[];   %Tabular Form
if(f(xu)*f(xl)<0)   %Checking if Root Exists
    for i = 1:iter
        xr = (xl + xu)/2;  %Bisection Method
        if(xr~=0)
            ea = abs((xr - xro)/xr) * 100;
            Rel_err(i) = ea;  %Array of Relative Errors
            et = abs ((tr - xr)/tr) * 100;
            Abs_err(i) = et;  %Array of True Errors
        end
        table=[table; i xl xu xr ea et];   %Tabular Data
        if(f(xl)*f(xr)<0)
            xu = xr; %Root lies in the lower subinterval
        elseif(f(xl)*f(xr)>0)
            xl = xr; %Root lies in the upper subinterval
        elseif(f(xl)*f(xr)==0)
            break;   %Root itself
        end
        if(ea<es)
            break;   %Error below specified
        end
        xro = xr;
    end
    disp(double(table));  %Display Table
    plot(t(1:length(Rel_err)), Rel_err); hold on; xlabel('Iterations');
    plot(t(1:length(Abs_err)), Abs_err); grid on; ylabel('True Percent Relative Error');
    figure(3);plot(t(1:length(Rel_err)), Rel_err); hold on; xlabel('Iterations');
    plot(t(1:length(Abs_err)), Abs_err); grid on; ylabel('True Percent Relative Error');
    title('Bisection Method'); legend('Bisection Relative Error','Bisection True Error');
    figure(2); else
    fprintf('No Root in this Vicinity');
end
%______________False Position___________________
xro = 0; clear Rel_err Abs_err; xl = tempxl; xu = tempxu;
fprintf('       \t\t\t\t\t\t False Position Method\n');
fprintf('\t\tIterations \t\tXl \t\t\tXu \t\t\tXr \t\t\tEa \t\t\tEt\n');
table=[];  %Tabular Form
if(f(xu)*f(xl)<0)  %Check for Root
    for i = 1:iter
        xr = xu - ((f(xu)*(xl-xu))/(f(xl)-f(xu))); %False Position
        if(xr~=0)
            ea = abs((xr - xro)/xr) * 100;
            Rel_err(i) = ea; %Array of Relative Errors
            et = abs ((tr - xr)/tr) * 100;
            Abs_err(i) = et; %Array of True Errors
        end
        table=[table; i xl xu xr ea et];    %Tabular Data
        if(f(xl)*f(xr)<0)
            xu = xr; %Root lies in the lower subinterval
        elseif(f(xl)*f(xr)>0)
            xl = xr; %Root lies in the upper subinterval
        elseif(f(xl)*f(xr)==0)
            break; %Root itself
        end
        if(ea<es)
            break;  %Error below specified
        end
        xro = xr;
    end
    disp(double(table));  %Display Table
    plot(t(1:length(Rel_err)), Rel_err); title('Bracketing Methods');
    plot(t(1:length(Abs_err)), Abs_err); figure(4);
    plot(t(1:length(Rel_err)), Rel_err); xlabel('Iterations'); hold on;
    plot(t(1:length(Abs_err)), Abs_err); ylabel('True Percent Relative Error');
    grid on; title('False Position');legend('False Position Relative Error',...
        'False Position True Error'); figure(2);
else
    fprintf('No Root in this Vicinity');
end
%__________Modified False Position Method________________________
xro = 0; clear Rel_err Abs_err; xl = tempxl; xu = tempxu;
fprintf('      \t\t\t\t\t\t Modified False Position Method\n');
fprintf('\t\tIterations \t\tXl \t\t\tXu \t\t\tXr \t\t\tEa \t\t\tEt\n');
table=[]; %Tabular Form
iu = 0; il = 0; fu = f(xu); fl = f(xl); %Modification
if(fu*fl<0)  %Check for Root
    for i = 1:iter
        xr = xu - ((fu*(xl-xu))/(fl-fu)); %False Position
        if(xr~=0)
            ea = abs((xr - xro)/xr) * 100;
            Rel_err(i) = ea; %Array of Relative Errors
            et = abs ((tr - xr)/tr) * 100;
            Abs_err(i) = et;  %Array of True Errors
        end
        table=[table; i xl xu xr ea et];  %Tabular Data
        if(fl*f(xr)<0)
            xu = xr;  %Root lies in the lower subinterval
            fu = f(xu);
            iu = 0;
            il = il + 1;
            if(il>=2)
                fl = fl/2;  %Stuck check condition
            end
        elseif(fl*f(xr)>0)
            xl = xr;   %Root lies in the upper subinterval
            fl = f(xl);
            il = 0;
            iu = iu + 1;
            if(iu>=2)
                fu = fu/2;  %Stuck check condition
            end
        elseif(fl*f(xr) == 0)
            break;  %Root itself
        end
        if(ea<es)
            break;  %%Error below specified
        end
        xro = xr;
    end
    disp(double(table));  %Display Table
    plot(t(1:length(Rel_err)), Rel_err); plot(t(1:length(Abs_err)), Abs_err);
    legend('Bisection Relative Error','Bisection True Error',...
        'False Position Relative Error','False Position True Error',...
        'Modified False Positon Relative Error',...
        'Modified False Positon True Error'); figure(5);
    plot(t(1:length(Rel_err)), Rel_err); xlabel('Iterations'); hold on;
    plot(t(1:length(Abs_err)), Abs_err);ylabel('True Percent Relative Error');
    grid on; title('Modified False Position');legend...
        ('Modified False Positon Relative Error',...
        'Modified False Positon True Error'); figure(2);
else
    fprintf('No Root in this Vicinity');
end