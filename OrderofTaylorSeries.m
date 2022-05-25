% Finding Order of Taylor Series

syms x; warning('off');
f = x - 1 - sin(x)/2;   %Function
aa = subs(f,x,pi/2);    %repacing x by pi/2
res2 = double(aa);      %Calculating true value
for i = 1:inf
    Tay = taylor(f, x, 'Order', i);     %For ith order approximation
    a = subs(Tay, x, pi/2);     %Substitution
    res1 = double(a);
    res = abs(res2 - res1);     %Absolute Error
    disp(res);
    disp(i);    
    if (res<=0.015)
        break;      %Required error achieved
    end
end