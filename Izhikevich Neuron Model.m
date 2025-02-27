% Initial Variables
a = 0.02;
b = 0.2;
c = -65;
d = 8;

dt = 1;
T = 100;
t = 0:dt:T;

%%

v = -65;
u = b*v;

I = zeros(size(t));
I(t >= 10 & t < 60) = 1;

V = zeros(size(t));
U = zeros(size(t));


for i = 1:length(t)
    V(1) = v;
    U(1) = u;

    if v >= 30
        V(i) = 30;
        v = c;
        u = u + d;
    else
        V(i) = v;
    end

    dv = ((0.04*v + 5).*v + 140 - u + I(i))*dt;
    du = a*(b*v - u).*dt;
    
    v = v + 0.5 * dv;
    v = v + 0.5 * dv;
    u = u + du;
    
    U(i) = u;

end
%%
