% Hodgkin-Huxley Model

% Initial Variables
C = 1;
E_L = -61; % mV
E_K = -77; 
E_Na = 55; 
g_L = 0.3; % mS/cm^2
g_K = 36; 
g_Na = 120; 
V_rest = -65; % mV
dt = 0.01; % ms
t_max = 25; % ms
time = 0:dt:t_max;

alpha_n0 = 0.01 * (10) / (exp(10 / 10) - 1);
beta_n0 = 0.125 * exp(0 / 80);
alpha_m0 = 0.1 * (25) / (exp((25) / 10) - 1);
beta_m0 = 4 * exp(0 / 18);
alpha_h0 = 0.07 * exp(0 / 20);
beta_h0 = 1 / (exp(30 / 10) + 1);


Iinj = zeros(size(time));
Iinj(time >= 5 & time < 10) = 5;
% Initial conditions
V = V_rest; % membrane potential
n = alpha_n0/(alpha_n0 + beta_n0); %0.3177
m = alpha_m0/(alpha_m0 + beta_m0); %0.0529
h = alpha_h0/(alpha_h0 + beta_h0); %0.5961


%For Graph
V_trace = zeros(size(time));
n_trace = zeros(size(time));
m_trace = zeros(size(time));
h_trace = zeros(size(time));
G_Na_trace = zeros(size(time));
G_K_trace = zeros(size(time));
I_Na_trace = zeros(size(time));
I_K_trace = zeros(size(time));
I_L_trace = zeros(size(time));
I_C_trace = zeros(size(time));



%Function
alpha_n = @(V) 0.01 * (10 + (-65-V)) / (exp((10 + (-65-V)) / 10) - 1);
beta_n = @(V) 0.125 * exp((-65-V) / 80);
alpha_m = @(V) 0.1 * (25 + (-65 - V)) / (exp((25 + (-65-V)) / 10) - 1);
beta_m = @(V) 4 * exp((-65-V) / 18);
alpha_h = @(V) 0.07 * exp((-65-V) / 20);
beta_h = @(V) 1 / (exp((30 + (-65-V)) / 10) + 1);

INa = @(V, m, h) g_Na * m^3 * h * (V - E_Na);
IK = @(V, n) g_K * n^4 * (V - E_K);
IL = @(V) g_L * (V - E_L);

% Euler Method

for i = 1:length(time)
    

    V_trace(i) = V;
    n_trace(i) = n;
    m_trace(i) = m;
    h_trace(i) = h;
    G_Na_trace(i) = g_Na* m^3 * h;
    G_K_trace(i) = g_K* n^4;
    I_Na_trace(i) = INa(V, m, h);
    I_K_trace(i) = IK(V, n);
    I_L_trace(i) =  IL(V);


    % Gating Variable 
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
    n = n + dt * dn;
    m = m + dt * dm;
    h = h + dt * dh;

    % HH equation

    Iion = INa(V, m, h) + IK(V, n) + IL(V);
    I_C = Iinj(i) - Iion;

    dV = I_C / C;
    
    V = V + dt * dV;
    
    I_C_trace(i) = I_C;
end
%%
