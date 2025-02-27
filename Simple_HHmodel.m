% Hodgkin-Huxley Model

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
h = alpha_h0/(alpha_h0 + beta_h0);% 0.5961


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


%% (i) 

figure;
plot(time, m_trace, 'r', time, n_trace, 'g', time, h_trace, 'b', 'LineWidth', 1.5);
title('Evolution of the HH variables m,n,h following an action potential');

xlabel('Time (ms)');
ylabel('HH Variable Value');
ylim([-0.5 1.5])

yyaxis right
plot(time, Iinj, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([-10 80]);
legend('m', 'n', 'h','Inject Current');

grid off;

%%




%% (iii)

figure;
plot(time, I_Na_trace, time, I_K_trace, time, I_L_trace, time, I_C_trace, 'LineWidth', 1.5);
title('Relative evolution of capacitive, leak, sodium and potassium currents');

xlabel('Time (ms)');
ylabel('Current (uA/cm^2)');
ylim([-1200 1200]);

yyaxis right
plot(time, Iinj, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([-5 100]);
legend('I_{Na}', 'I_{K}', 'I_{L}', 'I_{C}','I_{Inject}');
grid off;

%%




%% (ii)
figure;
plot(time, G_Na_trace, time, G_K_trace, 'LineWidth', 1.5);
title('Relative evolution of sodium and potassium conductances');

xlabel('Time (ms)');
ylabel('Conductance (mS/cm^2)');
ylim([-5 35]);

yyaxis right
plot(time, Iinj, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([-5 100]);
legend('G_{Na}', 'G_{K}','I_{Inject}');

grid off;

%%


%% 
figure;
plot(time, V_trace, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Hodgkin-Huxley Model');
yline(-55, '--r', 'Threshold');
yline(-68, '--k', 'Resting Potential');

ylim([-100 50]);


yyaxis right
plot(time, Iinj, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
legend('Membrane Potential','I_{Inject}');

ylim([-5 100]);

grid off;



%% (iV)

% Threshold behavior
I_ext_values = [3.5,3.995, 3.9969,3.997]; % Different current amplitudes
time1 = 0:dt:30;
V_results = zeros(length(I_ext_values), length(time1)); % Store membrane potential for each I_ext

figure;
hold on;
for i = 1:length(I_ext_values)
    V = V_rest;
    n = 0.3177;
    m = 0.0529;
    h = 0.5961;

    

    for t = 1:length(time1)
        dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
        dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
        dh = alpha_h(V) * (1 - h) - beta_h(V) * h;

        n = n + dt * dn;
        m = m + dt * dm;
        h = h + dt * dh;


        Iion = INa(V, m, h) + IK(V, n) + IL(V);
        I_C = I_ext_values(i) - Iion;

        dV = I_C / C;
    
        V = V + dt * dV;
    
        V_results(i, t) = V;
    
    end
    plot(time1, V_results(i,:),  'DisplayName', ['I_{inj} = ', num2str(I_ext_values(i)), ' μA/cm^2'], 'LineWidth', 1.5);
end

title('Threshold Spiking Behavior');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
yline(-53.7, '--k', 'Threshold');
legend;
grid off;
hold off;

%%


%% (iV)

% Rebound Behavior

V = V_rest; % Reset initial condition
n = 0.3177;
m = 0.0529;
h = 0.5961;
time2 = 0:dt:50;
I_ext = zeros(size(time2));
V_traceR = zeros(size(time2));

I_ext(time2 >= 5 & time2 <= 10) = -20;

for i = 1:length(time2)
    V_traceR(i) = V;

    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;

    % Update gating variables
    n = n + dt * dn;
    m = m + dt * dm;
    h = h + dt * dh;

    % Calculate total current and update membrane potential
    Iion = INa(V, m, h) + IK(V, n) + IL(V);
    I_C = I_ext(i) - Iion;

    dV = I_C / C;
    V = V + dt * dV;
    
end

figure;
plot(time2, V_traceR, 'LineWidth', 1.5);
title('Rebound Spiking Behavior');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
yline(-65, '--k', 'Resting Potential');
ylim([-150 55]);

yyaxis right
plot(time2, I_ext, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([-50 600]);
legend('Membrane Potential', 'Resting Potential','I_{Inject}');
hold off;

grid off;






%% (Vi)
%Reduced HH Model
time4 = 0:dt:t_max;

%Original One

V = V_rest; % Reset initial condition
n = 0.3177;
m = 0.0529;
h = 0.5961;
Iinj_2 = zeros(size(time4));
Iinj(time4 >= 5 & time4 < 10) = 5;


V_trace_Original = zeros(size(time4));

for i = 1:length(time4)
    

    V_trace_Original(i) = V;

    % Calculate gating variable derivatives
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;

    % Update gating variables
    n = n + dt * dn;
    m = m + dt * dm;
    h = h + dt * dh;

    % Calculate total current and update membrane potential
    Iion = INa(V, m, h) + IK(V, n) + IL(V);
    I_C = Iinj(i) - Iion;

    dV = I_C / C;
    
    V = V + dt * dV;
    
end

% Reduced model


V = V_rest; % Reset initial condition
n = 0.3177;
m = 0.0529;

V_trace_reduced = zeros(size(time4));


for i = 1:length(time4) 
    V_trace_reduced(i) = V;

    
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;

    
    n = n + dt * dn;
    m = m + dt * dm;

   
    h_reduced = 0.89 - (1.1 * n);

    
    Iion = INa(V, m, h_reduced) + IK(V, n) + IL(V);
    I_C = Iinj(i) - Iion;  % Use Iinj (not Iinj_2)

    dV = I_C / C;
    V = V + dt * dV;
end
    




figure;
plot(time4, V_trace_reduced, time4, V_trace_Original, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
ylim([-100 60]);


yyaxis right
plot(time4, Iinj, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([-5 100]);
legend('Reduced HH model', 'Original HH model', 'I_{Inject}');

title('Reduced VS. Original Hodgkin-Huxley Model');
grid off;

%%



%% change alpha_n (V)
V = V_rest; % Reset initial condition
n = 0.3177;
m = 0.0529;
h = 0.5961;
time3 = 0:dt:40;

Iinj_new = zeros(size(time3));
Iinj_new(time3 >= 5 & time3 < 10) = 5;

V_trace_Original = zeros(size(time3));

for i = 1:length(time3)
    

    V_trace_Original(i) = V;

    % Calculate gating variable derivatives
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;

    % Update gating variables
    n = n + dt * dn;
    m = m + dt * dm;
    h = h + dt * dh;

    % Calculate total current and update membrane potential
    Iion = INa(V, m, h) + IK(V, n) + IL(V);
    I_C = Iinj_new(i) - Iion;

    dV = I_C / C;
    
    V = V + dt * dV;
    
end


alpha_n = @(V) 0.01 * (22.5 + (-65-V))^2 / 110;

V = V_rest; % Reset initial condition
n = 0.3177;
m = 0.0529;
h = 0.5961;
V_trace_Change = zeros(size(time3));

for i = 1:length(time3)
    V_trace_Change(i) = V;

    % Calculate gating variable derivatives
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;

    % Update gating variables
    n = n + dt * dn;
    m = m + dt * dm;
    h = h + dt * dh;

    % Calculate total current and update membrane potential
    Iion = INa(V, m, h) + IK(V, n) + IL(V);
    I_C = Iinj_new(i) - Iion;

    dV = I_C / C;
    V = V + dt * dV;
    
end

figure;
plot(time3, V_trace_Change, time3, V_trace_Original, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');

title('Hodgkin-Huxley Model with alpha change');
ylim([-100 60]);

yyaxis right
plot(time3, Iinj_new, 'k', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([-5 100]);

legend('HH model with Alpha change', 'Original HH model', 'I_{Inject}');
grid off;
%% 
