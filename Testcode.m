a = 0.02;
b = 0.2;
c = -65;
d = 8;

dt = 1;
T = 100;
t = 0:dt:T;

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

figure;

% Plot membrane potential over time

plot(t, V, 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Izhikevich Neuron Model');

yyaxis right
plot(t, I, 'K', 'LineWidth', 1.5);
ylabel('Injected Current','Color','k');
ylim([-5 250]);

grid on;


a = 0.02;
b = 0.2;
c = -65;
d = 8;
dt = 1;
T = 1000;
t = 0:dt:T;
freq_list = [1 2 5 10 20 50 100];
spike_counts = zeros(size(freq_list));
for idx = 1:length(freq_list)
    f = freq_list(idx);
    v = -65;
    u = b*v;
    V = zeros(size(t));
    U = zeros(size(t));
    I_inj = 10 * sin(2*pi*f*(t/1000));
    spike_count = 0;
    V(1) = v;
    U(1) = u;
    for i = 2:length(t)
        dv = (0.04*v^2 + 5*v + 140 - u + I_inj(i))*dt;
        du = a*(b*v - u)*dt;
        v = v + 0.5 * dv;
        v = v+ 0.5 * dv;
        u = u + du;
        if v >= 30
            spike_count = spike_count + 1;
            V(i) = 30;
            v = c;
            u = u + d;
        else
            V(i) = v;
        end
        U(i) = u;
    end
    spike_counts(idx) = spike_count;
    if f == 10
        V_10 = V;
        I_10 = I_inj;
    end
end
figure;
subplot(2,1,1)
plot(t, I_10, 'b', 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Input Current (nA)')
title('Sinusoidal Input Current (10 Hz)')
grid on
subplot(2,1,2)
plot(t, V_10, 'r', 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
title('Spiking Pattern (10 Hz)')
grid on
figure;
plot(freq_list, spike_counts, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r','MarkerSize', 8)
xlabel('Stimulus Frequency (Hz)')
ylabel('Spike Count')
title('Spike Count vs. Stimulus Frequency')
set(gca, 'XScale', 'log');  % Log scale for frequency

grid on

