%%


R = 10;        %  MΩ
C = 1;         % C in nF
Vthr = 5;      %  mV
Vspk = 70;     %  mV
dt = 0.1;        %  ms

%% PART 1


time = 0:dt:100;   % Time 


V = zeros(size(time));  
I = zeros(size(time)); 
I(time >= 10 & time < 60) = 1;  % 1 nA
t = 1;
while t < length(time)-1 
    
    dV = dt * (-V(t)/(R*C) + I(t)/C);
    V(t+1) = V(t) + dV;

    if V(t+1) >= Vthr  
        V(t+1) = Vspk;  % Fire i+1
        V(t+2) = 0;     % Reset membrane at i+2
        t = t + 2;      
        continue;
     else
        t = t + 1;  
    end

end

figure;
yyaxis left
plot(time, V, 'r', 'LineWidth', 1.5);
ylabel('Membrane Voltage (mV)');
ylim([-10 75]); 

yyaxis right
plot(time, I, 'b', 'LineWidth', 1.5);
ylabel('Injected Current (nA)');
ylim([0 20]);

xlabel('Time (ms)');
title('IAF Model Response to Step Current Injection');
grid on;
legend('Membrane Voltage', 'Injected Current');

%%

Tmax = 1000;   
Ipeak = 1;      

% Frequencies (Hz)
freqs = [1, 2, 5, 10, 20, 50, 100];


time = 0:dt:Tmax;  
num_steps = length(time);
spike_counts = zeros(size(freqs));

%% 
for f_idx = 1:length(freqs)
    freq = freqs(f_idx);
    

    I_sin = Ipeak * sin(2 * pi * freq * time / 1000);  

    
    V_sin = zeros(num_steps,1);  
    spike_count = 0;

    % Euler integration loop using YOUR CODE structure
    t = 1;
    while t < length(time)-1  
        dV = dt * (-V_sin(t)/(R*C) + I_sin(t)/C);
        V_sin(t+1) = V_sin(t) + dV;

      
        if V_sin(t+1) >= Vthr  
            V_sin(t+1) = Vspk;  % Fire spike
            V_sin(t+2) = 0;     
            spike_count = spike_count + 1;  
            t = t + 2;           
            continue;
        else
            t = t + 1;  
        end
    end

  
    spike_counts(f_idx) = spike_count;

    % Plot 
    if freq == 10
        figure;
        yyaxis left
        plot(time, V_sin, 'r', 'LineWidth', 1.5);
        ylabel('Membrane Voltage (mV)');
        ylim([-10 75]);

        yyaxis right
        plot(time, I_sin, 'b', 'LineWidth', 1.5);
        ylabel('Injected Current (nA)');
        ylim([-1.5 1.5]);

        xlabel('Time (ms)');
        title(['IAF Response to Sinusoidal Stimulation (', num2str(freq), ' Hz)']);
        grid on;
        legend('Membrane Voltage', 'Injected Current');
    end
end



%% Spike Count vs. Frequency
figure;
plot(freqs, spike_counts, 'ko-', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
xlabel('Stimulus Frequency (Hz)');
ylabel('Spike Count (1s interval)');
title('Spike Count vs. Stimulus Frequency');
grid on;
set(gca, 'XScale', 'log');  % Log scale for frequency



