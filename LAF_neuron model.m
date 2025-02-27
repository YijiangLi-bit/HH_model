%%
% Initial Paramater
R = 10;        %  MΩ
C = 1;         % C in nF
Vthr = 5;      %  mV
Vspk = 70;     %  mV
dt = 0.1;        %  ms

%%


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

%%




