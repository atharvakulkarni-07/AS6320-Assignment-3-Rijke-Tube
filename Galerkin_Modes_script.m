% simulation Configuration
steps = 1000; %Total number of steps (N)
duration = 100; %Duration of experiment (T)
delta_t = duration / steps;  % (dt)

% physical Constants
gamma = 1.4; % (g)
c_sound = 399.6; % (c0)
velocity_mean = 0.5; % (u0)
Mach = velocity_mean / c_sound; % Calculating Mach number: M = u0/c0 % (M)
% gain = 0.05;
% gain = 0.6;
% gain = 1; 
gain = 0.3;



x_flame = 0.29; % (xf)
% lag_time = 0.2; % (tau)
% lag_time = 0.5;
lag_time = 0.45 * pi;

% initializations
modes = 3; % (J)
eta = zeros(modes, steps); % (y1)
eta_dot = zeros(modes, steps); % (y2)
vel_prime = zeros(1, steps); % (u)
pres_prime = zeros(1, steps); % (p)


% Initial conditions
eta(:,1) = 0; 
% eta(1,1) = 0.5
% eta(1,1) = 0.2
eta(1,1) = 0.18;     % Excitation for mode-1
eta_dot(:,1) = 0;
vel_prime(1) = eta(1,1) * cos(pi * x_flame); % Setting initial u'(1) using y1(1,1) and flame position xf
pres_prime(1) = 0; % Initial pressure perturbation p'(1) is set to zero
delay_idx = round(lag_time / delta_t); % Find how many steps correspond to the time lag tau

% pre combustion
for t = 1:delay_idx
    for m = 1:modes
        [dy1, dy2] = preCombustion(eta(m,t), eta_dot(m,t), m);
        [k2, l2] = preCombustion(eta(m,t) + 0.5*delta_t*dy1, eta_dot(m,t) + 0.5*delta_t*dy2, m);
        [k3, l3] = preCombustion(eta(m,t) + 0.5*delta_t*k2, eta_dot(m,t) + 0.5*delta_t*l2, m);
        [k4, l4] = preCombustion(eta(m,t) + delta_t*k3, eta_dot(m,t) + delta_t*l3, m);
        eta(m,t+1) = eta(m,t) + (delta_t/6)*(dy1 + 2*k2 + 2*k3 + k4);
        eta_dot(m,t+1) = eta_dot(m,t) + (delta_t/6)*(dy2 + 2*l2 + 2*l3 + l4);

        vel_prime(t+1) = vel_prime(t+1) + eta(m,t+1) * cos(m*pi*x_flame); % Summing modal contributions to velocity perturbation u'(n+1)
        pres_prime(t+1) = pres_prime(t+1) + eta_dot(m,t+1) * ((-gamma*Mach)/(m*pi)) * sin(m*pi*x_flame); % Summing modal contributions to pressure perturbation p'(n+1)
    end
end

% post combustion with feedback 

for t = delay_idx+1:steps-1
    for m = 1:modes
        [dy1, dy2] = postCombustion(eta(m,t), eta_dot(m,t), m, gain, vel_prime(t - delay_idx), x_flame);
        [k2, l2] = postCombustion(eta(m,t) + 0.5*delta_t*dy1, eta_dot(m,t) + 0.5*delta_t*dy2, m, gain, vel_prime(t - delay_idx), x_flame);
        [k3, l3] = postCombustion(eta(m,t) + 0.5*delta_t*k2, eta_dot(m,t) + 0.5*delta_t*l2, m, gain, vel_prime(t - delay_idx), x_flame);
        [k4, l4] = postCombustion(eta(m,t) + delta_t*k3, eta_dot(m,t) + delta_t*l3, m, gain, vel_prime(t - delay_idx), x_flame);

        eta(m,t+1) = eta(m,t) + (delta_t/6)*(dy1 + 2*k2 + 2*k3 + k4);
        eta_dot(m,t+1) = eta_dot(m,t) + (delta_t/6)*(dy2 + 2*l2 + 2*l3 + l4);

        vel_prime(t+1) = vel_prime(t+1) + eta(m,t+1) * cos(m*pi*x_flame); % Summing modal contributions to velocity perturbation u'(n+1)
        pres_prime(t+1) = pres_prime(t+1) + eta_dot(m,t+1) * ((-gamma*Mach)/(m*pi)) * sin(m*pi*x_flame); % Summing modal contributions to pressure perturbation p'(n+1)
    end
end

%energy Analysis 

energy = (0.5 * (pres_prime.^2) + 0.5 * (gamma * Mach * vel_prime).^2) / ((gamma * Mach)^2); % Calculating acoustic energy: e(n) = (p'(n)^2 + (g*M*u'(n))^2)/(g*M)^2 % (e) 
env_energy = envelope(energy, 200, 'peak'); % Extracting the peak envelope from energy to see growth/decay trends % (ee)

%let's look at the figuures...
time_vec = linspace(0, duration, steps);

figure(1);
plot(time_vec, vel_prime, 'b', 'LineWidth', 2);
title("u' Over Time"); xlabel('Time (s)'); ylabel("u'"); grid on;

figure(2);
plot(time_vec, eta(1,:), 'm', 'LineWidth', 2);
title('\eta_1(t)'); xlabel('Time (s)'); ylabel('\eta_1'); grid on;

figure(3);
plot(time_vec, eta(2,:), 'm', 'LineWidth', 2);
title('\eta_2(t)'); xlabel('Time (s)'); ylabel('\eta_2'); grid on;

figure(4);
plot(time_vec, eta(3,:), 'm', 'LineWidth', 2);
title('\eta_3(t)'); xlabel('Time (s)'); ylabel('\eta_3'); grid on;

figure(5);
plot(time_vec, env_energy, 'r', 'LineWidth', 2);
title('Amplitude of Energy Envelope'); xlabel('Time (s)');
ylabel('Normalized Energy'); legend('Nonlinear'); grid on;

figure(6);
tiledlayout(2,2);
nexttile; plot(time_vec, vel_prime, 'b', 'LineWidth', 2); title("u'(t)");
nexttile; plot(time_vec, eta(1,:), 'g', 'LineWidth', 2); title('\eta_1');
nexttile; plot(time_vec, eta(2,:), 'g', 'LineWidth', 2); title('\eta_2');
nexttile; plot(time_vec, eta(3,:), 'g', 'LineWidth', 2); title('\eta_3');
ylim([-0.4 0.4]);

%governing Equations
function [d1, d2] = preCombustion(a, b, idx)
    k = idx*pi;
    omega = k;
    base_freq = pi;
    c1 = 0.1; c2 = 0.06;
    damping = (1/(2*pi))*(c1*(omega/base_freq) + c2*sqrt(base_freq/omega));
    d1 = b; % dy1/dt = y2, modal displacement rate equals modal velocity
    d2 = -k^2 * a - 2*damping*omega*b; % dy2/dt = -k^2*y1 - 2*zeta*omega*y2, dynamics without heat feedback
end

function [d1, d2] = postCombustion(a, b, idx, gain, delayed_u, pos_x)
    k = idx*pi;
    omega = k;
    base_freq = pi;
    c1 = 0.1; c2 = 0.06;
    damping = (1/(2*pi))*(c1*(omega/base_freq) + c2*sqrt(base_freq/omega));
    source_term = (2*k*gain)*(sqrt(abs(1/3 + delayed_u)) - sqrt(1/3)) * sin(idx*pi*pos_x); % Nonlinear heat feedback source term added after delay
    d1 = b; % dy1/dt = y2, modal displacement rate equals modal velocity
    d2 = -k^2 * a - 2*damping*omega*b - source_term; % dy2/dt with additional feedback forcing from source_term
end
