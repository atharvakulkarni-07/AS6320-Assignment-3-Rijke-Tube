% Simulation parameters
numSteps= 1000;% number of timesteps
% total duration of experiment
totalTime = 40;% total time of experiment
timeStep= totalTime /numSteps; % dt
flamePos= 0.3;  % m (original xf)

% Physical constants
gammaVal = 1.4; % g
speedOfSound = 399.6;  % c0
baseVelocity= 0.5; % u0
MachNum = baseVelocity / speedOfSound; % M
numModes = 3; % number of Galerkin modes (J)

% Parameter sweeps
gainList = 0.2:0.01:1.4; % non-dimensional heater power (ks)
numGains = length(gainList); % n1
delayList = 0.2:0.01:0.8; % time-lag (tau)
numDelays = length(delayList);% n3

% Preallocate output arrays
rmsForward= zeros(numDelays,numGains); % urms1
dynStateA= zeros(numDelays, numModes, numSteps); %y1
dynStateB= zeros(numDelays,numModes, numSteps);  % y2

% Forward integration loop
for delayIdx = 1:numDelays
    % Initialize state for this delay
    dynStateA(delayIdx, :,: ) = 0; %y1
    dynStateB(delayIdx,:, :) = 0;  %y2
    dynStateA(delayIdx, 1, 1) = 0.18;% initial mode-1 excitation

    velHistF = zeros(numGains, numSteps); %u1 history

    for gainIdx = 1:numGains
        % Compute initial velocity perturbation
        velHistF(gainIdx,1) = 0;
        for modeIdx = 1:numModes
            velHistF(gainIdx,1) = velHistF(gainIdx,1) + dynStateA(delayIdx,modeIdx,1) * ...
                                   cos(modeIdx *pi *flamePos);
        end

        lagSteps = round(delayList(delayIdx) /timeStep); % n2

        % Phase 1: t < tau (no heat feedback)
        for t = 1:lagSteps
            for modeIdx = 1:numModes
                [k1a,l1a] = dynNoHeat(dynStateA(delayIdx,modeIdx,t),   dynStateB(delayIdx,modeIdx,t),   modeIdx);
                [k2a, l2a] = dynNoHeat(dynStateA(delayIdx,modeIdx,t)+0.5*timeStep*k1a, dynStateB(delayIdx,modeIdx,t)+0.5*timeStep*l1a, modeIdx);
                [k3a, l3a] = dynNoHeat(dynStateA(delayIdx,modeIdx,t)+0.5*timeStep*k2a, dynStateB(delayIdx,modeIdx,t)+0.5*timeStep*l2a, modeIdx);
                [k4a ,l4a] = dynNoHeat(dynStateA(delayIdx,modeIdx,t) + timeStep*k3a,       dynStateB(delayIdx,modeIdx,t)+timeStep*l3a, modeIdx);

                dynStateA(delayIdx,modeIdx,t+1) = dynStateA(delayIdx,modeIdx,t) + (timeStep/6)*(k1a + 2*k2a + 2*k3a + k4a);
                dynStateB(delayIdx,modeIdx,t+1) = dynStateB(delayIdx,modeIdx,t) + (timeStep/6)*(l1a + 2*l2a + 2*l3a + l4a);

                velHistF(gainIdx,t+1) = velHistF(gainIdx,t+1) + dynStateA(delayIdx,modeIdx,t+1) * cos(modeIdx*pi*flamePos);
            end
        end

        % Phase 2: t > tau (with heat feedback)
        for t = lagSteps+1:numSteps-1
            for modeIdx = 1:numModes
                [k1b, l1b] = dynWithHeat(dynStateA(delayIdx,modeIdx,t),   dynStateB(delayIdx,modeIdx,t),   modeIdx, gainList(gainIdx), velHistF(gainIdx,t-lagSteps), flamePos);
                [k2b, l2b] = dynWithHeat(dynStateA(delayIdx,modeIdx,t)+0.5*timeStep*k1b, dynStateB(delayIdx,modeIdx,t)+0.5*timeStep*l1b, modeIdx, gainList(gainIdx), velHistF(gainIdx,t-lagSteps), flamePos);
                [k3b, l3b] = dynWithHeat(dynStateA(delayIdx,modeIdx,t)+0.5*timeStep*k2b, dynStateB(delayIdx,modeIdx,t)+0.5*timeStep*l2b, modeIdx, gainList(gainIdx), velHistF(gainIdx,t-lagSteps), flamePos);
                [k4b, l4b] = dynWithHeat(dynStateA(delayIdx,modeIdx,t)+timeStep*k3b,       dynStateB(delayIdx,modeIdx,t)+timeStep*l3b,       modeIdx, gainList(gainIdx), velHistF(gainIdx,t-lagSteps), flamePos);

                dynStateA(delayIdx,modeIdx,t+1) = dynStateA(delayIdx,modeIdx,t) + (timeStep/6)*(k1b + 2*k2b + 2*k3b + k4b);
                dynStateB(delayIdx,modeIdx,t+1) = dynStateB(delayIdx,modeIdx,t) + (timeStep/6)*(l1b + 2*l2b + 2*l3b + l4b);

                velHistF(gainIdx,t+1) = velHistF(gainIdx,t+1) + dynStateA(delayIdx,modeIdx,t+1) * cos(modeIdx*pi*flamePos);
            end
        end

        
        dynStateA(delayIdx,:,1) = dynStateA(delayIdx,:,numSteps);
        dynStateB(delayIdx,:,1) = dynStateB(delayIdx,:,numSteps);
        rmsForward(delayIdx,gainIdx) = rms(velHistF(gainIdx,:));
    end
end

%backward integration (this will just look like the mirror code of the forrward one)
rmsBackward = zeros(numDelays, numGains);
for delayIdx = 1:numDelays
    velHistB = zeros(numGains, numSteps);
    for gainIdx = 1:numGains
        revIdx = numGains - gainIdx + 1; % reversed index
        % Initialize reversed velocity
        velHistB(revIdx,1) = 0;
        for modeIdx = 1:numModes
            velHistB(revIdx,1) = velHistB(revIdx,1) + dynStateA(delayIdx,modeIdx,1) * cos(modeIdx*pi*flamePos);
        end
        lagSteps = round(delayList(delayIdx) / timeStep);

        % t < tau (no heat)
        for t = 1:lagSteps
            for modeIdx = 1:numModes
                [k1c, l1c] = dynNoHeat(dynStateA(delayIdx,modeIdx,t),   dynStateB(delayIdx,modeIdx,t),   modeIdx);
                [k2c, l2c] = dynNoHeat(dynStateA(delayIdx,modeIdx,t) + 0.5*timeStep*k1c, dynStateB(delayIdx,modeIdx,t) + 0.5*timeStep*l1c, modeIdx);
                [k3c, l3c] = dynNoHeat(dynStateA(delayIdx,modeIdx,t) + 0.5*timeStep*k2c, dynStateB(delayIdx,modeIdx,t) + 0.5*timeStep*l2c, modeIdx);% corrected line
                [k4c, l4c] = dynNoHeat(dynStateA(delayIdx,modeIdx,t) + timeStep*k3c,       dynStateB(delayIdx,modeIdx,t)+ timeStep*l3c, modeIdx);

                dynStateA(delayIdx,modeIdx,t+1) = dynStateA(delayIdx,modeIdx,t) + (timeStep/6)*(k1c + 2*k2c + 2*k3c + k4c);
                dynStateB(delayIdx,modeIdx,t+1) = dynStateB(delayIdx,modeIdx,t) + (timeStep/6)*(l1c + 2*l2c + 2*l3c + l4c);

                velHistB(revIdx,t+1) = velHistB(revIdx,t+1) + dynStateA(delayIdx,modeIdx,t+1) * cos(modeIdx*pi*flamePos);
            end
        end

        % t > tau (with heat)
        for t = lagSteps+1:numSteps-1
            for modeIdx = 1:numModes
                [k1d, l1d] = dynWithHeat(dynStateA(delayIdx,modeIdx,t),dynStateB(delayIdx,modeIdx,t),   modeIdx, gainList(revIdx), velHistB(revIdx,t-lagSteps), flamePos);
                [k2d, l2d] = dynWithHeat(dynStateA(delayIdx,modeIdx,t)+0.5*timeStep*k1d, dynStateB(delayIdx,modeIdx,t)+0.5*timeStep*l1d, modeIdx, gainList(revIdx),velHistB(revIdx,t-lagSteps), flamePos);
                [k3d, l3d] = dynWithHeat(dynStateA(delayIdx,modeIdx,t)+0.5*timeStep*k2d, dynStateB(delayIdx,modeIdx,t)+0.5*timeStep*l2d, modeIdx,gainList(revIdx), velHistB(revIdx,t-lagSteps), flamePos);% corrected line
                [k4d, l4d] = dynWithHeat(dynStateA(delayIdx,modeIdx,t)+timeStep*k3d,       dynStateB(delayIdx,modeIdx,t)+timeStep*l3d,       modeIdx, gainList(revIdx),velHistB(revIdx,t-lagSteps), flamePos);

                dynStateA(delayIdx,modeIdx,t+1) = dynStateA(delayIdx,modeIdx,t) + (timeStep/6)*(k1d + 2*k2d + 2*k3d + k4d);
                dynStateB(delayIdx,modeIdx,t+1) = dynStateB(delayIdx,modeIdx,t) + (timeStep/6)*(l1d + 2*l2d + 2*l3d + l4d);

                velHistB(revIdx,t+1) = velHistB(revIdx,t+1) + dynStateA(delayIdx,modeIdx,t+1) * cos(modeIdx*pi*flamePos);
            end
        end

        
dynStateA(delayIdx,:,1) = dynStateA(delayIdx,:,numSteps);
dynStateB(delayIdx,:,1) = dynStateB(delayIdx,:,numSteps);
rmsBackward(delayIdx,revIdx) = rms(velHistB(revIdx,:));
    end
end

% plots
timeVec = linspace(0, totalTime, numSteps);
for delayIdx = 1:numDelays
    if abs(delayList(delayIdx) - 0.2) < 1e-6
        figure(1);
        scatter(gainList, rmsForward(delayIdx,:), 20, 'b', 'filled'); hold on;
        scatter(gainList, rmsBackward(delayIdx,:), 20, 'r', 'filled');
        xlabel('K');
        ylabel('|u_1|');
        legend('Forward', 'Backward');
        title('Bifurcation Plot');
        annotation('textarrow', [0.7 0.6], [0.3 0.15], 'String', 'Hopf Bifurcation');
        annotation('textarrow', [0.28 0.3], [0.8 0.55], 'String', 'Fold Bifurcation');
        grid on;
    end
end


figure(2);
mesh(rmsForward); hold on;
mesh(rmsBackward);
waterfall(gainList, delayList, rmsForward); hold on;
waterfall(gainList, delayList, rmsBackward);
xlabel('K');
ylabel('tau');
zlabel('|u_1|');
legend('Forward','Backward');
title('3D Bifurcation Plot');
annotation('textarrow', [0.6 0.45], [0.3 0.15], 'String', 'Hopf Bifurcation');
annotation('textarrow', [0.2 0.25], [0.6 0.5], 'String', 'Fold Bifurcation');
colorbar;


% functions
function [fA,fB] = dynNoHeat(aState, bState, modeIdx)
    kVal = modeIdx*pi; omega = kVal; base = pi;
    c1 = 0.1; c2 = 0.06;
    damp = (1/(2*pi)) * (c1*(omega/base) + c2*sqrt(base/ omega));
    fA = bState;
    fB = -2*damp*omega*bState - kVal^2 * aState;
end

function [fA,fB] = dynWithHeat(aState, bState, modeIdx, Kgain, delayedVel, pos)
    kVal = modeIdx*pi; omega = kVal; base = pi;
    c1 = 0.1; c2 = 0.06;
    damp = (1/(2*pi)) * (c1*(omega/base) + c2*sqrt(base/omega));
    source = modeIdx*pi*Kgain * (abs(sqrt(1/3 + delayedVel)) - sqrt(1/3)) * sin(modeIdx* pi *pos);
    fA = bState;
    fB = -1 * 2*damp*omega*bState - kVal^2 * aState - source;
end
