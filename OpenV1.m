
% Initilaize Quantities
N = 101;
M = 1; % total mass
m = M / N; % mass per unit length (uniform)
T = 1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curve: straight line between endpoints

x = zeros(3,N); % position vector
v = zeros(3,N); % velocity vector
R = (N - 1) / 2;
Z = 2 * N; %depth

x(:,1) = [0,0,0]; % stationary endpoint
x(:,N) = [R,0,Z];

for i = 2: (N-1)
    x(:,i) = (i-1) / (N-1) * x(:,N);
end

w = 0.5; %angular velocity of endpoint


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arc Length at Rest
l = 0.5 * norm(x(:,1)-x(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate Motion

tmax = 200; 
dt = 0.05; 
clockmax = tmax / dt;


%Video Setup
writerObj1 = VideoWriter('Open_V1.avi');
open(writerObj1);

% Initial Frame
scatter3(x(1,:),x(2,:),x(3,:))
title('Insert Title')
axis([- 3 * N, 3 * N,- 3 * N, 3 * N, 0, Z])
axis manual
currentFrame = getframe(gcf);
writeVideo(writerObj1, currentFrame);

for clock = 1: clockmax
    t = clock * dt;

    % Update Interior Velocities 

    for i = 2: (N-1)
        D1 = x(:,i+1) - x(:,i);
        D2 = x(:,i-1) - x(:,i);
        E1 = max([norm(D1) - l, 0]); % excess arc length 1
        E2 = max([norm(D2) - l, 0]); % excess arc length 2
        F1 = T * E1 * D1 / norm(D1);
        F2 = T * E2 * D2 / norm(D2);
        A = (F1 + F2) / m;
        v(:,i) = v(:,i) + A * dt; 
    end

    % Update Endpoint Position
    theta = w * t;
    x(:,N) = [ R * cos(theta), R * sin(theta) , Z];

    % Update Interior Positions
    for i = 2: (N-1)
        x(:,i) = x(:,i) + v(:,i) * dt;
    end
  
    % Plot and Save
    scatter3(x(1,:),x(2,:),x(3,:))
    axis([- 3 * N, 3 * N,- 3 * N, 3 * N, 0, Z])

    currentFrame = getframe(gcf);
    writeVideo(writerObj1, currentFrame);

end

close(writerObj1);
