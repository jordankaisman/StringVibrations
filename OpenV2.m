
% Endpoints move in same direction


% Initilaize Quantities
N = 101;
M = 1; % total mass
m = M / N; % mass per unit length (uniform)
T = 1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curve: straight line between endpoints

x = zeros(3,N); % position vector
v = zeros(3,N); % velocity vector
R1 = (N - 1) / 25;
R2 = R1;
Z = 2 * N; %depth

x(:,1) = [-R2,0,0]; 
x(:,N) = [R1,0,Z];

for i = 2: (N-1)
    x(:,i) = x(:,1) + (i-1) / (N-1) * (x(:,N)- x(:,1));
end

w = 2; % angular velocity of endpoints 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arc Length at Rest
l = 0.5 * norm(x(:,1)-x(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate Motion

tmax = 200; 
dt = 0.05; 
clockmax = tmax / dt;


%Video Setup
writerObj1 = VideoWriter('Open_V2.avi');
open(writerObj1);

% Initial Frame
scatter3(x(1,:),x(2,:),x(3,:))
title('Insert Title')
k = 1; % axis scale factor
axis([- k * N, k * N,- k * N, k * N, 0, 2 * N])
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
    x(:,1) = [ - R2 * cos(theta), - R2 * sin(theta) , 0];
    x(:,N) = [ R1 * cos(theta), R1 * sin(theta) , Z]; 


    % Update Interior Positions
    for i = 2: (N-1)
        x(:,i) = x(:,i) + v(:,i) * dt;
    end
  
    % Plot and Save
    scatter3(x(1,:),x(2,:),x(3,:))
    axis([- k * N, k * N,- k * N, k * N, 0, 2 * N])

    currentFrame = getframe(gcf);
    writeVideo(writerObj1, currentFrame);

end

close(writerObj1);
