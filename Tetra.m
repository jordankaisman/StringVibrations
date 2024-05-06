clear all; 

% Initilaize Quantities
N = 200; 
M = 1; % total mass
m = M / N; % mass per unit length
T = 5; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tetracuspid: we start with a circle and then invert it

x = zeros(2,N); % position vector
v = zeros(2,N); % velocity vector

R = 30;

for i = 1: N
    x(1,i) = R * cos(2 * pi * i / N);
    x(2,i) = R * sin(2 * pi * i / N);
end

% Reflection Matrices
reflect1 = [0 1; 1 0];
reflect2 = [0 -1; -1 0];


% Invert the Circle
for i = 1: N / 4
    x(:,i) = reflect2 * x(:,i);
    x(1,i) = x(1,i) + R;
    x(2,i) = x(2,i) + R;
end

for i = N / 4: N / 2
    x(:,i) = reflect1 * x(:,i);
    x(1,i) = x(1,i) - R;
    x(2,i) = x(2,i) + R;
end

for i = N / 2: 3 /4 * N
    x(:,i) = reflect2 * x(:,i);
    x(1,i) = x(1,i) - R;
    x(2,i) = x(2,i) - R;
end

for i = 3 /4 * N: N 
    x(:,i) = reflect1 * x(:,i);
    x(1,i) = x(1,i) + R;
    x(2,i) = x(2,i) - R;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arc Length at Rest

arc = 2 * pi * R / N;
l = arc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 1: single point outward on all sides

v(1, N / 8) = 20;
v(2, N / 8) = 20;

v(1, 3 * N / 8) = -20;
v(2, 3 * N / 8) = 20;

v(1, 5 * N / 8) = -20;
v(2, 5 * N / 8) = -20;

v(1, 7 * N / 8) = 20;
v(2, 7 * N / 8) = -20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 2: single point inward on all sides

% v(1, N / 8) = -20;
% v(2, N / 8) = -20;
% 
% v(1, 3 * N / 8) = 20;
% v(2, 3 * N / 8) = -20;
% 
% v(1, 5 * N / 8) = 20;
% v(2, 5 * N / 8) = 20;
% 
% v(1, 7 * N / 8) = -20;
% v(2, 7 * N / 8) = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate Motion

tmax = 50; 
dt = 0.005; 
clockmax = tmax / dt;


writerObj1 = VideoWriter('Tetra.avi');
open(writerObj1);

scatter(x(1,:),x(2,:))
axis([- N / 4, N / 4,- N / 4, N / 4])
axis manual
currentFrame = getframe(gcf);
writeVideo(writerObj1, currentFrame);

for clock = 1: clockmax
    t = clock * dt;

    % Update Velocity 

    for i = 1: N
        % Find Neighbors (modulo in case we are at index 1 or N)
        n1 = mod(i,N) + 1; % N->1
        n2 = mod(i-2,N) + 1; % 1->N
        
        D1 = x(:,n1) - x(:,i);
        D2 = x(:,n2) - x(:,i);
        E1 = max([norm(D1) - l, 0]); % excess arc length 1
        E2 = max([norm(D2) - l, 0]); % excess arc length 2
        F1 = T * E1 * D1 / norm(D1);
        F2 = T * E2 * D2 / norm(D2);
        A = (F1 + F2) / m;

        v(:,i) = v(:,i) + A * dt; 
    end

    % Update Position

    for i = 1: N
        x(:,i) = x(:,i) + v(:,i) * dt;
    end
      
    % Plot and Save
    scatter(x(1,:),x(2,:))
    axis([- N / 4, N / 4,- N / 4, N / 4])

    currentFrame = getframe(gcf);
    writeVideo(writerObj1, currentFrame);

end

close(writerObj1);


