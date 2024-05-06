clear all; 

% Initilaize Quantities
N = 200; 
M = 1; % total mass
m = M / N; % mass per unit length
T = 1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Circle

x = zeros(2,N); % position vector
v = zeros(2,N); % velocity vector

R = 30;

for i = 1: N
    x(1,i) = R * cos(2 * pi * i / N);
    x(2,i) = R * sin(2 * pi * i / N);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arc Length at Rest

arc = 2 * pi * R / N;
l= arc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 1: single point lateral

% v(1, N) = 10;
% v(1, N/2) = -10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 2: distributed larger

% for i = 40 : 60 
%     v(2,i) = R / 80 * (i - 40) * (60 - i);
%     v(2, N/2 + i) = - R / 80 * (i - 40) * (60 - i);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 3: distributed smaller

for i = 40 : 60 
    v(2,i) = R / 200 * (i - 40) * (60 - i);
    v(2, N/2 + i) = - R / 200 * (i - 40) * (60 - i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate Motion

tmax = 50; 
dt = 0.05; 
clockmax = tmax / dt;


writerObj1 = VideoWriter('Circle.avi');
open(writerObj1);

scatter(x(1,:),x(2,:))
axis([- 0.4 * N, 0.4 * N,- 0.4 * N, 0.4 * N])
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
    axis([- 0.4 * N, 0.4 * N,- 0.4 * N, 0.4 * N])

    currentFrame = getframe(gcf);
    writeVideo(writerObj1, currentFrame);

end

close(writerObj1);


