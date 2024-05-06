clear all;

% Initilaize Quantities
N = 200; 
M = 1; % total mass
m = M / N; % mass per unit length
T = 1; % tension force per unit distance between nodes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shape: square with N/4 - 1 edge nodes per side 
% Indexing starts at bottom right corner

x = zeros(2,N); % position vector
v = zeros(2,N); % velocity vector
for i = 1: N
    if (i < N / 4)
        x(1,i) = N / 4;
        x(2,i) = i;
    end
    if (N / 4 <= i && i < N / 2)
        x(1,i) = N / 4 - (i - N / 4);
        x(2,i) = N / 4;
    end
    if (N / 2 <= i && i < N * 3 / 4 )
        x(1,i) = 0;
        x(2,i) = N / 4 - (i - N / 2);
    end
    if (N * 3 / 4 <= i)
        x(1,i) = i- (3 * N/ 4);
        x(2,i) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arc Length at Rest
l = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 1: single point lateral only

% v(1,ceil(N / 8)) = 10;
% v(1,ceil(5 * N / 8)) = -10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 2: single point lateral and vertical

% v(1,ceil(N / 8)) = 10;
% v(2,ceil(3 * N / 8)) = 10;
% v(1,ceil(5 * N / 8)) = -10;
% v(2,ceil(7 * N / 8)) = -10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 3: lateral distributed

% for i = 16: 34 
%     v(1,i) = v(1,1) + 1 / 20 * (i - 16) * (34 - i);
%     v(1, N/2 + i) = v(1, N/2 + i) - 1 / 20 * (i - 16) * (34 - i);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbation Choice 4: single point corners

v(1,N/4) = 10;
v(2,N/4) = 10;
v(1,N/2) = -10;
v(2,N/2) = 10;
v(1,3 * N/4) = -10;
v(2,3 * N/4) = -10;
v(1,N) = 10;
v(2,N) = -10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate Motion

tmax = 50; 
dt = 0.05; 
clockmax = tmax / dt;


writerObj1 = VideoWriter('Square.avi');
open(writerObj1);

scatter(x(1,:),x(2,:))
axis([- 0.25 * N, 0.5 * N,- 0.25 * N, 0.5 * N])
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
    axis([- 0.25 * N, 0.5 * N,- 0.25 * N, 0.5 * N])

    currentFrame = getframe(gcf);
    writeVideo(writerObj1, currentFrame);

end

close(writerObj1);


