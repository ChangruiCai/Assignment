%% Related to Assignment 1
% Number of nodes
N = 3;
ndof = N * 2; % number of degrees of freedom
dt = 0.00001; % second - Time step size
RodLength = 0.1; % l = 1 m
deltaL = RodLength / (N-1);
% Radii of spheres
R1 = 0.005; % meter
R2 = 0.025; % meter
R3 = 0.005; % meter
% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
totalTime = 0.4; % second - total simulation time
% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
EA = Y * pi * r0^2; % Newton
% Geometry - initial configuration
nodes = zeros(N,2);
for c=1:N % Loop over all the nodes
    nodes(c,1) = (c-1) * deltaL; % x coordinates
    nodes(c,2) = 0;
end
% Mass, M
M = zeros(ndof, ndof);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;
M_inv = inv(M); % Precompute the inverse of mass matrix
% Viscous damping matrix, C
C = zeros(6,6);
C1 = 6 * pi * visc * R1;
C2 = 6 * pi * visc * R2;
C3 = 6 * pi * visc * R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;
% Weight vector, W
W = zeros(ndof, 1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;
% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end
u0 = zeros(ndof, 1); % initial velocity
% Time marching scheme
Nsteps = round(totalTime/dt);
% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
for c = 2:Nsteps
    fprintf('Time = %f\n', (c-1) * dt);
    
    % Calculate forces
    f = zeros(size(q0));  % Reset forces
    
    % Elastic forces
    % Linear spring between nodes 1 and 2
    xk = q0(1);
    yk = q0(2);
    xkp1 = q0(3);
    ykp1 = q0(4);
    l_k = deltaL;
    dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
    f(1:4) = f(1:4) + dF;

    % Linear spring between nodes 2 and 3
    xk = q0(3);
    yk = q0(4);
    xkp1 = q0(5);
    ykp1 = q0(6);
    l_k = deltaL;
    dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
    f(3:6) = f(3:6) + dF;

    % Bending spring at node 2
    xkm1 = q0(1);
    ykm1 = q0(2);
    xk = q0(3);
    yk = q0(4);
    xkp1 = q0(5);
    ykp1 = q0(6);
    curvature0 = 0;
    l_k = deltaL;
    dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
    f(1:6) = f(1:6) + dF;

    % Viscous force
    f = f - C * u0;

    % Weight (gravity or external force)
    f = f + W;

    % Explicit update using forward Euler method
    % Calculate acceleration a = M_inv * f
    a = M_inv * f;

    % Update velocity: u(t+dt) = u(t) + dt * a
    u = u0 + dt * a;

    % Update position: q(t+dt) = q(t) + dt * u
    q = q0 + dt * u;

    % Store some information
    all_mid_v(c) = u(4);
    % Plot
    figure(1);
    plot(q(1:2:end), q(2:2:end), 'ro-');
    axis equal
    xlabel('x [meter]');
    ylabel('y [meter]');
    drawnow

    % Update (new becomes old)
    q0 = q;
    u0 = u;
end
% Plot middle node downward velocity
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of middle node, v [m/s]');
