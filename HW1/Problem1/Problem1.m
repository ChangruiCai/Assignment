%% Related to Assignment 1
% Number of nodes
N = 3;
ndof = N * 2; % number of degrees of freedom
dt = 0.01; % second - Time step size
RodLength = 0.1; % l = 1 m
deltaL = RodLength / (N-1);
% Radii of spheres
R1 = 0.005; % meter
R2 = 0.025; % meter
R3 = 0.005; % meter

%initial setting of 2nd ball
cycle=2;
ball2position(1)=0;

% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
totalTime = 10; % second - total simulation time
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
u0 = zeros(ndof, 1); % old velocity (initial velocity)
% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected
% Time marching scheme
Nsteps = round(totalTime/dt);
% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
for c = 2:Nsteps
fprintf('Time = %f\n', (c-1) * dt);
% Guess
q = q0; % New DOFs are initialized to be equal to old DOFs

% Newton Raphson
err = 10 * tol;
while err > tol
f = M / dt * ( (q-q0)/dt - u0 );
J = M / dt^2;
%
% Elastic forces
%
% Linear spring between nodes 1 and 2
xk = q(1);
yk = q(2);
xkp1 = q(3);
ykp1 = q(4);
l_k = deltaL;
dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
f(1:4) = f(1:4) + dF;
J(1:4,1:4) = J(1:4,1:4) + dJ;
% Linear spring between nodes 2 and 3
xk = q(3);
yk = q(4);
xkp1 = q(5);
ykp1 = q(6);
l_k = deltaL;
dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
f(3:6) = f(3:6) + dF;
J(3:6, 3:6) = J(3:6, 3:6) + dJ;
% Bending spring at node 2
xkm1 = q(1);
ykm1 = q(2);
xk = q(3);
yk = q(4);
xkp1 = q(5);
ykp1 = q(6);
curvature0 = 0;
l_k = deltaL;
dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
f(1:6) = f(1:6) + dF;
J(1:6, 1:6) = J(1:6, 1:6) + dJ;
% Viscous force
f = f + C * (q-q0) / dt;
J = J + C / dt;
% Weight
f = f - W;
% Update
q = q - J \ f;
err = sum ( abs(f) );
end



% New velocity
u = (q - q0) / dt;
% Store some information
all_mid_v(c) = u(4);
% Plot
figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
xlabel('x [meter]');
ylabel('y [meter]');
drawnow
ball2position(cycle)=q(4);
cycle=cycle+1;


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


% plot position of middle node
figure(3);
timeArray = (1:Nsteps) * dt;
plot(timeArray, ball2position, 'k-');
xlabel('Time, t [sec]');
ylabel('Position (vertical) of middle node, d [m]');