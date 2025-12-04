function [] = geodesic_hyperbolic(theta0,z0,theta_dot0,z_dot0,tf)
% Define initial conditions for theta, z, and their derivatives
% theta0 = -pi/2;     % initial theta
% z0 = -0.5;          % initial z
% theta_dot0 = 1.5;   % initial theta_dot (dtheta/dtau)
% z_dot0 = 0.8;       % initial z_dot (dz/dtau)

% Pack the initial conditions into a vector
initial_conditions = [theta0, z0, theta_dot0, z_dot0];

% Time span
tau_span = [0, tf];

% Solve the system using ode45 with event detection for boundary condition
options = odeset('Events', @bounce_event, 'AbsTol',1e-10);
tau = [];
sol = [];

% Loop for handling multiple bounces
while tau_span(1) < tau_span(2)
    [tau_segment, sol_segment, te, ye, ie] = ode45(@geodesic_equations, tau_span, initial_conditions, options);
    tau = [tau; tau_segment];
    sol = [sol; sol_segment]; 

    % No more events, exit the loop
    if isempty(te)
        break;
    end

    % Reflect z_dot and adjust z if it goes out of bounds
    z = ye(2);
    z_dot = ye(4);
    if z >= 1
        z = 0.99999999;
    elseif z <= -1
        z = -0.99999999; 
    end
    z_dot = -z_dot; 

    % Set new initial conditions for the next segment
    initial_conditions = [ye(1), z, ye(3), z_dot];

    % Update the time span to continue the integration
    tau_span(1) = tau_segment(end);
end

% Extract solutions
theta = sol(:, 1);
z = sol(:, 2);
theta_dot = sol(:, 3);
z_dot = sol(:, 4);

% Plot the solutions in 3D space
hold on;

% Define the radius
R = 1; %

% Plot the R = 1 surface as a transparent cylinder
[Theta, Z] = meshgrid(linspace(0, 2*pi, 50), linspace(-1.0, 1.0, 50));
X_surf = sqrt(R^2 + Z.^2) .* cos(Theta);
Y_surf = sqrt(R^2 + Z.^2) .* sin(Theta);
surf(X_surf, Y_surf, Z, ...
    'FaceColor', [0, 0.4470, 0.7410], ... 
    'FaceAlpha', 0.3, ...
    'EdgeAlpha', 0.5, 'EdgeColor', 'black', ...
    'DisplayName', 'R = 1 Surface');

% Label axes and set title
xlabel('x [arb.]');
ylabel('y [arb.]');
zlabel('z [arb.]');
%title('Geodesics, Refl. Bound. at z = 1, (Red t = 0, Blue t = tf) ');
grid off;
axis equal;

% Convert to x,y 
x = sqrt(R^2 + z.^2) .* cos(theta);
y = sqrt(R^2 + z.^2) .* sin(theta);

% Create a colormap from red to blue
num_points = length(tau);
colors = [linspace(1, 0, num_points)', zeros(num_points, 1), linspace(0, 1, num_points)'];

% Plot the particle's trajectory with color change over time
for i = 1:num_points-1
    line([x(i), x(i+1)], [y(i), y(i+1)], [z(i), z(i+1)], ...
        'Color', colors(i,:), 'LineWidth', 2);
end

hold off;
xlim([-1.5,1.5])
zlim([-1.5,1.5])
ylim([-1.5,1.5])

% Set the view angle
view(45, 30);
end



function dXdt = geodesic_equations(tau, X)
% Unpack the variables from the input vector X
theta = X(1);
z = X(2);
theta_dot = X(3);
z_dot = X(4);

% Constants (set R according to the problem)
R = 1;

% Equations of motion
theta_dot_dot = -(2 * theta_dot * z * z_dot) / (R^2 + z^2);
z_dot_dot = (theta_dot^2 * z * (R^2 + z^2)) / (R^2 + 2 * z^2) ...
    - (R^2 * z * z_dot^2) / (R^4 + 3 * R^2 * z^2 + 2 * z^4);

% Pack the derivatives into a vector
dXdt = [theta_dot; z_dot; theta_dot_dot; z_dot_dot];
end

function [value, isterminal, direction] = bounce_event(tau, X)
% This event function detects when z exceeds 1 or goes below -1
z = X(2);

% Event is triggered when z is either greater than 1 or less than -1
value = 1 - abs(z);

% Stop the integration to handle the bounce
isterminal = 1; 

% The event is detected when z moves outside the boundary in either direction
direction = 0;
end