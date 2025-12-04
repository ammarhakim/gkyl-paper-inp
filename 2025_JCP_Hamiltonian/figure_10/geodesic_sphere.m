function [] = geodesic_sphere(R,theta0,phi0,theta_dot0,phi_dot0,tf)
% Initial conditions
% R = 1;
% theta0 = pi/2;        % Initial theta (latitude)
% phi0 = -pi/2;         % Initial phi (longitude)
% theta_dot0 = 1.5;     % Initial theta velocity
% phi_dot0 = 1;         % Initial phi velocity

% Combine initial conditions into a vector
y0 = [theta0; phi0; theta_dot0; phi_dot0];

% Time span for the simulation
tspan = [0 tf];

% Solve the system using ode45 with event detection for boundary condition
options = odeset('Events', @bounce_event,'AbsTol',1e-10);
tau = [];
sol = [];

% Loop for handling multiple bounces
while tspan(1) < tspan(2)
    [tau_segment, sol_segment, te, ye, ie] = ode45(@geodesic_ode, tspan, y0, options);
    tau = [tau; tau_segment];
    sol = [sol; sol_segment]; 

    % No more events, exit the loop
    if isempty(te)
        break;
    end

    % Reflect theta_dot and adjust theta if it goes out of bounds
    theta = ye(1);
    theta_dot = ye(3);
    if theta >= pi/2 + pi/4
        theta = (pi/2 + pi/4) - 1e-6;
    elseif theta <= pi/2 - pi/4
        theta = (pi/2 - pi/4) + 1e-6;
    end
    theta_dot = -theta_dot; 

    % Set new initial conditions for the next segment
    y0 = [theta, ye(2), theta_dot, ye(4) ];

    % Update the time span to continue the integration
    tspan(1) = tau_segment(end);
end

%The solutions
t = tau;
y = sol;

% Extract theta and phi from the solution
theta = y(:, 1);
phi = y(:, 2);

% Plot the solutions in 3D space
hold on;


% Overplot a semi-transparent spherical contour
[theta_grid, phi_grid] = meshgrid(linspace(pi/2 - pi/4, pi/2 + pi/4, 25), linspace(0, 2*pi, 50));
x = R * sin(theta_grid) .* cos(phi_grid);
y = R * sin(theta_grid) .* sin(phi_grid);
z = R * cos(theta_grid);

% Plot the semi-transparent surface
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.5, 'EdgeColor', 'black', 'FaceColor', 'yellow');
%title('Geodesics, Refl. Bound. at \theta = \pi/2 +/- \pi/4, (Red t = 0, Blue t = t_f) ');

x = R * sin(theta) .* cos(phi);
y = R * sin(theta) .* sin(phi);
z = R * cos(theta);

% Create a colormap from red to blue corresponding to sim-time
num_points = length(tau);
colors = [linspace(1, 0, num_points)', zeros(num_points, 1), linspace(0, 1, num_points)'];


% Plot the particle's trajectory with color change over time
for i = 1:num_points-1
    line([x(i), x(i+1)], [y(i), y(i+1)], [z(i), z(i+1)], ...
        'Color', colors(i,:), 'LineWidth', 2); 
end
xlabel('x [arb.]');
ylabel('y [arb.]');
zlabel('z [arb.]');
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])
grid off;
axis equal;
hold off;

% Set a view angle
view(45, 30); 

end

% Define the geodesic equations as a system of first-order ODEs
function dydt = geodesic_ode(t, y)

    % y(1) = theta, y(2) = phi, y(3) = theta_dot, y(4) = phi_dot
    theta = y(1);
    phi = y(2);
    theta_dot = y(3);
    phi_dot = y(4);
    
    % Radius of the sphere
    R = 1; 
    
    % The geodesic equations in first-order form
    dtheta_dt = theta_dot;
    dphi_dt = phi_dot;
    dtheta_dot_dt = phi_dot^2 * cos(theta) * sin(theta);
    dphi_dot_dt = -2 * phi_dot * theta_dot * cos(theta) / sin(theta);
    
    % Return the derivatives as a column vector
    dydt = [dtheta_dt; dphi_dt; dtheta_dot_dt; dphi_dot_dt];
end

% For reflecting BC's
function [value, isterminal, direction] = bounce_event(tau, X)
% This event function detects when z exceeds 1 or goes below -1
theta = X(1);

% Define the boundaries
lower_bound = pi/2 - pi/4;
upper_bound = pi/2 + pi/4;

% Event is triggered when theta is outside the boundaries
value = min(theta - lower_bound, upper_bound - theta);

% Stop the integration to handle the bounce
isterminal = 1; 

% The event is detected when z moves outside the boundary in either direction
direction = 0;
end
