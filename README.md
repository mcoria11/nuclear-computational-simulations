# nuclear-computational-simulations
function star_structure

    % Constants
    G = 6.673e-8; % Gravitational constant in dyne-cm^2/g^2
    c = 2.997e10; % Speed of light in cm/s
    k_B = 1.380649e-16; % Boltzmann constant in cgs units
    T = 1e7; % Temperature in K
    mu = 1.0; % Mean molecular weight
    m_p = 1.6726219e-24; % Proton mass in g

    % Initial conditions
    R0 = 1e-2; % Initial radius in cm
    P0 = 1.0e19; % Central pressure in dyne/cm^2
    M0 = 1e-10; % Small initial mass in g to avoid singularities

    % Solar units for output
    R_sun = 6.96e10; % Solar radius in cm
    M_sun = 1.99e33; % Solar mass in g

    % Solve the ODEs using ode45
    opts = odeset('Events', @(r, y) pressure_zero(r, y));
    [r, y] = ode45(@(r, y) odefun(r, y, G, k_B, T, mu, m_p, c), [R0, 1e9], [M0, P0], opts);

    % Extract solutions
    M = y(:, 1);
    P = y(:, 2);

    % Convert to solar units
    r_solar_units = r / R_sun;
    M_solar_units = M / M_sun;

    % Display results
    fprintf('Maximum Mass: %.2f M_sun\n', max(M_solar_units));
    fprintf('Maximum Radius: %.2f R_sun\n', max(r_solar_units));

    % Plot results 2D
    figure;
    subplot(2, 1, 1);
    plot(r_solar_units, P, 'r', 'LineWidth', 2);
    xlabel('Radius ($R_\odot$)', 'Interpreter', 'latex');
    ylabel('Pressure (dyne/cm$^2$)', 'Interpreter', 'latex');
    title('Pressure vs Radius');
    grid on;

    subplot(2, 1, 2);
    plot(r_solar_units, M_solar_units, 'LineWidth', 2);
    xlabel('Radius ($R_\odot$)', 'Interpreter', 'latex');
    ylabel('Mass ($M_\odot$)', 'Interpreter', 'latex');
    title('Mass vs Radius');
    grid on;

    % Plot results 3D
    figure;
    scatter3(r_solar_units, P, M_solar_units, 36, r_solar_units, 'filled');
    xlabel('Radius ($R_\odot$)', 'Interpreter', 'latex');
    ylabel('Pressure (dyne/cm$^2$)', 'Interpreter', 'latex');
    zlabel('Mass ($M_\odot$)', 'Interpreter', 'latex');
    title('3D Plot of Star Structure');
    colormap jet;
    colorbar;
    caxis([min(r_solar_units) max(r_solar_units)]);
    grid on;
    set(gcf, 'Color', 'w');
    view(3); % 3D view
    rotate3d on; % Enable interactive rotation

    % Enhance 3D visualization
    lighting phong;
    camlight('headlight');
    material shiny;
    set(gca, 'Projection', 'perspective');

end

function dydr = odefun(r, y, G, k_B, T, mu, m_p, c)
    % Extract variables
    M = y(1);
    P = y(2);

    % Equation of state: energy_density = P / (k_B * T / (mu * m_p))
    if P > 0
        energy_density = P / (k_B * T / (mu * m_p));
        rho = energy_density / c^2;

        dMdr = 4 * pi * r^2 * rho;
        dPdr = -G * (rho + P/c^2) * (M + 4 * pi * r^3 * P/c^2) / (r^2 * (1 - 2 * G * M / (c^2 * r)));
    else
        dMdr = 0;
        dPdr = 0;
    end

    % Output derivative
    dydr = [dMdr; dPdr];
end

function [value, isterminal, direction] = pressure_zero(r, y)
    % Event function to stop integration when pressure reaches a small value
    value = y(2) - 1e-3; % Stop when pressure is close to zero
    isterminal = 1; % Stop the integration
    direction = -1; % Negative direction
end
