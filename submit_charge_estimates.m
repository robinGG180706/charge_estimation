% % EPFL EE466 - Fabrizio Sossan
% % 
% % Demo script to obtain cell measurements (voltage and current) and submit your
% % charge (Ah) estimates through the endpoint "url_base/utt/submit/[...]".
% %
% clc;
% clear all;
% close all;
% 
% 
% % Base url (STRING)
% url_base = 'http://72.60.84.102:8080/';
% % Your team number (INTEGER) and password (STRING)
% team_number = 16;
% pwd = 'E7F8A9B0';
% 
% 
% % Dataset number (INTEGER). Allowed values are: 1, 2, or 3.
% % ---
% % VERY IMPORTANT NOTICE ON DATASETS - reminder
% % ---
% % There are three datasets: 1, 2, 3. Dataset 1 is for testing: you can
% % read and write data as many times as you want. Datasets 2 and 3 can be
% % written only once, in a sequential fashion (*), implying that you
% % can't take a peek at values in the future. 
% 
% % (*) Sequential fashion means that the variable <sample_index> here below 
% % should be submitted with the following values (once at a time): 0, 1, 2, 3, ...,
% % exactly as in the example provided below. With Dataset 1, you can rewind,
% % with datasets 2 and 3 you can't.
% 
% % Test your charge estimator on Dataset 1. Once you are ready, move to the
% % two other datasets. As it should be crystal clear by now, with datasets
% % 1, and 2 you have one shot per time interval. To adjust the dataset you
% % are using, adapt the variable <dataset_number> in the following code line.
% % ---
% dataset_number = 1;
% 
% 
% 
% % Inline function to generate a suitably formatted URL for the endpoint:
% % url_base/utt/submit/<team_number>/<pwd>/<dataset_number>/<sample_index>/<charge_estimate>
% compose_url = @(dataset_number, sample_index, charge_estimate) sprintf('%s/utt/submit/%d/%s/%d/%d/%.2f', url_base, team_number, pwd, ...
%     dataset_number, sample_index, charge_estimate);
% 
% 
% 
% % Empty current and voltage vector 
% current_vector = [];
% voltage_vector = [];
% % Charge estimate vector
charge_estimate_vector = [-1];
error_vector = [];
voltage_estimator_vector = [];
% Number of samples per dataset
% Dataset 1: 405 available samples.
% Dataset 2: 661 available samples.
% Dataset 3: 110 available samples.

% The first sample index is 0, the last index is what reported here above.

% Parameter
% =============================================================

Ts = 2.0/3600;         % Abtastzeit [h]
tau = 5.223e4;              % RC-Zeitkonstante [s]
R0  = 0.0874;            % ohmscher Widerstand
R1  = 4.4;            % RC-Widerstand

% Diskrete Modellierung
a = exp(-Ts/tau);          % Diskreter a-Parameter
b = R1 * (1 - a);          % Diskreter b-Parameter
% ============================================================
% Luenberger Observer Setup
% =============================================================

% Zustände: x = [SOC; Vrc]
x_hat = zeros(2,1);
x_hat(1) = -1;

% Beobachter-Gain L
L1 = 1e-3;        % kleiner SOC-Gain
lambda_RC = 0.9; % gewünschter RC-Pole
L2 = lambda_RC - a;


L = [L1; L2];

R0est = [];
for sample_index=1:405
    % Compose appropriate endpoint
    % url_endpoint = compose_url(dataset_number, sample_index, charge_estimate_vector(end));
    % Read the endpoint, and decode json data
    % data_txt = webread(url_endpoint, weboptions("ContentType","text"));
    % try
    %     data = jsondecode(data_txt);
    % catch ME
    %     error(sprintf("Couldn't parse json. Error message: %s", data_txt));
    % end
    % 
    % fprintf('Sample number %d requested OK.\n', sample_index);

    % Extract the current and voltage measurements delivered by the endpoint
    voltage_value = voltage_vector(sample_index);
    current_value = current_vector(sample_index);

    if(x_hat(1)<0)
        x_hat(1) = voltage_to_charge(voltage_value - R0*current_value);
    end
    if(abs(current_vector(end-1)-current_vector(end)))
        R0est = vercat(R0est,abs((voltage_vector(end-1)-voltage_vector(end))/(current_vector(end-1)-current_vector(end))));
    end
    % --- Messfehler / Innovation ---
    V_est = charge_to_voltage(x_hat(1)) + R0*current_value - x_hat(2);
    % if(sample_index==1)
    %     V_est = V_est_new;
    % end
    % if (V_est_new - V_est > 0.05)
    %     V_est = V_est_new;
    % end
    e = voltage_value - V_est;
    error_vector = vertcat(error_vector, e);
    voltage_estimator_vector = vertcat(voltage_estimator_vector, V_est);
    % --- Zustandspdate ---
    SOC_next = x_hat(1) + Ts*current_value + L(1)*e;
    Vrc_next = a*x_hat(2) - b*current_value + L(2)*e;
    
    if ((sample_index > 111) && (sample_index < 150))
        fprintf("x_hat(2) : %d\n", x_hat(2));
    end
    % --- Speichern + begrenzen ---
    x_hat(1) = SOC_next;
    x_hat(2) = Vrc_next;

    % SOC_est(k) = x_hat(1);
    Vrc_est(sample_index) = Vrc_next;
    

    % ----
    % TODO Computation of the charge estimate. Your show, now.
    % ----
    charge_estimate = SOC_next;

    % --- End of Computation of the charge estimate


    % Append the charge estimate to the vector. It will be sent out at the
    % next iteration of the for loop.
    charge_estimate_vector = vertcat(charge_estimate_vector, charge_estimate);
end

voltage_estimator_vector2 = zeros(110,1);
error_vector2 = zeros(110,1);
charge_estimate_vector2 = zeros(110,1);
x_hat = zeros(2,1);
x_hat(1) = -1;

L1 = 1e-3;        % kleiner SOC-Gain
lambda_RC = 0.9; % gewünschter RC-Pole
L2 = lambda_RC - a;

L = [L1; L2];

for sample_index=111:405
    % Compose appropriate endpoint
    % url_endpoint = compose_url(dataset_number, sample_index, charge_estimate_vector(end));
    % Read the endpoint, and decode json data
    % data_txt = webread(url_endpoint, weboptions("ContentType","text"));
    % try
    %     data = jsondecode(data_txt);
    % catch ME
    %     error(sprintf("Couldn't parse json. Error message: %s", data_txt));
    % end
    % 
    % fprintf('Sample number %d requested OK.\n', sample_index);

    % Extract the current and voltage measurements delivered by the endpoint
    voltage_value = voltage_vector(sample_index);
    current_value = current_vector(sample_index);

    if(x_hat(1)<0)
        x_hat(1) = voltage_to_charge(voltage_value - R0*current_value);
    end

    % --- Messfehler / Innovation ---
    V_est = charge_to_voltage(x_hat(1)) + R0*current_value - x_hat(2);
    % if(sample_index==1)
    %     V_est = V_est_new;
    % end
    % if (V_est_new - V_est > 0.05)
    %     V_est = V_est_new;
    % end
    e = voltage_value - V_est;
    error_vector2 = vertcat(error_vector2, e);
    voltage_estimator_vector2 = vertcat(voltage_estimator_vector2, V_est);
    % --- Zustandspdate ---
    SOC_next = x_hat(1) + Ts*current_value + L(1)*e;
    Vrc_next = a*x_hat(2) - b*current_value + L(2)*e;
    if(sample_index == 111)
        fprintf("test 2\n");
    end
    if ((sample_index > 111) && (sample_index < 150))
        fprintf("x_hat(2) : %d\n", x_hat(2));
    end

    % --- Speichern + begrenzen ---
    x_hat(1) = SOC_next;
    x_hat(2) = Vrc_next;

    % SOC_est(k) = x_hat(1);
    Vrc_est(sample_index) = Vrc_next;
    

    % ----
    % TODO Computation of the charge estimate. Your show, now.
    % ----
    charge_estimate = SOC_next;

    % --- End of Computation of the charge estimate


    % Append the charge estimate to the vector. It will be sent out at the
    % next iteration of the for loop.
    charge_estimate_vector2 = vertcat(charge_estimate_vector2, charge_estimate);
end


% Plot data
% Sampling time is 2 s! 
Ts = 2;
time_vector = (1:numel(voltage_vector))' * Ts;

subplot(4, 1, 1)
plot(time_vector, voltage_vector)
hold on
plot(time_vector, voltage_estimator_vector, '-r')
hold on
plot(time_vector, voltage_estimator_vector2, '-g')
hold off
xlabel('Time [seconds]')
ylabel('Voltage [V]')

subplot(4, 1, 2)
plot(time_vector, current_vector)
xlabel('Time [seconds]')
ylabel('Current [A]')


subplot(4, 1, 3)
plot(time_vector, charge_estimate_vector(2:end))
hold on
plot(time_vector, charge_estimate_vector2, '-r')
hold off
xlabel('Time [seconds]')
ylabel('Charge level [Ah]')

% Delete avant de rendre

subplot(4, 1, 4)
plot(time_vector, error_vector)
hold on
plot(time_vector, error_vector2, '-g')
hold on
yline(0, 'r');
hold off
xlabel('Time [seconds]')
ylabel('error')

error_1 = (charge_estimate_vector(end)-voltage_to_charge(voltage_vector(end)-R0*current_vector(end)+x_hat(2)))*1000;
fprintf("Fin : Erreur entre le modèle et l'estimation : %.2f mAh \n", error_1);
R0estfin = mean(R0est);

