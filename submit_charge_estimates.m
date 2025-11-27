% EPFL EE466 - Fabrizio Sossan
% 
% Demo script to obtain cell measurements (voltage and current) and submit your
% charge (Ah) estimates through the endpoint "url_base/utt/submit/[...]".
%
clc;
clear all;
close all;


% Base url (STRING)
url_base = 'http://72.60.84.102:8080/';
% Your team number (INTEGER) and password (STRING)
team_number = 16;
pwd = 'E7F8A9B0';


% Dataset number (INTEGER). Allowed values are: 1, 2, or 3.
% ---
% VERY IMPORTANT NOTICE ON DATASETS - reminder
% ---
% There are three datasets: 1, 2, 3. Dataset 1 is for testing: you can
% read and write data as many times as you want. Datasets 2 and 3 can be
% written only once, in a sequential fashion (*), implying that you
% can't take a peek at values in the future. 

% (*) Sequential fashion means that the variable <sample_index> here below 
% should be submitted with the following values (once at a time): 0, 1, 2, 3, ...,
% exactly as in the example provided below. With Dataset 1, you can rewind,
% with datasets 2 and 3 you can't.

% Test your charge estimator on Dataset 1. Once you are ready, move to the
% two other datasets. As it should be crystal clear by now, with datasets
% 1, and 2 you have one shot per time interval. To adjust the dataset you
% are using, adapt the variable <dataset_number> in the following code line.
% ---
dataset_number = 1;



% Inline function to generate a suitably formatted URL for the endpoint:
% url_base/utt/submit/<team_number>/<pwd>/<dataset_number>/<sample_index>/<charge_estimate>
compose_url = @(dataset_number, sample_index, charge_estimate) sprintf('%s/utt/submit/%d/%s/%d/%d/%.2f', url_base, team_number, pwd, ...
    dataset_number, sample_index, charge_estimate);



% Empty current and voltage vector 
current_vector = [];
voltage_vector = [];
% Charge estimate vector
charge_estimate_vector = [-1];

% Number of samples per dataset
% Dataset 1: 405 available samples.
% Dataset 2: 661 available samples.
% Dataset 3: 110 available samples.

% The first sample index is 0, the last index is what reported here above.

% Parameter
% =============================================================

Ts = 2.0/3600;         % Abtastzeit [h]
tau = 60;              % RC-Zeitkonstante [s]
R0  = 0.08;            % ohmscher Widerstand
R1  = 0.02;            % RC-Widerstand

% Diskrete Modellierung
a = exp(-Ts/tau);          % Diskreter a-Parameter
b = R1 * (1 - a);          % Diskreter b-Parameter
% ============================================================
% Luenberger Observer Setup
% =============================================================

% Zustände: x = [SOC; Vrc]
x_hat = zeros(2,1);

% Beobachter-Gain L
L1 = 1e-3;        % kleiner SOC-Gain
lambda_RC = 0.95; % gewünschter RC-Pole
L2 = lambda_RC - a;

L = [L1; L2];


for sample_index=1:405
    % Compose appropriate endpoint
    url_endpoint = compose_url(dataset_number, sample_index, charge_estimate_vector(end));
    % Read the endpoint, and decode json data
    data_txt = webread(url_endpoint, weboptions("ContentType","text"));
    try
        data = jsondecode(data_txt);
    catch ME
        error(sprintf("Couldn't parse json. Error message: %s", data_txt));
    end

    fprintf('Sample number %d requested OK.\n', sample_index);

    % Extract the current and voltage measurements delivered by the endpoint
    voltage_value = data.voltage;
    current_value = data.current;

    if(sample_index==1)
        x_hat(1) = voltage_to_charge(voltage_value);
    end

    % and append them to a vector to keep an history
    voltage_vector = vertcat(voltage_vector, data.voltage);
    current_vector = vertcat(current_vector, data.current);

    % --- Messfehler / Innovation ---
    V_est = voltage_to_charge(x_hat(1)) + R0*current_value - x_hat(2);
    e = voltage_value - V_est;

    % --- Zustandspdate ---
    SOC_next = x_hat(1) + Ts*current_value + L(1)*e;
    Vrc_next = a*x_hat(2) - b*current_value + L(2)*e;

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


% Plot data
% Sampling time is 2 s! 
Ts = 2;
time_vector = (1:numel(voltage_vector))' * Ts;

subplot(3, 1, 1)
plot(time_vector, voltage_vector)
xlabel('Time [seconds]')
ylabel('Voltage [V]')

subplot(3, 1, 2)
plot(time_vector, current_vector)
xlabel('Time [seconds]')
ylabel('Current [V]')


subplot(3, 1, 3)
plot(time_vector, charge_estimate_vector(2:end))
xlabel('Time [seconds]')
ylabel('Charge level [Ah]')


