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



function battery_joint_ekf()
    % --- 1. Initialisierung ---
    dt = 2; % Abtastrate in Sekunden
    
    % Initialer Zustandsvektor x: 
    % [SOC; Vc1; Vc2; Rs; R1; C1; R2; C2]
    % WICHTIG: Gute Startwerte für Parameter raten, sonst divergiert der Filter!
    x_hat = [0; 0; 0; 0.01; 0.005; 100; 0.005; 300]; 
    
    % Kovarianzmatrix P (Unsicherheit der Schätzung)
    P = eye(8) * 1e-3; 
    
    % Prozessrauschen Q (Wie sehr vertrauen wir dem Modell?)
    % Parameter haben kleines Rauschen (Random Walk), SOC folgt Strom
    diag_Q = [1e-6, 1e-5, 1e-5, 1e-8, 1e-8, 1e-6, 1e-8, 1e-6]; 
    Q = diag(diag_Q);
    
    % Messrauschen R (Wie sehr vertrauen wir dem Spannungssensor?)
    R = 0.01; % z.B. 0.01 V^2 Varianz

    % Lade OCV Kurve (Beispiel)
    % soc_axis = 0:0.01:1;
    % ocv_values = ...; % Ihre Daten hier

    % --- Simulation Loop (oder Echtzeit-Daten) ---
    % Angenommen 'measurements' ist eine Matrix [Zeit, Strom, Spannung]
    % for k = 1:length(measurements)
    %     uk = measurements(k, 2); % Strom
    %     zk = measurements(k, 3); % Gemessene Spannung (v)
        
        % In der Praxis rufen Sie hier die Step-Funktion auf:
        % [x_hat, P] = ekf_step(x_hat, P, uk, zk, dt, Q, R);
        
        % Speichern der Ergebnisse...
    % end
end

function [x_new, P_new] = ekf_step(x, P, i, v_meas, dt, Q, R)
    % Entpacken der Zustände zur besseren Lesbarkeit
    soc = x(1); v1 = x(2); v2 = x(3);
    Rs = exp(x(4)); R1 = exp(x(5)); C1 = exp(x(6)); R2 = exp(x(7)); C2 = exp(x(8));

    % --- A. Prediction Step (Time Update) ---
    
    % 1. State Prediction (A Priori)
    % SOC Integration
    Ah_pred = soc + dt * i; 
    
    % Analytische Lösung für RC-Glieder (genauer als Euler)
    exp1 = exp(-dt / (R1 * C1));
    v1_pred = v1 * exp1 + R1 * (1 - exp1) * i;
    
    exp2 = exp(-dt / (R2 * C2));
    v2_pred = v2 * exp2 + R2 * (1 - exp2) * i;
    
    % Parameter bleiben konstant (Random Walk)
    % x_pred = [Ah_pred; v1_pred; v2_pred; Rs; R1; C1; R2; C2];
    x_pred = state_transition(x, i, dt);

    % 2. Jacobian F (Linearisierung des Zustandsübergangs) berechnen
    % Wir brauchen df/dx. Da Parameter im Zustandsvektor sind, ist das komplex.
    % Hier vereinfacht oder numerisch berechnen.
    F = compute_jacobian_F(x, i, dt); 
    
    % 3. Covariance Prediction
    P_pred = F * P * F' + Q;

    % --- B. Correction Step (Measurement Update) ---
    
    % 1. Vorhersage der Messung (Spannung)
    ocv_val = voltage_to_charge(Ah_pred); % Ihre OCV-Lookup Funktion
    v_pred = ocv_val - v1_pred - v2_pred - i * Rs;
    
    % 2. Jacobian H (Linearisierung der Messgleichung)
    % H = [dOCV/dSOC, -1, -1, -i, 0, 0, 0, 0]
    docv_dsoc = get_ocv_derivative(Ah_pred); % Ableitung der OCV Kurve
    H = [docv_dsoc, -1, -1, -i, 0, 0, 0, 0];
    
    % 3. Kalman Gain
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    
    % 4. Update State Estimate
    measurement_residual = v_meas - v_pred;
    x_new = x_pred + K * measurement_residual;

    % Parameter Constraints erzwingen (Widerstände und Kapazitäten müssen positiv sein)
    % Indizes: 4=Rs, 5=R1, 6=C1, 7=R2, 8=C2
    % min_R = 1e-4;    % 0.1 mOhm
    % min_C = 1.0;     % 1 Farad
    
    % x_new(4) = max(x_new(4), min_R); % Rs
    % x_new(5) = max(x_new(5), min_R); % R1
    % x_new(6) = max(x_new(6), min_R); % C1
    % x_new(7) = max(x_new(7), min_R); % R2
    % x_new(8) = max(x_new(8), min_R); % C2
    
    % keine ladung <0
    x_new(1) = max(x_new(1), 0);
    
    % 5. Update Covariance
    P_new = (eye(8) - K * H) * P_pred;
end

% --- Hilfsfunktionen ---
function F = compute_jacobian_F(x, i_input, dt)
    % Anzahl der Zustände (hier 8)
    nx = length(x);
    F = zeros(nx, nx);
    
    % Wahl der Schrittweite (epsilon)
    % Wichtig: Nicht zu klein (Rundungsfehler) und nicht zu groß (Linearisierungsfehler)
    % 1e-6 bis 1e-8 ist meist gut für 'double' Precision.
    epsilon = 1e-7;
    
    % 1. Berechnung des nominellen nächsten Zustands (ohne Störung)
    x_next_nominal = state_transition(x, i_input, dt);
    
    % 2. Schleife über jeden Zustandseintrag, um partielle Ableitungen zu bilden
    for j = 1:nx
        % Temporären gestörten Vektor erstellen
        x_perturbed = x;
        
        % Wir addieren epsilon auf den j-ten Zustand
        % Trick: Falls x(j) sehr groß ist, sollte epsilon skalieren. 
        % Hier reicht oft additives epsilon, solange x nicht riesig ist.
        delta = epsilon * max(abs(x(j)), 1); % Adaptive Skalierung (optional, aber robuster)
        x_perturbed(j) = x(j) + delta;
        
        % Berechnung des nächsten Zustands mit gestörtem Eingang
        x_next_perturbed = state_transition(x_perturbed, i_input, dt);
        
        % Differenzenquotient (Numerische Ableitung)
        % Spalte j der Matrix F
        F(:, j) = (x_next_perturbed - x_next_nominal) / delta;
    end
end

function docv = get_ocv_derivative(Ah)
    % Steigung der OCV Kurve beim aktuellen SOC
    % Wichtig für den Filter, um zu wissen, wie stark er den SOC korrigieren muss
    h = 0.05;
    docv = (voltage_to_charge(Ah+h)-voltage_to_charge(Ah-h))/(2*h);
end

function x_next = state_transition(x, i, dt)
    % Entpacken
    soc = x(1); 
    v1 = x(2); 
    v2 = x(3);
    Rs = exp(x(4)); R1 = exp(x(5)); C1 = exp(x(6)); R2 = exp(x(7)); C2 = exp(x(8));
    
    % --- Physik ---
    
    % 1. SOC Update (Coulomb Counting)
    % Vorzeichenkonvention beachten: Entladen (i>0) verringert SOC? 
    % Üblich: i>0 ist Entladen -> soc_new = soc - ...
    % Oder i>0 ist Laden -> soc_new = soc + ... 
    % Hier nehme ich an: i ist positiv beim Laden (oder Vorzeichen in u geregelt)
    soc_next = soc + dt * i;
    
    % 2. RC-Glieder (Diskretisierung exakt via Exponential)
    % Schutz gegen Division durch Null, falls C oder R vom Filter auf 0 geschätzt werden
    tau1 = R1 * C1;
    if tau1 < 1e-4, tau1 = 1e-4; end % Numerische Sicherung
    
    tau2 = R2 * C2;
    if tau2 < 1e-4, tau2 = 1e-4; end 
    
    alpha1 = exp(-dt / tau1);
    alpha2 = exp(-dt / tau2);
    
    v1_next = v1 * alpha1 + R1 * (1 - alpha1) * i;
    v2_next = v2 * alpha2 + R2 * (1 - alpha2) * i;
    
    % 3. Parameter Update (Random Walk)
    % Parameter bleiben konstant, Unsicherheit wird über Matrix Q addiert
    % rs_next = Rs;
    % r1_next = R1;
    % c1_next = C1;
    % r2_next = R2;
    % c2_next = C2;
    % 
    % Zusammenpacken
    x_next = [soc_next; v1_next; v2_next; x(4); x(5); x(6); x(7); x(8)];
end

% function [x_new, P_new] = ekf_step(x, P, i, v_meas, dt, Q, R)
% 
%     % --- A. Prediction ---
% 
%     % 1. State Prediction (nutzt jetzt die ausgelagerte Funktion)
%     x_pred = state_transition(x, i, dt);
% 
%     % 2. Jacobian F berechnen (NUMERISCH)
%     F = compute_jacobian_F(x, i, dt);
% 
%     % 3. Covariance Prediction
%     P_pred = F * P * F' + Q;
% 
%     % --- B. Correction ---
%     % ... (Rest bleibt wie zuvor, wobei Sie H ggf. auch numerisch machen können)
% 
% end

% Inline function to generate a suitably formatted URL for the endpoint:
% url_base/utt/submit/<team_number>/<pwd>/<dataset_number>/<sample_index>/<charge_estimate>
compose_url = @(dataset_number, sample_index, charge_estimate) sprintf('%s/utt/submit/%d/%s/%d/%d/%.2f', url_base, team_number, pwd, ...
    dataset_number, sample_index, charge_estimate);



% Empty current and voltage vector 
current_vector = [];
voltage_vector = [];
Rs_vector = [];
x_hat_vector = [];
resistance_estimate = 0.08;
% Charge estimate vector
charge_estimate_vector = [-1];

% Number of samples per dataset
% Dataset 1: 405 available samples.
% Dataset 2: 661 available samples.
% Dataset 3: 110 available samples.

% The first sample index is 0, the last index is what reported here above.
deltaT = 2/3600; % Sec -> Hour
Ttot = 2*405;
Ah_cnt = 0;

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

    % and append them to a vector to keep an history
    voltage_vector = vertcat(voltage_vector, data.voltage);
    current_vector = vertcat(current_vector, data.current);

    % ----
    % TODO Computation of the charge estimate. Your show, now.
    % ----

    if (sample_index>1)
        deltaV = voltage_vector(end)-voltage_vector(end-1);
        deltaI = current_vector(end)-current_vector(end-1);
        if abs(deltaI) > 1
            resistance_estimate = deltaV / deltaI;
        end
        Ah_cnt = Ah_cnt + current_vector(end)*deltaT;
        
    else
        Ah_cnt = voltage_to_charge(voltage_value);
        % Initialer Zustandsvektor x: 
        % [SOC; Vc1; Vc2; Rs; R1; C1; R2; C2]
        % WICHTIG: Gute Startwerte für Parameter raten, sonst divergiert der Filter!
        % x_hat = [Ah_cnt; 0; 0; 0.08; 0.0553; 1000; 0.0373; 3000];
        x_hat = [Ah_cnt; 0; 0; 0.0828; 0.0553; 1000; 0.0373; 5000]; 
        
        % Kovarianzmatrix P (Unsicherheit der Schätzung)
        P = eye(8) * 1e-3; 
        
        % Prozessrauschen Q (Wie sehr vertrauen wir dem Modell?)
        % Parameter haben kleines Rauschen (Random Walk), SOC folgt Strom
        diag_Q = [1e-6, 1e-5, 1e-5, 1e-8, 1e-8, 1e-6, 1e-8, 1e-6]; 
        Q = diag(diag_Q);
        
        % Messrauschen R (Wie sehr vertrauen wir dem Spannungssensor?)
        R = 0.01; % z.B. 0.01 V^2 Varianz
    end

   
    %charge = voltage_to_charge(v)
    OCV = voltage_vector(end) - resistance_estimate * current_vector(end);
    
    [x_hat, P] = ekf_step(x_hat, P, current_vector(end), voltage_vector(end), deltaT, Q, R);


    % A = [-1/(R1*C1) 0 0;0 -1/(R2*C2) 0;0 0 1];
    % B = [1/C1 0; 1/C2 0];
    % C = [1 1 0];
    % D = Rs;
    % x = [v*C1 v*c2 E];
    % u = [i; 1];
    % xdot = A*x + B*u;
    % v = C*x + D*u;

    charge_estimate = Ah_cnt;
    Rs_vector = vertcat(Rs_vector, resistance_estimate);

    x_hat_vector = vertcat(x_hat_vector, x_hat');
    % --- End of Computation of the charge estimate


    % Append the charge estimate to the vector. It will be sent out at the
    % next iteration of the for loop.
    charge_estimate_vector = vertcat(charge_estimate_vector, charge_estimate);
end
%% 

t = 0:2:Ttot-2;
% 1. Parameter schätzen
params = estimate_params(t, current_vector, voltage_vector);

% 2. EKF benutzen
dt = t(2)-t(1);
x_est = ekf_battery(current_vector, voltage_vector, dt, params);

% 3. Ergebnisse plotten
figure; 
plot(t, x_est(1,:), t, x_est(2,:), t, x_est(3,:));
legend('v_C1','v_C2', 'E');
xlabel('t'); ylabel('Spannung [V]');
title('Geschätzte RC-Spannungen (EKF)');

% 3. Ergebnisse plotten
Rs = params.Rs;


v_est = x_est(1,:) + x_est(2,:) + x_est(3,:) + current_vector'*Rs;


figure;
plot(t, voltage_vector, t, v_est');
legend('v_meas','v_est');
xlabel('t'); ylabel('Spannung [V]');
title('Geschätzte RC-Spannungen (EKF)');

%% 
x_hat_vector = [];
Ah_cnt = voltage_to_charge(voltage_vector(1));
% Initialer Zustandsvektor x: 
% [SOC; Vc1; Vc2; Rs; R1; C1; R2; C2]
% WICHTIG: Gute Startwerte für Parameter raten, sonst divergiert der Filter!
% x_hat = [Ah_cnt; 0; 0; 0.08; 0.0553; 1000; 0.0373; 3000];
x_hat = [Ah_cnt; 0; 0; log(0.0123); log(0.0678); log(1016); log(0.0547); log(3016)]; 

% Kovarianzmatrix P (Unsicherheit der Schätzung)
P = diag([
    1e-4, ... % SOC: Wir wissen Start-SOC meist auf +/- 1% genau
    1e-2, ... % Vc1: Unsicherheit beim Start (0V angenommen)
    1e-2, ... % Vc2
    0.5,  ... % ln(Rs): Erlaubt Fehlerfaktor ~2 (exp(sqrt(0.5)) ≈ 2.0)
    0.5,  ... % ln(R1)
    1.0,  ... % ln(C1): Erlaubt Fehlerfaktor ~2.7 (sehr unsicher)
    0.5,  ... % ln(R2)
    1.0   ... % ln(C2)
]);
% Prozessrauschen Q (Wie sehr vertrauen wir dem Modell?)
% Parameter haben kleines Rauschen (Random Walk), SOC folgt Strom
diag_Q = [11e-6,   ... % SOC: Folgt dem Strom, kleines Rauschen für Modellfehler
    1e-4,   ... % Vc1: Spannung darf sich dynamisch ändern
    1e-4,   ... % Vc2
    1e-8,   ... % ln(Rs): Sehr stabil. Sollte sich kaum ändern. (~0.01% step)
    1e-7,   ... % ln(R1): Etwas flexibler für Erwärmung etc.
    1e-5,   ... % ln(C1): !! Höher gesetzt !! Damit C sich bewegt (~0.3% step)
    1e-7,   ... % ln(R2)
    1e-5    ... % ln(C2): !! Höher gesetzt !!
    ]; 
Q = diag(diag_Q);

% Messrauschen R (Wie sehr vertrauen wir dem Spannungssensor?)
R = 0.01; % z.B. 0.01 V^2 Varianz

for sample_index=1:405
    [x_hat, P] = ekf_step(x_hat, P, current_vector(sample_index), voltage_vector(sample_index), deltaT, Q, R);
    x_hat_vector = vertcat(x_hat_vector, x_hat');
end
 Rs = exp(x_hat(4)); R1 = exp(x_hat(5)); C1 = exp(x_hat(6)); R2 = exp(x_hat(7)); C2 = exp(x_hat(8));
% Plot data
% Sampling time is 2 s! 
figure;
Ts = 2;
time_vector = (1:numel(voltage_vector))' * Ts;
subplot(4, 1, 1)
plot(time_vector, voltage_vector)
xlabel('Time [seconds]')
ylabel('Voltage [V]')

subplot(4, 1, 2)
plot(time_vector, current_vector)
xlabel('Time [seconds]')
ylabel('Current [V]')


subplot(4, 1, 3)
plot(time_vector, charge_estimate_vector(2:end))
xlabel('Time [seconds]')
ylabel('Charge level [Ah]')

subplot(4, 1, 4)
plot(time_vector, x_hat_vector(:,1))
xlabel('Time [seconds]')
ylabel('xhat [Ah]')


