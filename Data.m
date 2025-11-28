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
charge_estimate_vector = [-1];
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

    % and append them to a vector to keep an history
    voltage_vector = vertcat(voltage_vector, data.voltage);
    current_vector = vertcat(current_vector, data.current);
end
% a = linspace(4.35, 4.4, 35)  + 0.05*rand(1,35);
% a = a';
% b = linspace(2.85,2.55,115) + +0.05*rand(1,115);
% b = b';
% 
% c = 7.4 * ones(1, 35) + 0.1*rand(1,35);
% c = c';
% d = -12.5*ones(1, 115) + 0.1*rand(1,115);
% d = d';
% 
% voltage_vector = [a; b];
% current_vector = [c; d];