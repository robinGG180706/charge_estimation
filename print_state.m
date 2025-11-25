% EPFL EE466 - Fabrizio Sossan
% 
% This script is to query and plot data on the charge estimates that you
% have delivered.
% 
% Query the endpoint 'url_base/utt/state/<team_number>/<pwd>' and print the
% results.

clc;
clear all;
close all;


% Your team number (INTEGER) and password (STRING)
team_number = 99;
pwd = '99';
% Base url (STRING)
url_base = 'http://ec2-51-20-83-40.eu-north-1.compute.amazonaws.com';

url_endpoint = sprintf('%s/utt/state/%d/%s', url_base, team_number, pwd);

try
    data_txt = webread(url_endpoint, weboptions("ContentType", "text"));
catch ME
    fprintf(sprintf("Error. The url '%s' is likely malformed.\n", url_endpoint));
    rethrow(ME)
end
try
    data = jsondecode(data_txt);
catch ME
    error(sprintf("Couldn't parse json. Error message: %s", data_txt));
end



dss = {'d1', 'd2', 'd3'};
for d=1:numel(dss)
    ds = dss{d};

    fprintf('---- Dataset %s ----\n', ds);
    fprintf('Served %d out of %d available samples.\n', ...
        data.datasets.(dss{d}).last_consumed_sample, ...
        data.datasets.(dss{d}).total_samples)


    % Plot of the served data and delivered estimations
    voltage = data.datasets.(dss{d}).observed_voltage_samples;
    current = data.datasets.(dss{d}).observed_current_samples;
    charge = data.datasets.(dss{d}).submitted_charge_estimates;


    subplot(3, 3, (d-1)*3+ 1)
    plot(voltage)
    xlabel('Samples index [-]')
    ylabel('Voltage [V]')
    

    subplot(3, 3, (d-1)*3+ 2)
    plot(current)
    xlabel('Samples index [-]')
    ylabel('Current [V]')    
    title(sprintf('Dataset %s', ds))

    subplot(3, 3, (d-1)*3+ 3)
    plot(charge(2:end))
    xlabel('Samples index [-]')
    ylabel('Charge level [Ah]')
end



