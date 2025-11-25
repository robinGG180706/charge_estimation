clc;
clear all;
close all;


v = 2.85:0.1:4;
charge = voltage_to_charge(v);

plot(charge, v)
xlabel('Charge [Ah]')
ylabel('Voltage [V]')