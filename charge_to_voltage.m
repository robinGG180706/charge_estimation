function voltage_query = charge_to_voltage(chg_query)

voltage = [2.8500    2.9500    3.0500    3.1500    3.2500    3.3500    3.4500    3.5500    3.6500    3.7500];
charge = [0.0004    0.0209    0.0509    0.0911    0.1584    0.2445    0.3492    0.7328    1.1148    1.4882];
voltage_query = interp1(charge(:), voltage(:), chg_query(:), 'linear', 'extrap');

end