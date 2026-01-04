function PV = time_days_to_PV(t_days, pvph)
% time_days_to_PV converts time (days) to pore volumes using PV per hour (pvph)
PV = t_days * 24 * pvph;
end
