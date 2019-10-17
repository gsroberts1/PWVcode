function [pwv] = pwv_from_pressure(upstream, downstream, sampling_rate)

upstream = upstream';
downstream = downstream';

% Shift curves to baseline
upstream = upstream - min(upstream);
downstream = downstream - min(downstream);

figure;plot(upstream); hold on; plot(downstream)

interp_kernel = 10;

[upstream_peak, upstream_I] = max(upstream); clear upstream_peak;
[downstream_peak, downstream_I] = max(downstream); clear downstream_peak;

time = 0:sampling_rate:sampling_rate*(size(upstream,2)-1);

[time,upstream] = lowpassfilter(upstream,time);
[time,downstream] = lowpassfilter(downstream,time);

% time_interpolant = time(1):time(end)/200:time(end);
% figure;plot(time, upstream);hold on;
% upstream = spline(time, upstream, time_interpolant);
% plot(time_interpolant, upstream);
% downstream = spline(time, downstream, time_interpolant);



% counter = 1;
% for i = upstream_I:-1:interp_kernel+1
%     upstream_slope(counter) = (upstream(i)-upstream(i-interp_kernel))/(time(i)-time(i-interp_kernel));
%     counter = counter+1;
% end
% 
% [slope_max, slope_max_I] = max(upstream_slope); clear slope_max;
% slope_max_I = upstream_I - slope_max_I;
% upstream_points = upstream(slope_max_I-12:slope_max_I+12);
% upstream_points_time = time(slope_max_I-12:slope_max_I+12);
% p_upstream = polyfit(upstream_points_time,upstream_points,1);
% figure;plot(time, upstream); hold on;
% f = @(x) p_upstream(1)*x + p_upstream(2); % Define an anonymous function, y = mx + b
% fplot(f,[0 time(upstream_I)])
% 
% counter = 1;
% for i = downstream_I:-1:interp_kernel+1
%     downstream_slope(counter) = (downstream(i)-downstream(i-interp_kernel))/(time(i)-time(i-interp_kernel));
%     counter = counter+1;
% end
% 
% [slope_max, slope_max_I] = max(downstream_slope); clear slope_max;
% slope_max_I = upstream_I - slope_max_I;
% downstream_points = downstream(slope_max_I-12:slope_max_I+12);
% downstream_points_time = time(slope_max_I-12:slope_max_I+12);
% p_downstream = polyfit(downstream_points_time,downstream_points,1);
% figure;plot(time, downstream); hold on;
% f = @(x) p_downstream(1)*x + p_downstream(2); % Define an anonymous function, y = mx + b
% fplot(f,[0 time(downstream_I)])



% Hard code
upstream_x = time(79:99);
upstream_y = upstream(79:99);
p_upstream = polyfit(upstream_x,upstream_y,1);
figure;plot(time, upstream); hold on;
f = @(x) p_upstream(1)*x + p_upstream(2); % Define an anonymous function, y = mx + b
fplot(f,[0 time(end)])
ttf_upstream = -1*p_upstream(2)/p_upstream(1);

downstream_x = time(105:123);
downstream_y = downstream(105:123);
p_downstream = polyfit(downstream_x,downstream_y,1);
figure;plot(time, downstream); hold on;
f = @(x) p_downstream(1)*x + p_downstream(2); % Define an anonymous function, y = mx + b
fplot(f,[0 time(end)])
ttf_downstream = -1*p_downstream(2)/p_downstream(1);

pwv = 0.1/(ttf_downstream-ttf_upstream);






