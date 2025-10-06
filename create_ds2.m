clc; clear all; close all;

shot = 87899;
% shot = 79742;

t = liuqe_times(shot); % get times

% decimate times
t = t(1:20:end);

% %get the second half
% t = t(t>0.25);

disp(['Number of times: ', num2str(length(t))]);

[L, LX] = liuqe(shot,t,'LYall',false,'ndom',5);

LY = liut(L, LX);


% liuqeplot({{L,LX,LY}},{'MEQ'});


% plot Fx
figure; hold on; box on; axis equal;
t_idx = 0.25; % time index to plot
[~, idx] = min(abs(t - t_idx));
[rr, zz] = meshgrid(L.rx, L.zx);
Fx = LY.Fx(:,:,idx);
scatter(rr(:), zz(:), 20, Fx(:), 'filled');
xlabel('r [m]'); ylabel('Fx [N/m]');
title(['Fx at t = ', num2str(t(idx)), ' s']);
grid on;    


