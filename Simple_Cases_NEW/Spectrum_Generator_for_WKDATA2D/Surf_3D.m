function Surf_3D(theta,freq,Spectrum_Energy,theta_mean,mtheta)

% theta = linspace(-pi/2,pi/2,mtheta);
% Spectrum_Energy = reshape(Spectrum_Energy,[mtheta,length(freq)/mtheta]);
% freq = reshape(freq,[mtheta,length(freq)/mtheta]);
% freq = freq(1,:);
figure(51)
scatter3(theta,freq,Spectrum_Energy')
xlim([-pi/2+theta_mean*pi/180-0.2 pi/2+theta_mean*pi/180+0.2])
ylim([0.03 0.26])
set(gca,'fontsize', 18)
grid on
h = colorbar;
set(gca,'Box','on');
% caxis([0 0.2])
% caxis([0 2])
%% line on fp and theta mean
xline(theta_mean*pi/180,'-.k','linewidth',1)
yline(fp,'-.k','linewidth',1)

xlabel('$\theta (rad)$','Interpreter','latex')
ylabel('$f (Hz)$','Interpreter','latex')
ylabel(h,'$E (m^2/Hz/rad)$','Interpreter','latex')
saveas(gcf,fullfile(['2D_' num2str(Coherent_Percentage) 'p.png']))


% large view
ylim([0.066 0.086])
pbaspect([3 1 1])
saveas(gcf,fullfile(['2D_' num2str(Coherent_Percentage) 'p_detail.png']))