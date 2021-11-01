function Plot_1Df_v2(freq,freq_no_coh,Spectrum_Energy,freq_resolution, ...
    mfreq,Coherent_Percentage,idx_theta,mtheta,host_freqs_idx)

%% 1d spectrum, m2/Hz
k = 0;
for i = 1:mfreq
    if (i>1)&&(freq(i) == freq(i-1))
        reps_condensed(k) = reps_condensed(k)+1;
    else
        k = k+1;
        reps_condensed(k) = 1;
    end
end

%% number of wave components in an arrow
arrow_size = host_freqs_idx(1);
for i = 2:length(host_freqs_idx)
    arrow_size(i) = host_freqs_idx(i) - host_freqs_idx(i-1);
end

Syy1d_f = zeros(1,length(host_freqs_idx));
k = 0;
for i = 1:length(host_freqs_idx)
    for j = 1:arrow_size(i)
        k = k+1;
        Syy1d_f(i) = Syy1d_f(i)+Spectrum_Energy(k)*pi;
    end
    Syy1d_f(i) = Syy1d_f(i)/arrow_size(i);
end

Hmo_Syy1d_f = 4*sqrt(sum(Syy1d_f.*freq_resolution));
Hmo_Syy1d_f = 0.01*round(Hmo_Syy1d_f*100);

figure(11)
plot(Syy1d_f,freq_no_coh(host_freqs_idx),'linewidth',1.3), hold on
set(gca, 'XDir','reverse')

ylim([0.03 0.26]); 
xlim([0 0.3]);

ylabel('$f (Hz)$','Interpreter','latex')
xlabel('$E (m^2/Hz)$','Interpreter','latex')

% title(['Hmo = ' num2str(Hmo_Syy1d_f) ' m']);

grid on
grid minor

% xticks([0 0.2])

pbaspect([1 2 1])
set(gca,'fontsize', 18)

saveas(gcf,fullfile(['1Df_' num2str(Coherent_Percentage) 'p_vertical.png']))

close