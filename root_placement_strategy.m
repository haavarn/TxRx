M = 10;
w_EA = chebwin(2*M-1, 30);
w_EA = w_EA/sum(w_EA); % normalize to unit gain: sum(w)= 1
tic;
R = roots(w_EA); % complex roots of effective aperture

% Alternating root strategy for -30 dB Chebyshev window
% NB: if roots are not on the unit circle, consider pairing if reciprocal 

% Sort roots by angle, and distribute every other to Tx and Rx 
[~,sidx] = sort(wrapTo2Pi(angle(R))); 
Rsort = R(sidx); 
R_rx = Rsort(1:2:end);
R_tx = Rsort(2:2:end);

% Create normalized Tx and Rx windows 
w_rx = poly( R_rx ); 
w_rx = w_rx/ abs(sum(w_rx)); 
w_tx = poly( R_tx ); 
w_tx = w_tx/ abs(sum(w_tx));
toc;
% Plot roots, beampattern, and window
set(groot,'defaultAxesFontSize',16)
fig = figure(1);
fig.Position = [198 299 1035 335];

subplot(1,3,1)
plot( cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k:', 'linewidth', 1 ); hold on
plot( real(R), imag(R), 'ko', 'markersize', 10)
plot( real(R_rx), imag(R_rx), 'r^', 'markerfacecolor', 'r')
plot( real(R_tx), imag(R_tx), 'bv', 'markerfacecolor', 'b')
axis equal
ylabel("Imag")
xlabel("Real")
plot([0 0], [-2 2], 'k-','linewidth',1,"handlevisibility", 'off')
plot([-2 2], [0 0], 'k-','linewidth',1,"handlevisibility", 'off')
ylim([-1.5 1.5])
title('Roots')

subplot(1,3,2)
one_way_spectrum(w_rx, 'r-', 'linewidth',3,'displayname','Rx'); hold on
one_way_spectrum(w_tx, 'b-', 'linewidth',3,'displayname','Tx');
one_way_spectrum(w_EA, 'k-', 'linewidth',2,'displayname','EA');
ylim([-60, 10])
ylabel("Power (dB)")
xlabel("\itu\rm = sin(\theta)")
legend()
title('Beampatterns')

subplot(1,3,3)
plot(1:M, abs(w_rx), 'r^-', 'displayname', '$w_{rx}$', 'LineWidth',4, 'displayname', 'Rx ampl.', 'MarkerFaceColor','r'); hold on
plot(1:M, abs(w_tx), 'bv-', 'displayname', '$w_{tx}$', 'LineWidth',2, 'displayname', 'Tx ampl.', 'MarkerFaceColor','b');
plot(1:M, (unwrap(angle(w_rx))/(pi))*max(abs(w_tx))+0.5*max(abs(w_tx)),'r^--', 'LineWidth',3, 'displayname', '[Rx phase]', 'MarkerFaceColor','r'); hold on
plot(1:M, (unwrap(angle(w_tx))/(pi))*max(abs(w_tx))+0.5*max(abs(w_tx)),'bv--', 'LineWidth',2, 'displayname', '[Tx phase]', 'MarkerFaceColor','b');
ax = gca;
set(ax,'YTick',[0,max(abs(w_tx))])
set(ax,'YTickLabel',{"0 [-\pi/2]","0.13 [\pi/2]","0.26 [\pi/2]"})
ytickangle(65)
ylim([0, max(abs(w_tx))])
legend('location', 'southwest', 'NumColumns',2);
xlabel('Element \itm')
xlim([1,M])
ylabel('|w| [\angle \itw\rm]')
title('Windows')


%% Calculate WNGP of all root combinations

idx_Alternatives = nchoosek([1:2*M-2], M-1); % all combinations of M-1 roots from 2M-2 alternatives
AllIdx = 1:2*M-2; % all possible indices
nAlts = size(idx_Alternatives, 1); % number of alternative combinations to consider

WNGP = zeros(nAlts,1); 
for idx = 1:nAlts
    w_rx = poly( R(idx_Alternatives(idx,:)) );
    w_tx = poly( R( ~ismember(AllIdx, idx_Alternatives(idx,:)) ));
    
    w_rx = w_rx/ abs(sum(w_rx));
    w_tx = w_tx/ abs(sum(w_tx));
        
    WNGP(idx) = WNG(w_tx)*WNG(w_rx);
end 

[m,mind] = max(WNGP);

R_rx_found = R(idx_Alternatives(mind,:));
R_tx_found = R( ~ismember(AllIdx, idx_Alternatives(mind,:)) );

WNGP_of_selected = WNG( poly(R_rx_found)) * WNG( poly(R_tx_found));

% Histogram
fig = figure(1);
fig.Position = [560 177 1000 250];
set(gcf,'color','w');
set(fig,'DefaultLineLineWidth',3)

subplot(1,2,1)
histogram(WNGP)
set(gca,'YScale','log')
xlabel('WNGP')
ylabel('Count')
yticks([1 10 100 1000 10000])
ylim([1 50000])
ax = gca;
ax.FontSize = 15;


subplot(1,2,2)
plot( cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k:', 'linewidth', 1 ); hold on
plot( real(R), imag(R), 'ko', 'markersize', 10)
plot( real(R_rx_found), imag(R_rx_found), 'r^', 'markerfacecolor', 'r')
plot( real(R_tx_found), imag(R_tx_found), 'bv', 'markerfacecolor', 'b')
axis equal
ylabel("Imag")
xlabel("Real")
plot([0 0], [-2 2], 'k-','linewidth',1,"handlevisibility", 'off')
plot([-2 2], [0 0], 'k-','linewidth',1,"handlevisibility", 'off')
ylim([-1.5 1.5])
title(strcat('Roots, WNGP = ',num2str(WNGP_of_selected,4)))





%% Helper function

function [Y, x] = one_way_spectrum(w, varargin)
    Y = fft(w, 2048);
    Y = abs(fftshift(Y));
    Y = db(Y);
    x = linspace(-1,1, 2048);
   
    if nargout < 2
        plot(x,Y, varargin{:});
        ylim([-60,0]);
    end
end

function output = WNG(weights)
    w = weights(:);
    output = abs(sum(w)).^2 / sum(abs(w).^2);
end

