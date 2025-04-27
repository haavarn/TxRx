% Example code showing the implementation of Eq.(10)

M = 64; % elements in array


Delta = exp(1j*pi/(2*M-1)); % phase term in Eq.(10)
w = zeros(1,M);
w(1) = 1; 
for N = 2:(M/2)
    S = 0; % initialize sum
    for ind = 2:N-1
        S = S + w(ind)*w(N+1-ind)*Delta^(N+1 - 2*ind);
    end
    w(N) = (1 - S)/(w(1)*(Delta^(N-1) + Delta^(-N+1)));
end
w(M/2+1:M) = flip(w(1:M/2)); % flip when half due to symmetry

ph = exp(1j*(0:M-1)*pi/(2*M-1)); % steering phase
w_R = w.*ph;
w_T = w.*conj(ph);


subplot(1,2,1)
plot(abs(w), 'LineWidth',2);
xlabel("Element \itm")
ylabel("|w_{Rx/Tx}|")
ylim([0, 1.1*max(abs(w))])

subplot(1,2,2)
w_EA_recombined = conv(w_R,w_T);
plot(real(w_EA_recombined),'DisplayName','Real'); hold on
plot(imag(w_EA_recombined),'DisplayName','Imag')
ylim([0, 1.1*max( real(w_EA_recombined))])
xlim([0,M])
xlabel("Element in EA")
ylabel("w_{EA} = conv(w_{Tx}, w_{Rx})")
legend