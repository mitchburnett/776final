clearvars;
figC =1;     % figure counter
alpha = .5;  % excess bandwith
Tfsamp = 10; % samples/symbol
SPAN = 12;   % symbols spanned by g(t)

SL = -SPAN/2;
SR = SPAN/2;

% pulse shape
g = rcosdesign(alpha, SPAN, Tfsamp).';
gn = g(1:Tfsamp:end);
g_span = linspace(SL, SR, length(g));
gn_span = linspace(SL, SR, length(gn));

% channel model
problem = 3;
load_c = 'measured_channel/c7.mat';
if problem == 0
    c = zeros(11,1);  
    c(1) = 1;
    c(11) = -.9;
    cn = c(1:Tfsamp:end);
elseif problem == 1
    c = zeros(31, 1);
    c(1)  = 1;
    c(6)  = -.9;
    c(25) = .5;
    cn = c(1:Tfsamp:end);
elseif problem == 2
    c = zeros(31, 1);
    c(1)  = 1;
    c(4)  = -1/2;
    c(10) = 1/4;
    c(22) = -1/8;
    cn = c(1:Tfsamp:end);
elseif problem ==3
    load(load_c);
    % extend c to span integer multiple symbols.
    sym_diff = Tfsamp-mod(length(c), Tfsamp);
    c = [c; zeros(sym_diff+1,1)];
    cn = c(1:Tfsamp:end);
end

%% Channel Model
% h(t) = g(t)*c(t);
h = conv(g,c);
hn = h(1:Tfsamp:end);
h_span = linspace(SL,(SR+length(cn)-1),length(h));
hn_span = linspace(SL, (SR+length(cn)-1), length(hn));

%% Ungerboeck Observation Model
% x(t) = h(t)*h(-t) and downsampled xn = x(nT)
x = conv(h,conj(flip(h)));
xn = x(1:Tfsamp:end);
XL = -(length(xn)-1)/2;
XR = -XL;
x_span = linspace(XL,XR,length(x));
xn_span = linspace(XL,XR, length(xn));

% 'prune' xn to 99.99% energy
k = 0;
max_k = (length(xn)+1)/2;
dxn = xn(2)-xn(1); % this isnt right.... isnt dxn 1 because they are space 1 apart?
E = 1;
Etot = sum(abs(xn).^2)*dxn;
while E > .9999
    k = k+1;
    if k > max_k
        E = 0;
    else
        tmp = sum(abs(xn(k:end-(k-1))).^2)*dxn;
        E = tmp/Etot;
    end
end

if E == 0
    disp('No min xn possible');
else
    k = k-1;
    xn_short = xn(k:end-(k-1));
end

xn_short_span = -(max_k-k):(max_k-k);

%% Forney Observation Model
% perform spectral factorization
xn_roots = roots(xn_short).';
f_n_plus = -xn_roots(abs(xn_roots)<1);

% number of polynomials
k_roots = sum(abs(xn_roots)<1);
f_n = [1, f_n_plus(1)];
k = 1;
while k < k_roots
    k = k+1;
    f_n = conv(f_n, [1, f_n_plus(k)]);
end

% normalize f_n to unit energy.
sigma_sq = sum(abs(f_n).^2);
sigma = sqrt(sigma_sq);

f_n = f_n/sigma;

%% Pulse Shape Matched Filter Observation Model
p = conv(h, conj(flip(g)));
%p = [zeros(15,1); p; zeros(15,1)]; % pad on the end with zeros to get a signal that sample on correct intger sample periods
pn = p(6:Tfsamp:end);
% PL = -(length(pn)-1)/2;
% PR = -PL;
% 
% p_span = linspace(PL, PR, length(p));
% pn_span = linspace(PL, PR, length(pn));

%% 'prune' pn to 99.99% energy.
max_kp = (length(pn)+1)/2;
kp_left = 1;
kp_right = length(pn);

[Y, idx] = sort(pn, 'descend');

pn_short = pn(sort(idx(1:6)));

%normalize
sigma_pn_sq = 1/sum(abs(pn).^2);
sigma_pn = sqrt(sigma_pn_sq);
pn_short = sigma_pn*pn_short;

pn_short_span = linspace(-2, 3, length(pn_short));

%%
% E = [1, 1]; % left and right
% Etot_p = sum(abs(pn).^2);
% while E(1) > .9999 && E(2) > .9999
%     
%    %check left
%    if kp_left +1 > max_kp
%        E_left = 0;
%    else
%        tmp = sum(abs(pn(kp_left+1:kp_right)).^2); % advance left here as to not have to reset it.
%        E_left = tmp/Etot_p;
%    end
%    
%    %check right
%    if kp_right - 1 < max_kp
%        E_right = 0;
%    else
%        tmp = sum(abs(pn(kp_left:kp_right-1)).^2);
%        E_right = tmp/Etot_p;
%    end
%    
%    % advance
%    kp_left = kp_left + 1;
%    kp_right = kp_right -1;
%    
%    %save
%    E = [E_left, E_right];
% end
% 
% if E(2) > E(1)
%     kp_right = kp_right+1;
% else
%     kp_left = kp_left-1;
% end
% 
% pn_short = pn(kp_left:kp_right);
% pn_short_span = -(max_kp-kp_left):(kp_right-max_kp);
%%
% compute frequency domain for each signal.
nfft = 1024;
G  = fftshift(fft(g,nfft));
C  = fftshift(fft(c,nfft));
H  = fftshift(fft(h,nfft));
X  = fftshift(fft(x,nfft));
XN = fftshift(fft(xn_short,nfft));
F  = fftshift(fft(f_n, nfft));
P = fftshift(fft(p,nfft));
PN = fftshift(fft(pn_short,nfft));

f_min = -Tfsamp/2;
f_max = Tfsamp/2;
freq_span = linspace(f_min,f_max,nfft);
dtft_f_span = linspace(-.5, .5, nfft);

%% Save models
Forney.L1 = 0;
Forney.L2 = length(f_n)-1;
Forney.filter = f_n.';

% Ungerboeck with ways to make correlated noise.
newx = conv(f_n, conj(flip(f_n)));
Ungerboeck.L1 = -(length(xn_short)-1)/2;
Ungerboeck.L2 = (length(xn_short)-1)/2;
Ungerboeck.old_filter = xn_short;

Ungerboeck.filter = newx.';
Ungerboeck.noise_filter = f_n.';

PSMFOM.L1 = min(pn_short_span);
PSMFOM.L2 = max(pn_short_span);
PSMFOM.filter = pn_short;

SAVE = 1;
if SAVE
    save('Models/Forney.mat', 'Forney');
    save('Models/Ungerboeck.mat', 'Ungerboeck');
    save('Models/PSMFOM.mat', 'PSMFOM');
end

%% Create plots
PLOT = 0;
if PLOT
    width = 560; height = 160;
    xPos = 470; yPos = 980;

    %% Plot g(t), g(nT) Measured channel, h(t) = g(t)*c(t)
    hFig = figure(figC); clf; figC=figC+1;
    subplot(311);
    plot(g_span,g); grid on; hold on;
    stem(gn_span,gn); hold off;
    ylabel('Pulse Shape');
    xlabel('t/T (symbols)');
    legend('g(t)', 'g(nT)');

    subplot(312);
    stem((0:length(c)-1)/Tfsamp, abs(c)); grid on;
    ylabel({'Measured Channel', 'Magnitude'});
    xlabel('t/T (symbols)');

    subplot(313);
    plot(h_span, abs(h)); grid on; hold on;
    stem(hn_span, abs(hn)); hold off;
    ylabel('Channel Model');
    xlabel('t/T (symbols)');
    legend('h(t)', 'h(nT)');

    %% plot G(f) and C(f)
    hFig = figure(figC); clf; figC=figC+1;
    set(hFig, 'Position', [xPos, yPos, width, height]);
    plot(freq_span, 20*log10(abs(G))); grid on; hold on;
    plot(freq_span, 20*log10(abs(C))); grid on; hold off;
    ylim([-40, 20]);
    legend('G(f)', 'C(f)');
    ylabel('Magnitude (dB)');

    %% plot G(f) and H(f)
    hFig = figure(figC); clf; figC=figC+1;
    set(hFig, 'Position', [xPos, yPos, width, height]);
    plot(freq_span, 20*log10(abs(G))); grid on; hold on;
    plot(freq_span, 20*log10(abs(H))); grid on; hold off;
    ylim([-40, 20]);
    xlim([-5, 5]);
    legend('G(f)', 'H(f)');
    ylabel('Magnitude (dB)');

    %% plot G(f) and X(f)
    hFig = figure(figC); clf; figC=figC+1;
    set(hFig, 'Position', [xPos, yPos, width, height]);
    plot(freq_span, 20*log10(abs(G))); grid on; hold on;
    plot(freq_span, 20*log10(abs(X))); grid on; hold off;
    ylim([-20, 40]);
    xlim([-5, 5]);
    legend('G(f)', 'X(f)');
    ylabel('Magnitude (dB)');

    %% zplane for spectral factorization
    N = 40;
    drad = 0:2*pi/N:2*pi-2*pi/N;
    unit_circle = exp(-1j*drad);

    hFig = figure(figC); clf; figC=figC+1;
    set(hFig, 'Position', [xPos, yPos, 460, 460]);
    plot(real(unit_circle), imag(unit_circle), '--'); grid on; hold on;
    plot(real(xn_roots), imag(xn_roots), 'o'); hold off;
    LIM = ceil(max(abs(xn_roots)));
    xlim([-LIM, LIM]); ylim([-LIM, LIM]);
    axis('square');
    title('Spectral Factorization of X(z)');

    %% Plot x(t), and p(t) with their downsampled sequence
    hFig = figure(figC); clf; figC=figC+1;
    subplot(211);
    plot(x_span,abs(x)); grid on; hold on;
    stem(xn_span, abs(xn)); hold off;
    ylabel('Magnitude');
    legend('x(t)', 'x(nT)');

    subplot(212);
    plot(p_span, abs(p)); grid on; hold on;
    stem(pn_span, abs(pn));
    ylabel('Magnitude');
    legend('p(t)', 'p(nT)');

    %% Plot the 'pruned' sequences for the observation models
    hFig = figure(figC); clf; figC=figC+1;
    subplot(311);
    stem(xn_short_span, abs(xn_short)); grid on;
    title('Observation Model Sequences');
    ylabel({'UOM','Magnitude'});
    xlabel('Sample Index');
    ylim([0,1]);

    subplot(312);
    stem((0:length(f_n)-1), abs(f_n)); grid on; % get xaxis for fn
    ylabel({'FOM','Magnitude'});
    xlabel('Sample Index');
    ylim([0,1]);

    subplot(313);
    stem(pn_short_span, abs(pn_short)); grid on;
    ylabel({'PSMFOM','Magnitude'});
    xlabel('Sample Index');
    ylim([0,1]);

    %% Plot the frequency response of the observation models
    hFig = figure(figC); clf; figC=figC+1;
    subplot(311);
    plot(dtft_f_span, 20*log10(abs(PN))); grid on;
    ylim([-60, 20]);
    title('Frequency Response Comparison of Observation Model');
    xlabel('Frequency (cycles/symbol)');
    ylabel({'UOM', 'Magnitude (dB)'});

    subplot(312);
    plot(dtft_f_span, 20*log10(abs(F))); grid on;
    ylim([-60, 20]);
    ylabel({'FOM', 'Magnitude (dB)'});

    subplot(313);
    plot(dtft_f_span, 20*log10(abs(XN))); grid on;
    ylim([-60, 20]);
    xlabel('frequency (cycles/symbol)');
    ylabel({'PSMFOM', 'Magnitude (dB)'});

    %% print plots
    printPlots = 0;
    PATH = '/Users/mitcburnett/Dropbox/BYU/776/Homework/hw12/img_other/';
    if printPlots
        for k_plot = 1:(figC-1)
            fig = figure(k_plot);
            width = 560; height = 160;
            xPos = 470; yPos = 980;
            set(hFig, 'Position', [xPos, yPos, 460, 460]);
            fig.PaperPositionMode = 'auto';
            %print(figure(k_plot), [PATH, 'fig_', num2str(k_plot), '.eps'], '-depsc', '-r0');
            print(fig, [PATH, 'fig_', num2str(k_plot), '_1', '.eps'],'-depsc','-r0');
        end
    end
    
end % if PLOT


