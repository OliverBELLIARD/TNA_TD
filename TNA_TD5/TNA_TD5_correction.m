clearvars; close all; clc;

Nord    = 31;
wp      = .43;
[h0,h1,g0,g1] = firpr2chfb(Nord,wp);
%fvtool(h0,1,h1,1,g0,1,g1,1);
%legend('H_0(\nu)','H_1(\nu)', 'G_0(\nu)', 'G_1(\nu)');

H0      = fft(h0, 1024);
H1      = fft(h1, 1024);
G0      = fft(g0, 1024);
G1      = fft(g1, 1024);

nu      = (0:numel(H0)-1)/numel(H0);
Zh0     = roots(h0);
Zh1     = roots(h1);
Zg0     = roots(g0);
Zg1     = roots(g1);

t       = .5*(conv(h0,g0)+conv(h1,g1));
nt      = 0:length(t)-1;
[ph0,w] = grpdelay(h0,1,512);
[ph1,w] = grpdelay(h1,1,512);

turquoise=[0,.45,.74];      % h0
orange  = [.85,.33,.10];    % h1
moutarde= [.93,.69,.13];    % g0
violet  = [.49,.18,.56];    % g1

figure;
subplot(231);
    zplane(h0,1);
    hold on;
    plot(real(Zh1),imag(Zh1),'o','color',orange); 
    hold off;
    legend('Z_{H0}','P_{H0}','unit','Z_{H1}');
subplot(234);
    zplane(g0,1);
    hold on;
    plot(real(Zg1),imag(Zg1),'o','color',violet);
    plot(real(Zg0),imag(Zg0),'o','color',moutarde);
    hold off;
    legend('Z_{G0}','P_{G0}','unit','Z_{G1}');
subplot(232);
    hold on;
    stem(0:Nord, h0,'color',turquoise);
    stem(0:Nord, h1,'color',orange);
    hold off;
    xlabel('samples'); legend('h_0[n]','h_1[n]');
subplot(235);
    hold on;
    stem(0:Nord, g0,'color',moutarde); % BUG couleur + ordre
    stem(0:Nord, g1,'color',violet);
    hold off;
    xlabel('samples'); legend('g_0[n]','g_1[n]');
subplot(233);
    plot(nu, 20*log(abs([H0;H1;G0;G1])));
    xlim([0,.5]); grid; xlabel('\nu');
    legend('|H_0(\nu)|_{dB}','|H_1(\nu)|_{dB}', '|G_0(\nu)|_{dB}', '|G_1(\nu)|_{dB}');
subplot(236);
    plot(w/(2*pi), [ph0,ph1]); %[0,diff(angle(H0))./diff(nu)/(2*pi)], nu,[0,diff(angle(H1))./diff(nu)/(2*pi)]);
    xlim([0,.5]); grid; xlabel('\nu');
    legend('\tau_{H_0}(\nu)','\tau_{H_1}(\nu)');

figure;
subplot(311);
    stem(nt, t);
    xlabel('samples');
    legend('t[n]');
subplot(312);
    plot(nu, 20*log(abs([H0.*G0;H1.*G1])));
    xlim([0,.5]); grid; xlabel('\nu');
    legend('|H_0.G_0(\nu)|_{dB}','|H_1.G_1(\nu)|_{dB}');
subplot(313);
    plot(nu, 20*log(abs(.5*(H0.*G0+H1.*G1))));
    xlim([0,.5]); grid; xlabel('\nu');
    legend('|T(\nu)|_{dB}');

%% 8-channel tree structured filter bank
M       = 2;
h02     = zeros(1,M*length(h0));    % H02(z) = H0(z^2)
h02(1:M:length(h02)) = h0;  
h12     = zeros(1,M*length(h1));    % H12(z) = H1(z^2)
h12(1:M:length(h12)) = h1;
M       = 4; 
h04     = zeros(1,M*length(h0));    % H04(z) = H0(z^4)
h04(1:M:length(h04)) = h0;  
h14     = zeros(1,M*length(h1));    % H14(z) = H1(z^4)
h14(1:M:length(h14)) = h1;

H0  = freqz(h0,1,1024)';
H1	= freqz(h1,1,1024)';
H02 = freqz(h02,1,1024)';
H12 = freqz(h12,1,1024)';
H04 = freqz(h04,1,1024)';
H14 = freqz(h14,1,1024)';
nu  = (0:1023)/2048;

Ha  = H0.*H02.*H04;
Hb  = H0.*H02.*H14;
Hc  = H0.*H12.*H04;
Hd  = H0.*H12.*H14;
He  = H1.*H02.*H04;
Hf  = H1.*H02.*H14;
Hg  = H1.*H12.*H04;
Hh  = H1.*H12.*H14;

figure;
subplot(311);
    plot(nu, abs([H0;H1]));
    legend('|H_{0}(\nu)|','|H_{1}(\nu)|');
    title('Decimation of the orthogonal 2-channel filters');
subplot(312);
    plot(nu, abs([H02;H12]));
    legend('|H_{02}(\nu)|','|H_{12}(\nu)|');
subplot(313);
    plot(nu,abs([H04;H14]));
    legend('|H_{04}(\nu)|','|H_{14}(\nu)|');
    xlabel('\nu');
    
figure;
subplot(311);
    plot(nu, abs([Ha;Hb;Hc;Hd;He;Hf;Hg;Hh]));
    legend('|H_a(\nu)|','|H_b(\nu)|','|H_c(\nu)|','|H_d(\nu)|','|H_e(\nu)|','|H_f(\nu)|','|H_g(\nu)|','|H_h(\nu)|');
    title('Tree-structured analysis filter bank');
subplot(312);
    plot(nu, 20*log(abs([Ha;Hb;Hc;Hd;He;Hf;Hg;Hh])));
    legend('|H_a(\nu)|_{dB}','|H_b(\nu)|_{dB}','|H_c(\nu)|_{dB}','|H_d(\nu)|_{dB}','|H_e(\nu)|_{dB}','|H_f(\nu)|_{dB}','|H_g(\nu)|_{dB}','|H_h(\nu)|_{dB}');
    title('Tree-structured analysis filter bank');
subplot(313);
    plot(nu,abs(ones(1,8)*[Ha;Hb;Hc;Hd;He;Hf;Hg;Hh]));
    legend('|H(\nu)|');
    xlabel('\nu');
    
%% Decomposition of a rectangle signal
x1  = [ones(1,200),zeros(1,1000)];
x2  = [ones(1,20),zeros(1,100)];
x2  = [x2,x2,x2,x2,x2,x2,x2,x2,x2,x2];

X1  = freqz(x1,1,1024)';
X2	= freqz(x2,1,1024)';

Y1  = [Ha;Hb;Hc;Hd;He;Hf;Hg;Hh].*(ones(8,1)*X1);
Y2  = [Ha;Hb;Hc;Hd;He;Hf;Hg;Hh].*(ones(8,1)*X2);

figure;
subplot(311);
    plot(nu, abs([X1;X2]));
    legend('|X_{1}(\nu)|','|X_{2}(\nu)|');
    title('Rectangular signals through the tree-structured filter bank');
subplot(312);
    plot(nu, abs(Y1));
    legend('|H_{1a}(\nu)|','|H_{1b}(\nu)|','|H_{1c}(\nu)|','|H_{1d}(\nu)|','|H_{1e}(\nu)|','|H_{1f}(\nu)|','|H_{1g}(\nu)|','|H_{1h}(\nu)|');
subplot(313);
    plot(nu,abs(Y2));
    legend('|H_{2a}(\nu)|','|H_{2b}(\nu)|','|H_{2c}(\nu)|','|H_{2d}(\nu)|','|H_{2e}(\nu)|','|H_{2f}(\nu)|','|H_{2g}(\nu)|','|H_{2h}(\nu)|');
    xlabel('\nu');

