clear;
clc;
close all;

%% m 序列产生 7 个码元
BaseVal     = 2;        % 底数
PowerVal    = 3;        % 阶次（√）
N           = 10;        % 码元重复周期数
Shift       = [];       % 码元移位数目
WhichSeq    = 1;        % 本原多项式选择
Mseq        = Mseq_function(BaseVal, PowerVal, N, Shift, WhichSeq);
Code_N      = BaseVal ^ PowerVal - 1;   % 码元个数

%% 对M序列填充37.5kHz的正弦信号，每个码元填充10个周期的正弦信号
Fc          = 1 * 37.5e3;       % 信号的中心频率
Fs          = 4 * Fc;       % 信号的采样频率
Tc          = 1 / Fc;       % 信号周期
Ts          = 1 / Fs;       % 信号的重采样间隔
Period_N    = 5;           % 每个码元包含的信号周期数（√）
OneCodeSampNum      = Period_N * Fs / Fc;       % 单个码元所拥有的采样点数

%% 产生多普勒频移之后的回波信号
V       = 5.1;       % 目标速度
c       = 1500;     % 声速
Fd      = round((c + V) / (c - V) * Fc);            % 含有多普勒的频率
DeltaF  = (2*10) / (c - 10) * Fc;
PeriodMSeq  = N;                                    % M 序列的周期数
NumSampled  = OneCodeSampNum * Code_N * PeriodMSeq; % 三个周期M序列的采样点数
SignalSend  = zeros( NumSampled, 1 );               % 发送信号存储空间分配

% for nn = 1:1:Code_N * PeriodMSeq
%     t1 = (nn - 1) * OneCodeSampNum * Ts : Ts : (nn * OneCodeSampNum  - 1) * Ts; % t1从0开始的
%     SignalSend((nn - 1) * OneCodeSampNum  + 1 : 1 : nn * OneCodeSampNum ) = ...
%         sin(2 * pi * Fc * t1 + pi * (Mseq(nn)+ 1) / 2);
% end

for pp = 1:1:PeriodMSeq
    t0 = (pp - 1) * Code_N * Ts * OneCodeSampNum;
    for nn = 1:1:Code_N 
        t1 = t0 + (nn - 1) * OneCodeSampNum * Ts : Ts : t0 + (nn * OneCodeSampNum  - 1) * Ts; % t1从0开始的
        SignalSend((pp-1)*Code_N *OneCodeSampNum + (nn - 1) * OneCodeSampNum  + 1 : 1 : (pp-1)*Code_N *OneCodeSampNum + nn * OneCodeSampNum ) = ...
            sin(2 * pi * Fc * t1 + pi * (Mseq(nn)+ 1) / 2);
    end
end

figure
plot((0:NumSampled - 1) * Ts * 1000, SignalSend);
title('Siganl Transmit');
xlabel('Time/ms');
ylabel('Amplitude');
ylim([-1.3 1.3]);
grid on;

SignalReceiced      = resample(SignalSend, Fc, Fd);
NumOneCodeSamp      = length(SignalReceiced) / N / Code_N;

%% 发射信号和接收信号频谱分析
NFFT                = 2^nextpow2(length(SignalSend));
FFTSignalSend       = fft(SignalSend, NFFT);
FFTSignalReceiced   = fft(SignalReceiced, NFFT);
FrequencyDistri     = (0:1:NFFT / 2 - 1) * Fs / NFFT;

figure
plot(FrequencyDistri, abs(FFTSignalSend((1:1:NFFT / 2), 1)), '-o',...
     FrequencyDistri, abs(FFTSignalReceiced((1:1:NFFT / 2 ), 1)), '-*');
title('发送信号和接收信号的频谱');
xlabel('Frequency /Hz');
ylabel('Amplitude');
grid on;

%% 带通滤波设计，添加噪声用
B = 1 / (1 / Fd * Period_N);    % 编码信号的带宽

rp_passband  = 0.14;            % Passband ripple
rs_passband  = 75;              % Stopband ripple
fs_passband  = 600;             % Sampling frequency

% FreqPassNorm1   = (Fs - B ) /(2 * Fs);
% FreqPassNorm2   = 1 - FreqPassNorm1;
% FreqPass1       = FreqPassNorm1 * fs_passband / 2;
% FreqPass2       = FreqPassNorm2 * fs_passband / 2;
% FreqStop1       = FreqPass1 - 10;
% FreqStop2       = FreqPass2 + 10;

FreqPassNorm1   = 2 *(Fc - (B + DeltaF) / 2 ) /(Fs);
FreqPassNorm2   = 2 *(Fc + (B + DeltaF) / 2 ) /(Fs);
FreqPass1       = FreqPassNorm1 * fs_passband / 2;
FreqPass2       = FreqPassNorm2 * fs_passband / 2;
FreqStop1       = FreqPass1 - 10;
FreqStop2       = FreqPass2 + 10;

f_passband   = [FreqStop1 FreqPass1 FreqPass2 FreqStop2];     % Cutoff frequencies
a_passband   = [0 1 0];         % Desired amplitudes

dev_passband        = [10^(-rs_passband/20) (10^(rp_passband/20)-1)/(10^(rp_passband/20)+1)  10^(-rs_passband/20)];
[n,fo,ao,w]         = firpmord(f_passband,a_passband,dev_passband,fs_passband);
b_passband          = firpm(n,fo,ao,w);         % 求得滤波器系数

figure
freqz(b_passband,1,1024,fs_passband)

%% 给回波信号加噪声
SNR     = 10;           % 信噪比
Noise   = randn(2 * length(SignalReceiced),1);          % 产生长度为信号长度两倍的噪声
NoiseAfterFilter = filter(b_passband.', 1, Noise);      % 添加带限噪声
NoiseCut        = ...
    NoiseAfterFilter(length(Noise)/2 - length(SignalReceiced)/2 : length(Noise)/2 + length(SignalReceiced)/2 - 1);
% 选取噪声中间和信号点数相同的部分
EnergyNoise     = NoiseCut' * NoiseCut;                 % 噪声能量
EnergySignal    = SignalReceiced' * SignalReceiced;     % 信号的能量
NoiseNorm       = NoiseCut / sqrt(EnergyNoise);         % 噪声归一化
SignalNorm      = SignalReceiced / sqrt(EnergySignal);  % 信号归一化
CoeffSnr        = 10^(-SNR/20);
NoiseSnr        = CoeffSnr * NoiseNorm * sqrt(EnergyNoise);
SignalAddNoise  = SignalReceiced + NoiseSnr;

%% 信噪比验证
SnrVerify = 10 * log10( EnergySignal / (NoiseSnr' * NoiseSnr));


%% 低通滤波器设计
Rp_LowPass  = 0.05;         % Passband ripple
Rs_LowPass  = 55;          % Stopband ripple
Fs_LowPass  = 600;         % Sampling frequency

FreqPass3 = 2 * (B + DeltaF) / Fs * Fs_LowPass / 2;
FreqStop3 = FreqPass3 + 10;
F_LowPass   = [FreqPass3 FreqStop3];     % Cutoff frequencies
% F_LowPass   = [15 20];     % Cutoff frequencies
A_LowPass   = [1 0];       % Desired amplitudes

Dev_LowPass     = [(10^(Rp_LowPass/20)-1)/(10^(Rp_LowPass/20)+1),...
                    10^(-Rs_LowPass/20)];
[n,fo,ao,w]     = firpmord(F_LowPass,A_LowPass,Dev_LowPass,Fs_LowPass);
B_Lowband       = firpm(n,fo,ao,w);         % 求得滤波器系数

figure
freqz(B_Lowband,1,1024,Fs_LowPass)

%% 基带解调
t3              = (0:1:length(SignalReceiced)-1).'*Ts;
BaseSignal      = exp( -1i * 2 * pi * Fc * t3);
SignalDemodu    = BaseSignal .* SignalAddNoise;
SignalFiltered  = filter(B_Lowband, 1, SignalDemodu);

fftSignalFiltered = fft(SignalFiltered,Fs);
fftSignalDemodu   = fft(SignalDemodu,Fs);
figure
plot(abs(fftSignalDemodu))
hold on
plot(abs(fftSignalFiltered))
%% 复相关法测频
% 信号延迟8个码元

Delay       = round(NumOneCodeSamp * Code_N);
Segment1    = SignalFiltered(1:1:1 + length(SignalFiltered)/2);
Segment2    = SignalFiltered(Delay + 1 : Delay + 1 + length(SignalFiltered)/2);
autocorr    = sum(Segment1.* conj(Segment2));
theta       = -atan2(imag(autocorr), real(autocorr));
EstiFreq    = theta / (2 * pi * Delay * Ts);
EstiV       = EstiFreq * c / (2 * Fc);

fprintf('估计目标速度：%.3f\n', EstiV);
