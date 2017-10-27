clear all; close all; clc;
load('lp_fil_coeff.txt');
load('hp_fil_coeff.txt');
load('bp_fil_coeff.txt');

%%
w_in=16; %wordsize of input signal
w_h=16; %wordsize of filter coefficients
w_a=40; %wordsize of accumulator
w_out=16; %wordsize of output signal

%%----- Divide by 2 for second case--------%%
w_in=8; %wordsize of input signal
w_h=8; %wordsize of filter coefficients
w_a=20; %wordsize of accumulator
w_out=8; %wordsize of output signal

%%

b=[lp_fil_coeff,bp_fil_coeff,hp_fil_coeff];
B = 16; % Number of bits
bq = fi(b(:,1), true, B);  % signed = true, B = 16 bits
L = bq.FractionLength
bsc = b(:,1)*2^L;
hlp = fdesign.nyquist(16,bsc);
hlp.Arithmetic = 'fixed';
hlp.CoeffWordLength = 16;


bq = fi(b(:,2), true, B);  % signed = true, B = 16 bits
L = bq.FractionLength
bsc = b(:,2)*2^L;
hbp = dfilt.dffir(bsc);
hbp.Arithmetic = 'fixed';
hbp.CoeffWordLength = 16;

bq = fi(b(:,3), true, B);  % signed = true, B = 16 bits
L = bq.FractionLength
bsc = b(:,3)*2^L;
hhp = dfilt.dffir(bsc);
hhp.Arithmetic = 'fixed';
hhp.CoeffWordLength = 16;

% Check that the coefficients of h are all integers:

all(hlp.Numerator == round(hlp.Numerator))
all(hbp.Numerator == round(hbp.Numerator))
all(hhp.Numerator == round(hhp.Numerator))

fvtool(hlp, 'Color', 'white')
fvtool(hbp, 'Color', 'white')
fvtool(hhp, 'Color', 'white')

noiseIn=fi(randn(4000,1),true,w_in);
fnoiseIn=fftshift(abs(fft(noiseIn.data)));

%% Low-pass output
ylp=filter(hlp,noiseIn);
flp=fftshift(abs(fft(ylp.data)));

%% Band-pass output
ybp=filter(hbp,noiseIn);
fbp=fftshift(abs(fft(ybp.data)));

%% High-pass output
yhp=filter(hhp,noiseIn);
fhp=fftshift(abs(fft(yhp.data)));

%% Plotting combined outputs

figure
plot(flp,'b')
hold on
plot(fbp,'r')
hold on
plot(fhp,'m')
hold on
% plot(fnoiseIn,'g')

%% Check if signal is reconstructed
sum=flp+fbp+fhp;
residual=fnoiseIn-sum;
figure
plot(residual)


