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
% w_in=8; %wordsize of input signal
% w_h=8; %wordsize of filter coefficients
% w_a=20; %wordsize of accumulator
% w_out=8; %wordsize of output signal

%%

b=[lp_fil_coeff,bp_fil_coeff,hp_fil_coeff];
B = 16; % Number of bits
bq = fi(b(:,1), true, B);  % signed = true, B = 16 bits
L = bq.FractionLength
bsc = b(:,1)*2^L;
hlp = dfilt.dffir(bsc);
hlp.Arithmetic = 'fixed';
hlp.CoeffWordLength = 16;
hlp.FilterInternals = 'SpecifyPrecision';
hlp.AccumWordLength = w_a;


bq = fi(b(:,2), true, B);  % signed = true, B = 16 bits
L = bq.FractionLength
bsc = b(:,2)*2^L;
hbp = dfilt.dffir(bsc);
hbp.Arithmetic = 'fixed';
hbp.CoeffWordLength = 16;
hbp.FilterInternals = 'SpecifyPrecision';
hbp.AccumWordLength = w_a;

bq = fi(b(:,3), true, B);  % signed = true, B = 16 bits
L = bq.FractionLength
bsc = b(:,3)*2^L;
hhp = dfilt.dffir(bsc);
hhp.Arithmetic = 'fixed';
hhp.CoeffWordLength = 16;
hhp.FilterInternals = 'SpecifyPrecision';
hhp.AccumWordLength = w_a;

% Check that the coefficients of h are all integers:

all(hlp.Numerator == round(hlp.Numerator))
all(hbp.Numerator == round(hbp.Numerator))
all(hhp.Numerator == round(hhp.Numerator))

% fvtool(hlp, 'Color', 'white')
% fvtool(hbp, 'Color', 'white')
% fvtool(hhp, 'Color', 'white')

noiseIn=fi(-1+(1-(-1)).*rand(4000,1),true,w_in,15);% Following formula a+(b-a).*rand() for [a,b] interval
noiseL=noiseIn.FractionLength
noiseInSC=noiseIn*2^noiseL;
fnoiseIn=fftshift(abs(fft(noiseInSC.data)));

%% Low-pass output
ylp=conv(hlp.Numerator,noiseInSC);
flp=fftshift(abs(fft(ylp.data)));

%% Band-pass output
ybp=conv(hbp.Numerator,noiseInSC);
fbp=fftshift(abs(fft(ybp.data)));

%% High-pass output
yhp=conv(hhp.Numerator,noiseInSC);
fhp=fftshift(abs(fft(yhp.data)));

%% Truncating outputs

% ylp=fi(ylp, true, 16, -15);
% ybp=fi(ybp, true, 16, -17);
% yhp=fi(yhp, true, 16, -16);

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
residual=sum-fnoiseIn;
figure
plot(residual)


