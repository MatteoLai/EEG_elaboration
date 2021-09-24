close all
clear 
clc

%%
patient = menu('Choose the patient:',1,2);
switch patient 
    case 1
    load Patient1.mat
    case 2
    load Patient2.mat
end

Fs = 128; % [Hz]

%% Visualize EEG and ECG signals:
[n,m] = size(X_EEG);
t = (0:m-1)/Fs;

X_EEG = X_EEG - mean(X_EEG,2); % centro i dati

choice = menu('What do you want to visualize?','The whole signal','Only 30 seconds');
switch choice 
    case 1
    tstart = 1;
    tend = m;
    case 2
    t1 = input('First second of visualization: ');
    tstart = t1*Fs;
    tend = (t1 + 30)*Fs;
end


figure
offset = 200;
for i = 1:n
    plot(t(tstart:tend)/60, X_EEG(i,tstart:tend) - offset*(i-1))
    hold on
end
xlabel('t [min]')
ylim([-n*offset, offset])
set(gca, 'yTick', (-offset*(n-1):offset:0)) % Nota: gca : Current Axes
labels = {'F3';'F4';'T7';'C3';'Cz';'C4';'T8';'PO7';'PO3';'PO4';'PO8';'O1';'O2'};
set(gca, 'yTickLabel', flipud(labels))
title('EEG')

figure
plot(t(tstart:tend)/60,X_ECG(tstart:tend))
title('ECG')
xlabel('t [min]')

%% Visualize power spectral density:
window = 5*Fs;
nFFT = 2*window; % per avere risoluzione di 0.1 Hz

[PSD_EEG,F] = pwelch(X_EEG',window,[],nFFT,Fs);

[PSD_ECG,F] = pwelch(X_ECG,window,[],nFFT,Fs);

Fstart = 0.5*10+1; % [Hz]
Fend = 40*10+1;    % [Hz]
figure
for i = 1:n
    subplot(3,5,i)
    plot(F(Fstart:Fend),PSD_EEG(Fstart:Fend,i))
    title(['PSD ', labels(i)])
end

figure
plot(F(Fstart:Fend),PSD_ECG(Fstart:Fend))
title('PSD of ECG')

%%
% Using EEGLAB it's possible to estimate the demixing matrix W and export it to a .txt file.
switch patient 
    case 1
    load W_patient1.txt
    W = W_patient1;
    case 2
    load W_patient2.txt
    W = W_patient2;
end

% Using the demixing matrix it's possible to reconstruct the 13 indipendent components:
S = W*X_EEG; % Estimated sources
A = inv(W);

figure
for i = 1:n
    plot(t(tstart:tend)/60,S(i,tstart:tend) - offset*(i-1))
    hold on
end
xlabel('t [min]')
title(['Sources (IC) of the patient ',num2str(patient)])

ylim([-n*offset, offset])
set(gca, 'yTick', (-offset*(n-1):offset:0)) % Nota: gca : Curren Axes
labels = {'IC1','IC2','IC3','IC4','IC5','IC6','IC7','IC8','IC9','IC10','IC11','IC12','IC13'};
set(gca, 'yTickLabel', fliplr(labels))

%% Power spectral density of the indipendent components:
[PSD_S,F] = pwelch(S',window,[],nFFT,Fs);

figure
for i = 1:n
    subplot(3,5,i)
    plot(F(Fstart:Fend),PSD_S(Fstart:Fend,i))
    title(['PSD IC',num2str(i)])
end
% Note that the spectrum relative to blinking has only low-frequency components.

%% Reconstruct EEG signals without artefacts:
Snew = S;

switch patient 
    case 1
    Snew([1,2,5],:) = 0;
    case 2
    Snew([1,2,3],:) = 0;
end

Xnew = A*Snew;

figure
for i = 1:n
    plot(t(tstart:tend)/60, Xnew(i,tstart:tend) - offset*(i-1))
    hold on
end
xlabel('t [min]')
ylim([-n*offset, offset])
set(gca, 'yTick', (-offset*(n-1):offset:0)) % Nota: gca : Current Axes
labels = {'F3';'F4';'T7';'C3';'Cz';'C4';'T8';'PO7';'PO3';'PO4';'PO8';'O1';'O2'};
set(gca, 'yTickLabel', flipud(labels))
title(['EEG denoised of the patient ',num2str(patient)])

figure
for i = 1:n
    subplot(3,5,i)
    plot(F(Fstart:Fend),PSD_EEG(Fstart:Fend,i))
    title(['PSD ', labels(i)])
end
% Note that alpha rhythm (around 10 Hz) appears mainly in posterior channels, parietal (P) and occipital (O).

%% Using EEG denoised (without artefacts), estimate separately the PSD of the initial relax phase (R1), of the task (T) and the final relax phase (R2):

% First relax phase (first 5 minutes)
start_r1 = 1;
end_r1 = 5*60*Fs;
[PSD_r1,F] = pwelch(Xnew(:,start_r1:end_r1)',window,[],nFFT,Fs);

% Intermediate task phase (from 5 to 10 minutes)
start_T = end_r1 +1;
end_T = 10*60*Fs;
[PSD_T,F] = pwelch(Xnew(:,start_T:end_T)',window,[],nFFT,Fs);

% Final relax phase (last 5 minutes)
start_r2 = end_T +1;
end_r2 = m;
[PSD_r2,F] = pwelch(Xnew(:,start_r2:end_r2)',window,[],nFFT,Fs);

figure
subplot(131)
plot(F(Fstart:Fend),PSD_r1(Fstart:Fend,:))
title('PSD phase R1')
legend(labels)
subplot(132)
plot(F(Fstart:Fend),PSD_T(Fstart:Fend,:))
title('PSD phase T')
legend(labels)
subplot(133)
plot(F(Fstart:Fend),PSD_r2(Fstart:Fend,:))
title('PSD phase R2')
legend(labels)
% Note: comparing the three plots, it's possible to see that during the task 
% there's a reduction of the alpha rhythm, because the patient is concentrated.

%% For each phase, compute the average PSD in frontal channels (F3, F4), in temporal/central channels (T7, C3, Cz, C4, T8), in parietal/occipital channels (PO7, PO3, PO4, PO8, O1, O2).

ind_F = [1,2]; % F3, F4
ind_CT = [3:7]; % T7, C3, Cz, C4, T8
ind_PO = [8:13]; % PO7, PO3, PO4, PO8, O1, O2

PSD_r1_F = mean(PSD_r1(:,ind_F),2);
PSD_r1_CT = mean(PSD_r1(:,ind_CT),2);
PSD_r1_PO = mean(PSD_r1(:,ind_PO),2);

PSD_T_F = mean(PSD_T(:,ind_F),2);
PSD_T_CT = mean(PSD_T(:,ind_CT),2);
PSD_T_PO = mean(PSD_T(:,ind_PO),2);

PSD_r2_F = mean(PSD_r2(:,ind_F),2);
PSD_r2_CT = mean(PSD_r2(:,ind_CT),2);
PSD_r2_PO = mean(PSD_r2(:,ind_PO),2);

figure
subplot(311)
plot(F(Fstart:Fend),PSD_r1_F(Fstart:Fend))
hold on
plot(F(Fstart:Fend),PSD_T_F(Fstart:Fend))
plot(F(Fstart:Fend),PSD_r2_F(Fstart:Fend))
plot(8*ones(1,20),linspace(0,max(PSD_r1_F(Fstart:Fend)),20),'g--')
plot(14*ones(1,20),linspace(0,max(PSD_r1_F(Fstart:Fend)),20),'g--')
title('PSD frontal (F)')
legend('Relax 1','Task','Relax 2','Alpha band')
xlabel('f [Hz]')

subplot(312)
plot(F(Fstart:Fend),PSD_r1_CT(Fstart:Fend))
hold on
plot(F(Fstart:Fend),PSD_T_CT(Fstart:Fend))
plot(F(Fstart:Fend),PSD_r2_CT(Fstart:Fend))
plot(8*ones(1,20),linspace(0,max(PSD_r1_CT(Fstart:Fend)),20),'g--')
plot(14*ones(1,20),linspace(0,max(PSD_r1_CT(Fstart:Fend)),20),'g--')
title('PSD cantral/temporal (CT)')
legend('Relax 1','Task','Relax 2','Alpha band')
xlabel('f [Hz]')

subplot(313)
plot(F(Fstart:Fend),PSD_r1_PO(Fstart:Fend))
hold on
plot(F(Fstart:Fend),PSD_T_PO(Fstart:Fend))
plot(F(Fstart:Fend),PSD_r2_PO(Fstart:Fend))
plot(8*ones(1,20),linspace(0,max(PSD_r2_PO(Fstart:Fend)),20),'g--')
plot(14*ones(1,20),linspace(0,max(PSD_r2_PO(Fstart:Fend)),20),'g--')
title('PSD posterior/occipital (PO)')
legend('Relax 1','Task','Relax 2','Alpha band')
xlabel('f [Hz]')

%% Compute the PSD in the alpha band (8-14 Hz) for each region (frontal, temporal/central, parietal/occipital) in each of the three phases (R1, T, R2).

imin = find(F==8);
imax = find(F==14);

ALPHA_r1_F = trapz(F(imin:imax),PSD_r1_F(imin:imax));
ALPHA_T_F = trapz(F(imin:imax),PSD_T_F(imin:imax));
ALPHA_r2_F = trapz(F(imin:imax),PSD_r2_F(imin:imax));

ALPHA_r1_CT = trapz(F(imin:imax),PSD_r1_CT(imin:imax));
ALPHA_T_CT = trapz(F(imin:imax),PSD_T_CT(imin:imax));
ALPHA_r2_CT = trapz(F(imin:imax),PSD_r2_CT(imin:imax));

ALPHA_r1_PO = trapz(F(imin:imax),PSD_r1_PO(imin:imax));
ALPHA_T_PO = trapz(F(imin:imax),PSD_T_PO(imin:imax));
ALPHA_r2_PO = trapz(F(imin:imax),PSD_r2_PO(imin:imax));

figure
plot([1,2,3],[ALPHA_r1_F,ALPHA_T_F,ALPHA_r2_F],'-og')
hold on
plot([1,2,3],[ALPHA_r1_CT,ALPHA_T_CT,ALPHA_r2_CT],'-ob')
plot([1,2,3],[ALPHA_r1_PO,ALPHA_T_PO,ALPHA_r2_PO],'-or')

legend('F','CT','PO')
xlim([0,4])
set(gca, 'xTick', (1:3)) % Note: gca : Current Axes
phase = {'R1','T','R2'};
set(gca, 'xTickLabel', phase)
title('Power in alpha band')
% Note: during the task the energy of the alpha wave is very low, 
% while during the phases R1 and R2 is higher. 
% Also note that there is more energy in the signals from the rear electrodes than in the front ones.

%% To emphasize the independent component that contributes to the EEG spectrum in alpha band, it recalculates the spectrum of the denoised EEG signals (without artifacts) by removing one or two components at a time

switch patient 
    case 1
        Snew(4,:) = 0; % Remove IC4 - for the first patient, I noticed that IC4 gives an important contribute to alpha rhythm.
    case 2
        Snew(7,:) = 0; % Remove IC7 - for the second patient, IC7 is important for alpha rhythm in central/temporal zone.
end


figure
for i = 1:n
    plot(t(tstart:tend)/60,Snew(i,tstart:tend) - offset*(i-1))
    hold on
end
xlabel('t [min]')
title(['Sorgenti (IC) of the patient ',num2str(patient)])

ylim([-n*offset, offset])
set(gca, 'yTick', (-offset*(n-1):offset:0))
labels = {'IC1','IC2','IC3','IC4','IC5','IC6','IC7','IC8','IC9','IC10','IC11','IC12','IC13'};
set(gca, 'yTickLabel', fliplr(labels))

Xnew = A*Snew;

[PSD_r1,F] = pwelch(Xnew(:,start_r1:end_r1)',window,[],nFFT,Fs);
[PSD_T,F] = pwelch(Xnew(:,start_T:end_T)',window,[],nFFT,Fs);
[PSD_r2,F] = pwelch(Xnew(:,start_r2:end_r2)',window,[],nFFT,Fs);

figure
subplot(131)
plot(F(Fstart:Fend),PSD_r1(Fstart:Fend,:))
title('PSD phase R1')
legend(labels)
subplot(132)
plot(F(Fstart:Fend),PSD_T(Fstart:Fend,:))
title('PSD phase T')
legend(labels)
subplot(133)
plot(F(Fstart:Fend),PSD_r2(Fstart:Fend,:))
title('PSD phase R2')
legend(labels)

PSD_r1_F = mean(PSD_r1(:,ind_F),2);
PSD_r1_CT = mean(PSD_r1(:,ind_CT),2);
PSD_r1_PO = mean(PSD_r1(:,ind_PO),2);

PSD_T_F = mean(PSD_T(:,ind_F),2);
PSD_T_CT = mean(PSD_T(:,ind_CT),2);
PSD_T_PO = mean(PSD_T(:,ind_PO),2);

PSD_r2_F = mean(PSD_r2(:,ind_F),2);
PSD_r2_CT = mean(PSD_r2(:,ind_CT),2);
PSD_r2_PO = mean(PSD_r2(:,ind_PO),2);

figure
subplot(311)
plot(F(Fstart:Fend),PSD_r1_F(Fstart:Fend))
hold on
plot(F(Fstart:Fend),PSD_T_F(Fstart:Fend))
plot(F(Fstart:Fend),PSD_r2_F(Fstart:Fend))
plot(8*ones(1,20),linspace(0,max(PSD_r1_F(Fstart:Fend)),20),'g--')
plot(14*ones(1,20),linspace(0,max(PSD_r1_F(Fstart:Fend)),20),'g--')
title('PSD frontal (F)')
legend('Relax 1','Task','Relax 2','Alpha band')
xlabel('f [Hz]')

subplot(312)
plot(F(Fstart:Fend),PSD_r1_CT(Fstart:Fend))
hold on
plot(F(Fstart:Fend),PSD_T_CT(Fstart:Fend))
plot(F(Fstart:Fend),PSD_r2_CT(Fstart:Fend))
plot(8*ones(1,20),linspace(0,max(PSD_r1_CT(Fstart:Fend)),20),'g--')
plot(14*ones(1,20),linspace(0,max(PSD_r1_CT(Fstart:Fend)),20),'g--')
title('PSD cantral/temporal (CT)')
legend('Relax 1','Task','Relax 2','Alpha band')
xlabel('f [Hz]')

subplot(313)
plot(F(Fstart:Fend),PSD_r1_PO(Fstart:Fend))
hold on
plot(F(Fstart:Fend),PSD_T_PO(Fstart:Fend))
plot(F(Fstart:Fend),PSD_r2_PO(Fstart:Fend))
plot(8*ones(1,20),linspace(0,max(PSD_r2_PO(Fstart:Fend)),20),'g--')
plot(14*ones(1,20),linspace(0,max(PSD_r2_PO(Fstart:Fend)),20),'g--')
title('PSD posterior/occipital (PO)')
legend('Relax 1','Task','Relax 2','Alpha band')
xlabel('f [Hz]')
% note that by removing an independent component (IC4 for the first patient,
% or IC7 for the second patient), it significantly decrease the power in 
% the alpha band, which means that that component gives an important contribution to the alpha wave.

ALFA_r1_F = trapz(F(imin:imax),PSD_r1_F(imin:imax));
ALPHA_T_F = trapz(F(imin:imax),PSD_T_F(imin:imax));
ALPHA_r2_F = trapz(F(imin:imax),PSD_r2_F(imin:imax));

ALPHA_r1_CT = trapz(F(imin:imax),PSD_r1_CT(imin:imax));
ALPHA_T_CT = trapz(F(imin:imax),PSD_T_CT(imin:imax));
ALPHA_r2_CT = trapz(F(imin:imax),PSD_r2_CT(imin:imax));

ALPHA_r1_PO = trapz(F(imin:imax),PSD_r1_PO(imin:imax));
ALPHA_T_PO = trapz(F(imin:imax),PSD_T_PO(imin:imax));
ALPHA_r2_PO = trapz(F(imin:imax),PSD_r2_PO(imin:imax));

figure
plot([1,2,3],[ALPHA_r1_F,ALPHA_T_F,ALPHA_r2_F],'-og')
hold on
plot([1,2,3],[ALPHA_r1_CT,ALPHA_T_CT,ALPHA_r2_CT],'-ob')
plot([1,2,3],[ALPHA_r1_PO,ALPHA_T_PO,ALPHA_r2_PO],'-or')

legend('F','CT','PO')
xlim([0,4])
set(gca, 'xTick', (1:3)) % Nota: gca : Curren Axes
phase = {'R1','T','R2'};
set(gca, 'xTickLabel', phase)
title('Power in alpha bandwidth')