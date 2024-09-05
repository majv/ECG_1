% Equipo: Los Laplacianos 
clc; close all; clear;

% Obtener los datos del archivo de teams
Datos_ECG = load('Prac3_ECG_PARTE_1.mat');

% Datos para las graficas
Fs_Data = Datos_ECG.Fs;
No_Muestras = Fs_Data * 50;
Tiempo = linspace(0, 50, No_Muestras);

% Para el primer sujeto, solo para los primeros 50 segundos
ECG_REST = Datos_ECG.ECG_REST(1:No_Muestras);
ECG_Mat = Datos_ECG.ECG_MAT(1:No_Muestras);

% Para el segundo sujeto, solo para los primeros 50 segundos
ECG_REST2 = Datos_ECG.ECG_REST_2(1:No_Muestras);
ECG_Mat2 = Datos_ECG.ECG_MAT_2(1:No_Muestras);

%% Parte 2.- Graficar los resultados
figure;
% Px 1
subplot(2,2,1);
    plot(Tiempo, ECG_REST, 'm');
    grid on;
    title('ECG Rest vs Tiempo Px_1');
    ylabel('Voltaje [mV]');
    xlabel('Tiempo [s]');
    xlim([0 5]); % Limitar el eje x a 5 segundos

subplot(2,2,2);
    plot(Tiempo, ECG_Mat, 'm');
    grid on;
    title('ECG Mat vs Tiempo Px_1');
    ylabel('Voltaje [mV]');
    xlabel('Tiempo [s]');
    xlim([0 5]);

% Px 2
subplot(2,2,3);
    plot(Tiempo, ECG_REST2, 'b');
    grid on;
    title('ECG Rest vs Tiempo Px_2');
    ylabel('Voltaje [mV]');
    xlabel('Tiempo [s]');
    xlim([0 5]); 

subplot(2,2,4);
    plot(Tiempo, ECG_Mat2, 'b');
    grid on;
    title('ECG Mat vs Tiempo Px_2');
    ylabel('Voltaje [mV]');
    xlabel('Tiempo [s]');
    xlim([0 5]); 

%% Parte 3.- Identificar picos
% Participante 1 (Px1 con ECG_REST)
[picos1_1, ubicaciones1_1] = findpeaks(ECG_REST, 'MinPeakHeight', 0.35); 
tiempo_picos1_1 = Tiempo(ubicaciones1_1);

% Participante 1 (Px1 con ECG_Mat)
[picos1_2, ubicaciones1_2] = findpeaks(ECG_Mat, 'MinPeakHeight', 0.35); 
tiempo_picos1_2 = Tiempo(ubicaciones1_2);

% Participante 2 (Px2 con ECG_REST)
[picos2_1, ubicaciones2_1] = findpeaks(ECG_REST2, 'MinPeakHeight', 0.35); 
tiempo_picos2_1 = Tiempo(ubicaciones2_1);

% Participante 2 (Px2 con ECG_Mat)
[picos2_2, ubicaciones2_2] = findpeaks(ECG_Mat2, 'MinPeakHeight', 0.35); 
tiempo_picos2_2 = Tiempo(ubicaciones2_2);

%% Parte 4.- Grafica en una sola figura   
figure();
subplot(2,2,1);
    plot(Tiempo, ECG_REST,'b');
    hold on;
    plot(tiempo_picos1_1, picos1_1, 'r*');  % Marca los picos con un *
    grid on;
    title('Participante 1: ECG Rest vs Tiempo con Picos');
    xlabel('Tiempo [s]');
    ylabel('Voltaje [mV]');
    legend('ECG REST', 'Picos');
    xlim([0 5]); 

subplot(2,2,2);
    plot(Tiempo, ECG_REST2,'b');
    hold on;
    plot(tiempo_picos2_1, picos2_1, 'r*'); 
    grid on;
    title('Participante 2: ECG Rest vs Tiempo con Picos');
    xlabel('Tiempo [s]');
    ylabel('Voltaje [mV]');
    legend('ECG REST', 'Picos');
    xlim([0 5]); 

subplot(2,2,3);
    plot(Tiempo, ECG_Mat,'b');
    hold on;
    plot(tiempo_picos1_2, picos1_2, 'r*'); 
    grid on;
    title('Participante 1: ECG Mat vs Tiempo con Picos');
    xlabel('Tiempo [s]');
    ylabel('Voltaje [mV]');
    legend('ECG Mat', 'Picos');
    xlim([0 5]); 

subplot(2,2,4);
    plot(Tiempo, ECG_Mat2,'b');
    hold on;
    plot(tiempo_picos2_2, picos2_2, 'r*');  
    grid on;
    title('Participante 2: ECG Mat vs Tiempo con Picos');
    xlabel('Tiempo [s]');
    ylabel('Voltaje [mV]');
    legend('ECG Mat', 'Picos');
    xlim([0 5]); 
    %% Parte 5 Tacograma
    
% Participante 1 (Px1 con ECG_REST)
RR_intervals1_1 = diff(tiempo_picos1_1);
tacograma1_1 = RR_intervals1_1 * 1000; 

% Participante 1 (Px1 con ECG_Mat)
RR_intervals1_2 = diff(tiempo_picos1_2);
tacograma1_2 = RR_intervals1_2 * 1000; 

% Participante 2 (Px2 con ECG_REST)
RR_intervals2_1 = diff(tiempo_picos2_1);
tacograma2_1 = RR_intervals2_1 * 1000; 

% Participante 2 (Px2 con ECG_Mat)
RR_intervals2_2 = diff(tiempo_picos2_2);
tacograma2_2 = RR_intervals2_2 * 1000; 

%% Parte 5.1- Graficar Taco en una sola figura
figure('Name', 'Tacogramas', 'NumberTitle', 'off');

subplot(2,2,1);
plot(tiempo_picos1_1(1:end-1), tacograma1_1, 'bx-');
title('Tacograma Participante 1: ECG REST');
xlabel('Tiempo [s]');
ylabel('Intervalo RR [ms]');
grid on;

subplot(2,2,2);
plot(tiempo_picos1_2(1:end-1), tacograma1_2, 'rx-');
title('Tacograma Participante 1: ECG MAT');
xlabel('Tiempo [s]');
ylabel('Intervalo RR [ms]');
grid on;

subplot(2,2,3);
plot(tiempo_picos2_1(1:end-1), tacograma2_1, 'gx-');
title('Tacograma Participante 2: ECG REST');
xlabel('Tiempo [s]');
ylabel('Intervalo RR [ms]');
grid on;

subplot(2,2,4);
plot(tiempo_picos2_2(1:end-1), tacograma2_2, 'mx-');
title('Tacograma Participante 2: ECG MAT');
xlabel('Tiempo [s]');
ylabel('Intervalo RR [ms]');
grid on;

sgtitle('Tacogramas de los Participantes en Reposo y Matemáticas');
%% Punto 6 Calculo de parametros

% Participante 1 - ECG REST
RR_promedio1_1 = mean(tacograma1_1);
SDNN1_1 = std(tacograma1_1);
diff_RR1_1 = diff(tacograma1_1);
PNN50_1_1 = sum(abs(diff_RR1_1) > 50) / length(diff_RR1_1) * 100;
RMSSD1_1 = sqrt(mean(diff_RR1_1.^2));
FC_promedio1_1 = 60 / (RR_promedio1_1 / 1000);

% Participante 1 - ECG MAT
RR_promedio1_2 = mean(tacograma1_2);
SDNN1_2 = std(tacograma1_2);
diff_RR1_2 = diff(tacograma1_2);
PNN50_1_2 = sum(abs(diff_RR1_2) > 50) / length(diff_RR1_2) * 100;
RMSSD1_2 = sqrt(mean(diff_RR1_2.^2));
FC_promedio1_2 = 60 / (RR_promedio1_2 / 1000);

% Participante 2 - ECG REST
RR_promedio2_1 = mean(tacograma2_1);
SDNN2_1 = std(tacograma2_1);
diff_RR2_1 = diff(tacograma2_1);
PNN50_2_1 = sum(abs(diff_RR2_1) > 50) / length(diff_RR2_1) * 100;
RMSSD2_1 = sqrt(mean(diff_RR2_1.^2));
FC_promedio2_1 = 60 / (RR_promedio2_1 / 1000);

% Participante 2 - ECG MAT
RR_promedio2_2 = mean(tacograma2_2);
SDNN2_2 = std(tacograma2_2);
diff_RR2_2 = diff(tacograma2_2);
PNN50_2_2 = sum(abs(diff_RR2_2) > 50) / length(diff_RR2_2) * 100;
RMSSD2_2 = sqrt(mean(diff_RR2_2.^2));
FC_promedio2_2 = 60 / (RR_promedio2_2 / 1000);

%% Parte 6.1 Mostramos la tabla
fprintf('Participante 1:\n');
fprintf('| Parámetro | ECG REST | ECG MAT |\n');
fprintf('|-----------|----------|--------|\n');
fprintf('| RR promedio (ms) | %.2f | %.2f |\n', RR_promedio1_1, RR_promedio1_2);
fprintf('| SDNN (ms) | %.2f | %.2f |\n', SDNN1_1, SDNN1_2);
fprintf('| PNN50 (%%) | %.2f | %.2f |\n', PNN50_1_1, PNN50_1_2);
fprintf('| RMSSD (ms) | %.2f | %.2f |\n', RMSSD1_1, RMSSD1_2);
fprintf('| FC promedio (bpm) | %.2f | %.2f |\n\n', FC_promedio1_1, FC_promedio1_2);

fprintf('Participante 2:\n');
fprintf('| Parámetro | ECG REST | ECG MAT |\n');
fprintf('|-----------|----------|--------|\n');
fprintf('| RR promedio (ms) | %.2f | %.2f |\n', RR_promedio2_1, RR_promedio2_2);
fprintf('| SDNN (ms) | %.2f | %.2f |\n', SDNN2_1, SDNN2_2);
fprintf('| PNN50 (%%) | %.2f | %.2f |\n', PNN50_2_1, PNN50_2_2);
fprintf('| RMSSD (ms) | %.2f | %.2f |\n', RMSSD2_1, RMSSD2_2);
fprintf('| FC promedio (bpm) | %.2f | %.2f |\n', FC_promedio2_1, FC_promedio2_2);


%% Punto 7 Interpolación
% Frecuencia nueva
Fs_new = 4; 

% Vector de tiempo 
t_new = 0:1/Fs_new:50-1/Fs_new;

% Interpolar las señales 
ECG_REST_interp = interp1(Tiempo, ECG_REST, t_new);
ECG_Mat_interp = interp1(Tiempo, ECG_Mat, t_new);

% Calcular el espectro 
[Pxx_rest, F_rest] = pwelch(ECG_REST_interp, [], [], [], Fs_new);
[Pxx_mat, F_mat] = pwelch(ECG_Mat_interp, [], [], [], Fs_new);

% Definir LF y HF
LF_band = [0.04, 0.15];  
HF_band = [0.15, 0.4];   

LF_idx_rest = find(F_rest >= 0.04 & F_rest <= 0.15);  
HF_idx_rest = find(F_rest > 0.15 & F_rest <= 0.4);  

%Graficar
figure;

% ECG REST
subplot(1, 2, 1);
plot(F_rest, Pxx_rest, 'b');
hold on;

fill([F_rest(LF_idx_rest); flipud(F_rest(LF_idx_rest))], ...
     [Pxx_rest(LF_idx_rest); zeros(size(LF_idx_rest))], 'm', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

fill([F_rest(HF_idx_rest); flipud(F_rest(HF_idx_rest))], ...
     [Pxx_rest(HF_idx_rest); zeros(size(HF_idx_rest))], 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

title('Espectro de Potencia Participante 1: ECG REST');
xlabel('Frecuencia [Hz]');
ylabel('Potencia [mV²]');
xline(0.04, 'm--', 'LF');
xline(0.4, 'c--', 'HF');
grid on;

% ECG MAT
subplot(1, 2, 2);
plot(F_mat, Pxx_mat, 'b');
hold on;

LF_idx_mat = find(F_mat >= 0.04 & F_mat <= 0.15);
HF_idx_mat = find(F_mat > 0.15 & F_mat <= 0.4);


fill([F_mat(LF_idx_mat); flipud(F_mat(LF_idx_mat))], ...
     [Pxx_mat(LF_idx_mat); zeros(size(LF_idx_mat))], 'm', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

fill([F_mat(HF_idx_mat); flipud(F_mat(HF_idx_mat))], ...
     [Pxx_mat(HF_idx_mat); zeros(size(HF_idx_mat))], 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

title('Espectro de Potencia Participante 1: ECG MAT');
xlabel('Frecuencia [Hz]');
ylabel('Potencia [mV²]');
xline(0.04, 'm--', 'LF');
xline(0.4, 'c--', 'HF');
grid on;

sgtitle('Espectro de Potencia para ECG REST y ECG MAT Participante 1');

%% Punto 8 Ruido

% Generar ruido blanco
ruido_blanco = 0.1 * randn(size(ECG_REST));

% Llamar a la función ecgnoise para generar la señal con ruido de 60 Hz
[t, ECG_ruido_60hz] = ecgnoise(); 

% Asegurar que el tamaño de las señales coincida
minLength = min(length(ECG_REST), length(ECG_ruido_60hz));
ECG_REST = ECG_REST(1:minLength);
ECG_ruido_60hz = ECG_ruido_60hz(1:minLength);
Tiempo = Tiempo(1:minLength);

% Señal contaminada con ruido blanco
ECG_ruido_blanco = ECG_REST + ruido_blanco;

% Calcular espectros de potencia
nfft = length(ECG_REST);
f = Fs_Data * (0:(nfft/2)) / nfft;

% Espectro de potencia para señal con ruido blanco
Y_blanco = fft(ECG_ruido_blanco, nfft);
P2_blanco = abs(Y_blanco / nfft);
P1_blanco = P2_blanco(1:nfft/2+1);
P1_blanco(2:end-1) = 2 * P1_blanco(2:end-1);

% Espectro de potencia para señal con ruido de 60 Hz
Y_60hz = fft(ECG_ruido_60hz, nfft);
P2_60hz = abs(Y_60hz / nfft);
P1_60hz = P2_60hz(1:nfft/2+1);
P1_60hz(2:end-1) = 2 * P1_60hz(2:end-1);

% Graficar en una sola figura
figure('Name', 'ECG con Ruido', 'NumberTitle', 'off');

% Lado izquierdo superior: ECG con ruido de 60 Hz
subplot(2,2,1);
plot(Tiempo, ECG_ruido_60hz);
title('ECG con ruido de 60 Hz');
xlabel('Tiempo [s]');
ylabel('Voltaje [mV]');
xlim([0 5]);
grid on;

% Lado izquierdo inferior: Espectro de potencia de ECG con ruido de 60 Hz
subplot(2,2,3);
plot(f, 10*log10(P1_60hz));
title('Espectro de Potencia: ECG con ruido de 60 Hz');
xlabel('Frecuencia [Hz]');
ylabel('Densidad de Potencia [dB]');
grid on;
xlim([0 100]); % Limitamos a 100 Hz para mejor visualización

% Lado derecho superior: ECG con ruido blanco
subplot(2,2,2);
plot(Tiempo, ECG_ruido_blanco);
title('ECG con ruido blanco');
xlabel('Tiempo [s]');
ylabel('Voltaje [mV]');
xlim([0 5]);
grid on;

% Lado derecho inferior: Espectro de potencia de ECG con ruido blanco
subplot(2,2,4);
plot(f, 10*log10(P1_blanco));
title('Espectro de Potencia: ECG con ruido blanco');
xlabel('Frecuencia [Hz]');
ylabel('Densidad de Potencia [dB]');
grid on;
xlim([0 100]); % Limitamos a 100 Hz para mejor visualización

% Ajustar el diseño de la figura
sgtitle('ECG con Ruido Participante 1: ECG REST');
