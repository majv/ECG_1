clc; close all; clear;

% Cargar datos
load('medio.mat');  
b_medio = data(:, 1);  
t_medio = data(:, 2);

load('tres.mat');  
b_tres = data(:, 1);  
t_tres = data(:, 2);

fs = 6250;  
tiempo_medio = (0:length(b_medio)-1) / fs;
tiempo_tres = (0:length(b_tres)-1) / fs;

segmentos_medio = [3600 12100; 17000 25500; 28500 37000; 40500 49000; 52200 60700];
segmentos_tres = [3000 18000; 21700 36700; 38400 53400; 56800 71800; 77500 92500];

% Graficar EMG
figure;
subplot(2,1,1); 
plot(tiempo_medio, b_medio);
title('EMG del Bíceps con Medio Kilo de Carga');  
xlabel('Tiempo [s]');  
ylabel('Amplitud [V]'); 

subplot(2,1,2);  
plot(tiempo_medio, t_medio);
title('EMG del Tríceps con Medio Kilo de Carga'); 
xlabel('Tiempo [s]');  
ylabel('Amplitud [V]');  

figure;
subplot(2,1,1); 
plot(tiempo_tres, b_tres);
title('EMG del Bíceps con Tres Kilos de Carga');  
xlabel('Tiempo [s]');  
ylabel('Amplitud [V]'); 

subplot(2,1,2);  
plot(tiempo_tres, t_tres);
title('EMG del Tríceps con Tres Kilos de Carga'); 
xlabel('Tiempo [s]');  
ylabel('Amplitud [V]');  

figure;
for i = 1:5
    subplot(2, 5, i);  
    plot(tiempo_medio(segmentos_medio(i,1):segmentos_medio(i,2)), ...
         b_medio(segmentos_medio(i,1):segmentos_medio(i,2)));
    title(['Segmento ' num2str(i) ' Bíceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    
    subplot(2, 5, i+5);  
    plot(tiempo_medio(segmentos_medio(i,1):segmentos_medio(i,2)), ...
         t_medio(segmentos_medio(i,1):segmentos_medio(i,2)));
    title(['Segmento ' num2str(i) ' Tríceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
end

figure;
for i = 1:5
    subplot(2, 5, i);  
    plot(tiempo_tres(segmentos_tres(i,1):segmentos_tres(i,2)), ...
         b_tres(segmentos_tres(i,1):segmentos_tres(i,2)));
    title(['Segmento ' num2str(i) ' Bíceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    
    subplot(2, 5, i+5);  
    plot(tiempo_tres(segmentos_tres(i,1):segmentos_tres(i,2)), ...
         t_tres(segmentos_tres(i,1):segmentos_tres(i,2)));
    title(['Segmento ' num2str(i) ' Tríceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
end

% Función para calcular RMS
calcRMS = @(x) sqrt(mean(x.^2));

% Inicializar matrices de resultados
rms_medio_biceps = zeros(1, 5);
rms_medio_triceps = zeros(1, 5);
rms_tres_biceps = zeros(1, 5);
rms_tres_triceps = zeros(1, 5);

area_medio_biceps = zeros(1, 5);
area_medio_triceps = zeros(1, 5);
area_tres_biceps = zeros(1, 5);
area_tres_triceps = zeros(1, 5);

cruces_medio_biceps = zeros(1, 5);
cruces_medio_triceps = zeros(1, 5);
cruces_tres_biceps = zeros(1, 5);
cruces_tres_triceps = zeros(1, 5);

% Calcular métricas para medio kilo
for i = 1:5
    segmento_biceps_medio = b_medio(segmentos_medio(i,1):segmentos_medio(i,2));
    rms_medio_biceps(i) = calcRMS(segmento_biceps_medio);
    
    segmento_triceps_medio = t_medio(segmentos_medio(i,1):segmentos_medio(i,2));
    rms_medio_triceps(i) = calcRMS(segmento_triceps_medio);
    
    env_biceps_medio = envelope(abs(segmento_biceps_medio), 100, 'peak');
    area_medio_biceps(i) = trapz(tiempo_medio(segmentos_medio(i,1):segmentos_medio(i,2)), env_biceps_medio);
    
    env_triceps_medio = envelope(abs(segmento_triceps_medio), 100, 'peak');
    area_medio_triceps(i) = trapz(tiempo_medio(segmentos_medio(i,1):segmentos_medio(i,2)), env_triceps_medio);
    
    cruces_medio_biceps(i) = sum(diff(segmento_biceps_medio > 0) ~= 0);
    cruces_medio_triceps(i) = sum(diff(segmento_triceps_medio > 0) ~= 0);
end

% Calcular métricas para tres kilos
for i = 1:5
    segmento_biceps_tres = b_tres(segmentos_tres(i,1):segmentos_tres(i,2));
    rms_tres_biceps(i) = calcRMS(segmento_biceps_tres);
    
    segmento_triceps_tres = t_tres(segmentos_tres(i,1):segmentos_tres(i,2));
    rms_tres_triceps(i) = calcRMS(segmento_triceps_tres);
    
    env_biceps_tres = envelope(abs(segmento_biceps_tres), 100, 'peak');
    area_tres_biceps(i) = trapz(tiempo_tres(segmentos_tres(i,1):segmentos_tres(i,2)), env_biceps_tres);
    
    env_triceps_tres = envelope(abs(segmento_triceps_tres), 100, 'peak');
    area_tres_triceps(i) = trapz(tiempo_tres(segmentos_tres(i,1):segmentos_tres(i,2)), env_triceps_tres);
    
    cruces_tres_biceps(i) = sum(diff(segmento_biceps_tres > 0) ~= 0);
    cruces_tres_triceps(i) = sum(diff(segmento_triceps_tres > 0) ~= 0);
end

% Mostrar resultados en la consola
disp('Resultados para medio kilo:');
disp('RMS del Bíceps:');
disp(rms_medio_biceps);
disp('RMS del Tríceps:');
disp(rms_medio_triceps);
disp('Área bajo la curva del Bíceps:');
disp(area_medio_biceps);
disp('Área bajo la curva del Tríceps:');
disp(area_medio_triceps);
disp('Número de Cruces por Cero del Bíceps:');
disp(cruces_medio_biceps);
disp('Número de Cruces por Cero del Tríceps:');
disp(cruces_medio_triceps);

disp('Resultados para tres kilos:');
disp('RMS del Bíceps:');
disp(rms_tres_biceps);
disp('RMS del Tríceps:');
disp(rms_tres_triceps);
disp('Área bajo la curva del Bíceps:');
disp(area_tres_biceps);
disp('Área bajo la curva del Tríceps:');
disp(area_tres_triceps);
disp('Número de Cruces por Cero del Bíceps:');
disp(cruces_tres_biceps);
disp('Número de Cruces por Cero del Tríceps:');
disp(cruces_tres_triceps);

% Calculamos el valor promedio y la mediana para cada medida
promedio_rms_biceps_medio = mean(rms_medio_biceps);
mediana_rms_biceps_medio = median(rms_medio_biceps);

promedio_rms_biceps_tres = mean(rms_tres_biceps);
mediana_rms_biceps_tres = median(rms_tres_biceps);

promedio_rms_triceps_medio = mean(rms_medio_triceps);
mediana_rms_triceps_medio = median(rms_medio_triceps);

promedio_rms_triceps_tres = mean(rms_tres_triceps);
mediana_rms_triceps_tres = median(rms_tres_triceps);

promedio_auc_biceps_medio = mean(area_medio_biceps);
mediana_auc_biceps_medio = median(area_medio_biceps);

promedio_auc_biceps_tres = mean(area_tres_biceps);
mediana_auc_biceps_tres = median(area_tres_biceps);

promedio_auc_triceps_medio = mean(area_medio_triceps);
mediana_auc_triceps_medio = median(area_medio_triceps);

promedio_auc_triceps_tres = mean(area_tres_triceps);
mediana_auc_triceps_tres = median(area_tres_triceps);

promedio_cruces_biceps_medio = mean(cruces_medio_biceps);
mediana_cruces_biceps_medio = median(cruces_medio_biceps);

promedio_cruces_biceps_tres = mean(cruces_tres_biceps);
mediana_cruces_biceps_tres = median(cruces_tres_biceps);

promedio_cruces_triceps_medio = mean(cruces_medio_triceps);
mediana_cruces_triceps_medio = median(cruces_medio_triceps);

promedio_cruces_triceps_tres = mean(cruces_tres_triceps);
mediana_cruces_triceps_tres = median(cruces_tres_triceps);

% Mostrar resultados estadísticos en la consola
disp('Resultados Estadísticos:');
disp('RMS del Bíceps con 1/2 KG:');
disp(['Promedio: ', num2str(promedio_rms_biceps_medio)]);
disp(['Mediana: ', num2str(mediana_rms_biceps_medio)]);

disp('RMS del Bíceps con 3 KG:');
disp(['Promedio: ', num2str(promedio_rms_biceps_tres)]);
disp(['Mediana: ', num2str(mediana_rms_biceps_tres)]);

disp('RMS del Tríceps con 1/2 KG:');
disp(['Promedio: ', num2str(promedio_rms_triceps_medio)]);
disp(['Mediana: ', num2str(mediana_rms_triceps_medio)]);

disp('RMS del Tríceps con 3 KG:');
disp(['Promedio: ', num2str(promedio_rms_triceps_tres)]);
disp(['Mediana: ', num2str(mediana_rms_triceps_tres)]);

disp('Área bajo la curva (AUC) del Bíceps con 1/2 KG:');
disp(['Promedio: ', num2str(promedio_auc_biceps_medio)]);
disp(['Mediana: ', num2str(mediana_auc_biceps_medio)]);

disp('Área bajo la curva (AUC) del Bíceps con 3 KG:');
disp(['Promedio: ', num2str(promedio_auc_biceps_tres)]);
disp(['Mediana: ', num2str(mediana_auc_biceps_tres)]);

disp('Área bajo la curva (AUC) del Tríceps con 1/2 KG:');
disp(['Promedio: ', num2str(promedio_auc_triceps_medio)]);
disp(['Mediana: ', num2str(mediana_auc_triceps_medio)]);

disp('Área bajo la curva (AUC) del Tríceps con 3 KG:');
disp(['Promedio: ', num2str(promedio_auc_triceps_tres)]);
disp(['Mediana: ', num2str(mediana_auc_triceps_tres)]);

disp('Número de Cruces por Cero del Bíceps con 1/2 KG:');
disp(['Promedio: ', num2str(promedio_cruces_biceps_medio)]);
disp(['Mediana: ', num2str(mediana_cruces_biceps_medio)]);

disp('Número de Cruces por Cero del Bíceps con 3 KG:');
disp(['Promedio: ', num2str(promedio_cruces_biceps_tres)]);
disp(['Mediana: ', num2str(mediana_cruces_biceps_tres)]);

disp('Número de Cruces por Cero del Tríceps con 1/2 KG:');
disp(['Promedio: ', num2str(promedio_cruces_triceps_medio)]);
disp(['Mediana: ', num2str(mediana_cruces_triceps_medio)]);

disp('Número de Cruces por Cero del Tríceps con 3 KG:');
disp(['Promedio: ', num2str(promedio_cruces_triceps_tres)]);
disp(['Mediana: ', num2str(mediana_cruces_triceps_tres)]);

% Crear tabla para los resultados
resultados = table((1:5)', rms_medio_biceps', rms_tres_biceps', rms_medio_triceps', rms_tres_triceps', ...
    area_medio_biceps', area_tres_biceps', area_medio_triceps', area_tres_triceps', ...
    cruces_medio_biceps', cruces_tres_biceps', cruces_medio_triceps', cruces_tres_triceps', ...
    'VariableNames', {'Repeticion', 'RMS_Biceps_1_2KG', 'RMS_Biceps_3KG', 'RMS_Triceps_1_2KG', 'RMS_Triceps_3KG', ...
    'AUC_Biceps_1_2KG', 'AUC_Biceps_3KG', 'AUC_Triceps_1_2KG', 'AUC_Triceps_3KG', ...
    'Cruces_Biceps_1_2KG', 'Cruces_Biceps_3KG', 'Cruces_Triceps_1_2KG', 'Cruces_Triceps_3KG'});

disp(resultados);

% Crear tabla para valores promedio y medianos
resultados_estadisticos = table({'Promedio'; 'Mediana'}, ...
    [promedio_rms_biceps_medio; mediana_rms_biceps_medio], ...
    [promedio_rms_biceps_tres; mediana_rms_biceps_tres], ...
    [promedio_rms_triceps_medio; mediana_rms_triceps_medio], ...
    [promedio_rms_triceps_tres; mediana_rms_triceps_tres], ...
    [promedio_auc_biceps_medio; mediana_auc_biceps_medio], ...
    [promedio_auc_biceps_tres; mediana_auc_biceps_tres], ...
    [promedio_auc_triceps_medio; mediana_auc_triceps_medio], ...
    [promedio_auc_triceps_tres; mediana_auc_triceps_tres], ...
    [promedio_cruces_biceps_medio; mediana_cruces_biceps_medio], ...
    [promedio_cruces_biceps_tres; mediana_cruces_biceps_tres], ...
    [promedio_cruces_triceps_medio; mediana_cruces_triceps_medio], ...
    [promedio_cruces_triceps_tres; mediana_cruces_triceps_tres], ...
    'VariableNames', {'Medida', 'RMS_Biceps_1_2KG', 'RMS_Biceps_3KG', ...
    'RMS_Triceps_1_2KG', 'RMS_Triceps_3KG', 'AUC_Biceps_1_2KG', 'AUC_Biceps_3KG', ...
    'AUC_Triceps_1_2KG', 'AUC_Triceps_3KG', 'Cruces_Biceps_1_2KG', ...
    'Cruces_Biceps_3KG', 'Cruces_Triceps_1_2KG', 'Cruces_Triceps_3KG'});

disp(resultados_estadisticos);


%% Punto 5

segmento = 3;
% Extraer los segmentos de interés
senal_medio = b_medio(segmentos_medio(segmento,1):segmentos_medio(segmento,2));
senal_tres = b_tres(segmentos_tres(segmento,1):segmentos_tres(segmento,2));

% Calcular el espectro 
nfft = 2^nextpow2(length(senal_medio)); 
[Pxx_medio, F_medio] = pwelch(senal_medio, [], [], nfft, fs);
[Pxx_tres, F_tres] = pwelch(senal_tres, [], [], nfft, fs);

% Limitar el rango de frecuencia a 0-600 Hz
idx_600_medio = find(F_medio <= 600, 1, 'last');
idx_600_tres = find(F_tres <= 600, 1, 'last');

% Calcular frecuencia media y mediana
freq_media_medio = sum(F_medio(1:idx_600_medio) .* Pxx_medio(1:idx_600_medio)) / sum(Pxx_medio(1:idx_600_medio));
freq_mediana_medio = F_medio(find(cumsum(Pxx_medio(1:idx_600_medio)) >= sum(Pxx_medio(1:idx_600_medio))/2, 1));
freq_media_tres = sum(F_tres(1:idx_600_tres) .* Pxx_tres(1:idx_600_tres)) / sum(Pxx_tres(1:idx_600_tres));
freq_mediana_tres = F_tres(find(cumsum(Pxx_tres(1:idx_600_tres)) >= sum(Pxx_tres(1:idx_600_tres))/2, 1));

% Calcular potencia total
potencia_total_medio = sum(Pxx_medio(1:idx_600_medio));
potencia_total_tres = sum(Pxx_tres(1:idx_600_tres));

% Graficar
figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
graficarEspectro(F_medio(1:idx_600_medio), Pxx_medio(1:idx_600_medio), freq_media_medio, freq_mediana_medio, potencia_total_medio, 'Espectro de Frecuencia Bíceps con Medio Kilo: Segmento 3');

subplot(1, 2, 2);
graficarEspectro(F_tres(1:idx_600_tres), Pxx_tres(1:idx_600_tres), freq_media_tres, freq_mediana_tres, potencia_total_tres, 'Espectro de Frecuencia Bíceps con Tres Kilos: Segmento 3');

sgtitle('Espectro de Frecuencia para EMG del Bíceps: Segmento 3', 'FontSize', 16);

% Mostrar resultados
fprintf('Resultados para Bíceps con Medio Kilo del Segmento 3:\n');
fprintf('Frecuencia Media: %.2f Hz\n', freq_media_medio);
fprintf('Frecuencia Mediana: %.2f Hz\n', freq_mediana_medio);
fprintf('Potencia Total: %.2e\n\n', potencia_total_medio);

fprintf('Resultados para Bíceps con Tres Kilos del Segmento 3:\n');
fprintf('Frecuencia Media: %.2f Hz\n', freq_media_tres);
fprintf('Frecuencia Mediana: %.2f Hz\n', freq_mediana_tres);
fprintf('Potencia Total: %.2e\n\n', potencia_total_tres);

% Función para graficar y marcar las áreas de F. media y mediana
function graficarEspectro(F, Pxx, freq_media, freq_mediana, potencia_total, titulo)
    plot(F, Pxx, 'b', 'LineWidth', 1.5);
    hold on;

    idx_media = find(F <= freq_media, 1, 'last');
    idx_mediana = find(F <= freq_mediana, 1, 'last');

    fill([F(1:idx_mediana); flipud(F(1:idx_mediana))], [Pxx(1:idx_mediana); zeros(size(Pxx(1:idx_mediana)))], 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    fill([F(idx_mediana:idx_media); flipud(F(idx_mediana:idx_media))], [Pxx(idx_mediana:idx_media); zeros(size(Pxx(idx_mediana:idx_media)))], 'm', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    xline(freq_media, 'm--', 'LineWidth', 1.5);
    xline(freq_mediana, 'c--', 'LineWidth', 1.5);
    
    title(titulo, 'FontSize', 14);
    xlabel('Frecuencia [Hz]', 'FontSize', 12);
    ylabel('Potencia [V²]', 'FontSize', 12);
    xlim([0 600]); 
    ylim([0 max(Pxx)*1.1]);

    text(550, max(Pxx)*0.5, sprintf('Potencia Total: %.2e', potencia_total), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    
    legend('Espectro','Frecuencia Media', 'Frecuencia Mediana', 'Location', 'northeast');
    grid on;
end



