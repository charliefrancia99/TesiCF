% Importa i dati
data = importLandmarksData_copia('.csv');

% Estrai le colonne dei keypoints destro
shoulderX = data.LEFT_SHOULDERX;
shoulderY = data.LEFT_SHOULDERY;
shoulderZ = data.LEFT_SHOULDERZ;

elbowX = data.LEFT_ELBOWX;
elbowY = data.LEFT_ELBOWY;
elbowZ = data.LEFT_ELBOWZ;

wristX = data.LEFT_WRISTX;
wristY = data.LEFT_WRISTY;
wristZ = data.LEFT_WRISTZ;

% Calcola l'angolo di abduzione-adduzione della spalla
numFrames = length(shoulderX);

angoloSpallaSinistro = zeros(numFrames, 1);

for i = 1:numFrames
    % Segmento spalla-gomito
    shoulderElbowVector = [elbowX(i) - shoulderX(i), elbowY(i) - shoulderY(i), elbowZ(i) - shoulderZ(i)];

    % Punto ottenuto dalle coordinate specificate
    customPoint = [shoulderX(i), shoulderY(i) + 20, shoulderZ(i)];

    % Segmento spalla-punto
    shoulderCustomPointVector = [customPoint(1) - shoulderX(i), customPoint(2) - shoulderY(i), customPoint(3) - shoulderZ(i)];

    % Calcolo dell'angolo di abduzione-adduzione della spalla
    dotProduct = dot(shoulderCustomPointVector, shoulderElbowVector);
    normShoulderCustomPoint = norm(shoulderCustomPointVector);
    normShoulderElbow = norm(shoulderElbowVector);

    angoloSpallaSinistro(i) = acosd(dotProduct / (normShoulderCustomPoint * normShoulderElbow));
end

% Plot dell'angolo della spalla rispetto ai frame
frames = data.Frame;
time1= frames/60;
figure
plot(time1, angoloSpallaSinistro, 'b-', 'LineWidth', 2);

% Personalizza il grafico
title('Angolo di abduzione-adduzione della spalla sinistra');
xlabel('Time(s)');
ylabel('Angolo (gradi)');
grid on;

%% FILTRAGGIO GO PRO

% Frequenza di campionamento
fs = 60; % 60 FPS

% Frequenza di taglio desiderata (4 Hz)
fc = 4; % Frequenza di taglio in Hz

% Calcola le frequenze normalizzate
Wn = fc / (fs / 2);

% Progetta il filtro Butterworth
order = 4; % Ordine del filtro
[b, a] = butter(order, Wn);

% Filtra i dati dell'angolo destro del gomito
angoloSpallaSinistroFiltrato = filtfilt(b, a, angoloSpallaSinistro);

% Plot dell'angolo destro filtrato rispetto ai frame
figure
plot(time1, angoloSpallaSinistroFiltrato, 'r-', 'LineWidth', 2);

% Personalizza il grafico
title('Angolo sinistro spalla (Filtrato)');
xlabel('Time(s)');
ylabel('Angolo (gradi)');
grid on;


%% XSENS

% Specifica il percorso del tuo file Excel
file_path = '.xlsx';

% Carica il file Excel
xls_data = xlsread(file_path, 'Joint Angles ZXY');

% Estrai le colonne di interesse
frame = xls_data(:, 1);
time2=frame/60;
left_shoulder_abduction_adduction = xls_data(:, 35);

% Crea il grafico
figure;
plot(time2, left_shoulder_abduction_adduction, 'b.-');
xlabel('Time (s)');
ylabel('Left shoulder abduction adduction');
title('Grafico left shoulder abduction adduction vs. Time');
grid on;
%% SINCRONIZZAZIONE

% Calcola la differenza nella lunghezza tra i due vettori tempo
lunghezza_diff = length(time1) - length(time2);

% Verifico se il secondo vettore tempo è più corto del primo
if lunghezza_diff < 0
    % Rimuovi gli elementi iniziali in eccesso dal secondo angolo
    leftt_shoulder_abduction_adduction = left_shoulder_abduction_adduction(-lunghezza_diff : end-1);
    time2 = time2(1:length(left_shoulder_abduction_adduction)); % Aggiorna il secondo vettore tempo
elseif lunghezza_diff > 0
    % Rimuovi gli elementi iniziali in eccesso dal primo angolo
    angoloSpallaSinistroFiltrato = angoloSpallaSinistroFiltrato(lunghezza_diff : end-1);
    time1 = time1(1:length(angoloSpallaSinistroFiltrato)); % Aggiorna il primo vettore tempo
end

% Ora i due vettori angolo e i due vettori tempo hanno la stessa lunghezza
% Crea un nuovo figure con un subplot
figure;

% Sottoplot per il primo angolo
subplot(2, 1, 1);
plot(time1, angoloSpallaSinistroFiltrato, 'b-', 'LineWidth', 2);
xlabel('Time(s)');
ylabel('Angolo sinistro spalla (gradi)');
title('Angolo sinistro spalla vs. Time');
grid on;

hold on
% Sottoplot per il secondo angolo
%subplot(2, 1, 2);
plot(time2, left_shoulder_abduction_adduction, 'r-', 'LineWidth', 2);
xlabel('Time(s)');
ylabel('Angolo sinistro spalla (gradi)');
title('Angolo sinistro spalla vs. Time');
grid on;
i=length(time1);

% Regola la spaziatura tra i sottoplot
%spacing = 0.05;
%set(subplot(2, 1, 1), 'Position', [0.1, 0.55, 0.8, 0.4]);
%set(subplot(2, 1, 2), 'Position', [0.1, 0.1, 0.8, 0.4]);


%% SINCRONIZZAZIONE 2 se piu lungo right knee
% Visualizza i grafici originali
figure
subplot(2, 1, 1);
plot(time1, angoloSpallaSinistroFiltrato, 'b-', 'LineWidth', 2);
title('angoloSpallaSinistroFiltrato');
xlabel('Time(s)');
ylabel('Angolo (gradi)');
grid on;

subplot(2, 1, 2);
plot(time2, left_shoulder_abduction_adduction, 'b.-');
title('left_shoulder_abduction_adduction');
xlabel('Time (s)');
ylabel('left_shoulder_abduction_adduction');
grid on;

% Scegli manualmente un punto su angoloDestroGinocchioFiltrato
disp('Seleziona un punto su angoloSpallaSinistroFiltrato');
[x1, ~] = ginput(1);

% Scegli manualmente un punto su right_knee_flexion_extension
disp('Seleziona un punto su left_hip_flexion_extension.');
[x2, ~] = ginput(1);

% Calcola la differenza temporale
temporalDifference = x2 - x1;

% Interpola angoloDestroGinocchioFiltrato per ottenere la sincronizzazione
time1_sync = time1 + temporalDifference;
angoloSpallaSinistroFiltratoSync = interp1(time1_sync, angoloSpallaSinistroFiltrato, time1, 'linear', 'extrap');
angoloSpallaSinistroFiltrato=angoloSpallaSinistroFiltratoSync;
left_shoulder_abduction_adduction=left_shoulder_abduction_adduction(1:length(angoloSpallaSinistroFiltrato));
% Aggiorna il grafico sincronizzato sovrapposto
figure
plot(time1, angoloSpallaSinistroFiltrato, 'r-', 'LineWidth', 2);
hold on;
plot(time1, left_shoulder_abduction_adduction, 'b.-');
hold off;

title('Grafici Sovrapposti - Sincronizzati');
xlabel('Time(s)');
ylabel('Angolo (gradi)');
legend('Angolo sinistro spalla (Filtrato)', 'left_shoulder_abduction_adduction');
grid on;
time2=time1;

%%
% Grafico finale dopo la sincronizzazione
figName = figure;
MATLABcol = colororder;

% Plot della flessione/estensione continua del ginocchio destro (Mediapipe)
plot(time1, angoloSpallaSinistroFiltrato, 'Color', MATLABcol(1,:), 'LineWidth', 2);
hold on;

% Plot della flessione/estensione continua del ginocchio destro (Xsens)
plot(time2, left_shoulder_abduction_adduction, 'Color', MATLABcol(2,:), 'LineWidth', 2);

hold off;

title('Left Shoulder Abduction/Adduction Comparison','interpreter', 'latex');
xlabel('Time [s]','interpreter', 'latex');
ylabel('Angle [$^{\circ}$]','interpreter', 'latex');
legend('Mediapipe', 'Xsens', 'TextColor', [0 0 0],'interpreter', 'latex');

% Personalizzazione dell'aspetto
set(gca, 'FontSize', 12); % Imposta la dimensione del carattere dell'asse
set(gca, 'LineWidth', 1.5); % Imposta lo spessore delle linee dell'asse

% Personalizza i colori e lo spessore delle linee
set(gca, 'YColor', [0 0 0]); % Colore dell'asse y nero
set(gca, 'XColor', [0 0 0]); % Colore dell'asse x nero

set(gcf, 'Color', [1 1 1]); % Sfondo bianco del grafico
set(gca, 'Color', [1 1 1]); % Colore di sfondo dell'asse

% Aggiorna i colori e lo spessore delle linee delle curve
set(findobj(gca, 'Type', 'Line', 'Color', MATLABcol(1,:)), 'LineWidth', 2); % Linea Mediapipe
set(findobj(gca, 'Type', 'Line', 'Color', MATLABcol(2,:)), 'LineWidth', 2); % Linea Xsens

% Modifica l'aspetto dell'asse x
xlim([10, max(time1)]); % Imposta il limite inferiore dell'asse x a 10 secondi
set(gca,'TickLabelInterpreter','latex');
%ylim([0,100]);
%set(gca, 'XMinorTick', 'on'); % Abilita le tick minori sull'asse x
grid on

%%
exportgraphics(figName, "grafico_left_shoulder.pdf");

%% TEST di Jarque-Bera
% Esegui il test di Jarque-Bera su angoloDestroGinocchioFiltrato
[h1, p1, jbstat1, critval1] = jbtest(angoloSpallaSinistroFiltrato);

% Stampare i risultati del test per angoloDestroGinocchioFiltrato
fprintf('Risultati del test di Jarque-Bera per angoloSpallaSinistroFiltrato:\n');
fprintf('Statistiche JB: %.4f\n', jbstat1);
fprintf('Valore p: %.4f\n', p1);
if h1
    fprintf('I dati NON seguono una distribuzione normale.\n');
else
    fprintf('I dati seguono una distribuzione normale.\n');
end

% Esegui il test di Jarque-Bera su right_knee_flexion_extension
[h2, p2, jbstat2, critval2] = jbtest(left_shoulder_abduction_adduction);

% Stampare i risultati del test per right_knee_flexion_extension
fprintf('\nRisultati del test di Jarque-Bera per left_shoulder_abduction_adduction:\n');
fprintf('Statistiche JB: %.4f\n', jbstat2);
fprintf('Valore p: %.4f\n', p2);
if h2
    fprintf('I dati NON seguono una distribuzione normale.\n');
else
    fprintf('I dati seguono una distribuzione normale.\n');
end


%Test di Jarque-Bera:
%Il test di Jarque-Bera è un test statistico che valuta se un insieme di dati ha 
%una distribuzione simile a una distribuzione normale. Si basa su due 
%momenti statistici: la curtosi (misura la pesantezza delle code della 
%distribuzione rispetto a una distribuzione normale) e lo 
%skewness (misura l'asimmetria dei dati). Se i dati 
%provengono da una distribuzione normale, la curtosi sarà vicina a 3 e lo skewness sarà vicino a 0. Il test confronta queste due statistiche con la distribuzione normale, aspettandosi valori prossimi a 3 e 0 rispettivamente per dati normalmente distribuiti.

%Nel test di Jarque-Bera, vengono restituiti due valori:

%Statistiche del test (JB): Questo valore rappresenta il risultato del calcolo statistico del test di Jarque-Bera. Indica la distanza dei momenti (skewness e curtosi) dei dati dalla distribuzione normale. Maggiore è la deviazione, maggiore è la discrepanza rispetto alla normalità.
%Valore p: Questo valore indica la significatività statistica del test. Se è inferiore al livello di significatività (comunemente 0.05), si può rifiutare l'ipotesi che i dati seguano una distribuzione normale. In altre parole:
%Se il valore p è basso, si ha evidenza a favore dell'ipotesi che i dati non siano normalmente distribuiti.
%Se il valore p è alto, non c'è sufficiente evidenza per affermare che i dati non siano normali.
%Per interpretare i risultati, generalmente si considera il valore p: se è basso, si può concludere che i dati non seguono una distribuzione normale.


%Test di Shapiro-Wilk:
%Il test di Shapiro-Wilk valuta se un insieme di dati proviene da una distribuzione normale. Si basa sulla differenza tra i valori osservati e quelli attesi in una distribuzione normale, analizzando come questa discrepanza si distribuisce. Il test confronta i valori osservati con quelli attesi se i dati provenissero da una distribuzione normal


%% Calcolo del coefficiente di Spearman
rho = corr(angoloSpallaSinistroFiltrato, left_shoulder_abduction_adduction, 'Type', 'Spearman');
fprintf('Il coefficiente di correlazione di Spearman tra i due segnali è: %.4f\n', rho);

%% ICC NUOVO
% Calcola l'ICC
[r, LB, UB, F, df1, df2, p] = ICC([angoloSpallaSinistroFiltrato, left_shoulder_abduction_adduction], 'A-1', 0.05, 0);

% Stampa i risultati
fprintf('ICC tra angoloSpallaSinistroFiltrato e left_shoulder_abduction_adduction:\n');
fprintf('Coefficient: %.4f\n', r);
%fprintf('Lower Bound: %.4f\n', LB);
%fprintf('Upper Bound: %.4f\n', UB);
%fprintf('F-Value: %.4f\n', F);
%fprintf('Degrees of Freedom 1: %d\n', df1);
%fprintf('Degrees of Freedom 2: %d\n', df2);
fprintf('P-Value: %.4f\n', p);


% ICC vecchio
% I tuoi dati
%data1 = angoloDestroGinocchioFiltrato;
%data2 = right_knee_flexion_extension;

% Calcolo della media
%meanData = mean([data1, data2], 2);

% Calcolo delle somme dei quadrati
%SSB = sum((meanData - data1).^2) + sum((meanData - data2).^2);
%SSW = sum((data1 - meanData).^2) + sum((data2 - meanData).^2);

% Calcolo dell'ICC
%n = length(data1);
%k = 2; % Numero di osservazioni
%MSB = SSB / (k - 1);
%MSW = SSW / ((n * k) - k);

%ICC = (MSB - MSW) / (MSB + (k - 1) * MSW);
%fprintf('Il valore dell''Intraclass Correlation Coefficient (ICC) è: %.4f\n', ICC);

%% Bland Altman 
%se vuoi togliere qualcosa dal grafico basta che vai nella parte di Display the stats of interest on the BA plot
%e nello switch case commenti quello che non vuoi visualizzare

% Calcolo del bias
bias = mean(left_shoulder_abduction_adduction - angoloSpallaSinistroFiltrato);
fprintf('Il bias tra i due segnali è: %.4f\n', bias);

% Calcolo della deviazione standard
std_dev = std(left_shoulder_abduction_adduction - angoloSpallaSinistroFiltrato);
fprintf('La deviazione standard per la differenza dei  due segnali è: %.4f\n', std_dev);

%CV
% Calcolo del Coefficiente di Variazione (CV)
cv = (std_dev / mean((angoloSpallaSinistroFiltrato + left_shoulder_abduction_adduction)/2)) * 100;
%alternativa cv = (std_dev / mean((angoloDestroGinocchioFiltrato;left_shoulder_abduction_adduction)/2)) * 100;
%tanto li mette tutti in fila è la stessa cosa
fprintf('Il Coefficiente di Variazione tra i due segnali è: %.4f%%\n', cv);

% Bland-Altman plot
BlandAltman(angoloSpallaSinistroFiltrato, left_shoulder_abduction_adduction, {'Angolo destro ginocchio', 'Angolo ginocchio'});
title('Bland-Altman Plot tra Angolo destro anca e Angolo ginocchio');

% Visualizza bias e deviazione standard nel grafico
legend(sprintf('Bias: %.4f', bias), sprintf('SD: %.4f', std_dev));

%%

% Creazione di un nuovo Bland-Altman plot solo con il bias
figName = figure;
MATLABcol = colororder;

% Imposta il colore di sfondo del grafico su bianco
set(gcf, 'Color', 'w');

% Calcola la differenza e la media
differenza = left_shoulder_abduction_adduction - angoloSpallaSinistroFiltrato;
media = mean([angoloSpallaSinistroFiltrato, left_shoulder_abduction_adduction], 2);

% Calcola il bias
bias = mean(differenza);

% Calcola i limiti di accordo
std_dev = std(differenza);
limits_1_96 = bias + 1.96 * std_dev;
limits_minus_1_96 = bias - 1.96 * std_dev;

% Plot della differenza come punti (pallini blu scuro)
plot(media, differenza, 'o', 'Color', MATLABcol(1,:), 'DisplayName', 'Differenza');
hold on;

% Plot del bias come linea continua orizzontale (linea continua rossa)
bias_line = plot([min(media), max(media)], [bias, bias], '-r', 'LineWidth', 2, 'DisplayName', ' ');

% Plot dei limiti di accordo come linee orizzontali tratteggiate (linee tratteggiate nere)
limits_line_1 = plot([min(media), max(media)], [limits_1_96, limits_1_96], '--k', 'DisplayName', ' ');
limits_line_2 = plot([min(media), max(media)], [limits_minus_1_96, limits_minus_1_96], '--k');

% Aggiungi etichetta con il valore del bias a destra della linea
bias_label = sprintf('%.4f', bias);
text(max(media) + 0.5, bias, bias_label, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', 'red', 'interpreter', 'latex');

% Aggiungi etichetta con il valore della deviazione standard in alto a sinistra
sd_label = sprintf('SD: %.4f', std_dev);
text(min(media), max(differenza)+4.8, sd_label, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'black', 'interpreter', 'latex' );

% Aggiungi etichette ai limiti di accordo a destra
limits_label_1_96 = sprintf(' %.4f', limits_1_96);
limits_label_minus_1_96 = sprintf(' %.4f', limits_minus_1_96);

text(max(media) + 0.5, limits_1_96, limits_label_1_96, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', 'black', 'interpreter', 'latex');
text(max(media) + 0.5, limits_minus_1_96, limits_label_minus_1_96, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', 'black', 'interpreter', 'latex');

% Legenda
legend([bias_line, limits_line_1], {'Bias', sprintf('Limits of agreement (Bias $\\pm$ 1.96SD)')}, 'Location', 'Best', 'interpreter', 'latex', 'interpreter', 'latex');

title('Bland-Altman Plot Left Shoulder Abduction/Adduction', 'interpreter', 'latex');
xlabel('Mean Xsens \& Mediapipe', 'interpreter', 'latex');
ylabel('Difference Xsens-Mediapipe', 'interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
ylim([-15, 15]);
hold off;
box off;
%% RMSE giusto

kinect_segment = angoloSpallaSinistroFiltrato(1:i);

opto_segment = left_shoulder_abduction_adduction(1:i);
 
% Calcola la differenza tra i due segnali

diff = opto_segment - kinect_segment;

 % Calcola il quadrato di ciascun elemento della differenza

diff_squared = diff.^2;

 % Calcola la somma degli elementi quadrati

sum_squared = sum(diff_squared);
 
% Calcola la media degli elementi quadrati

mean_squared = sum_squared / numel(diff_squared);

% Calcola la radice quadrata della media degli elementi quadrati per ottenere l'RMSE

rmse = sqrt(mean_squared)

%% ROM GO PRO

% Calcola il ROM per "angoloDestroGinocchioFiltrato
rom_angoloSinistro = max(angoloSpallaSinistroFiltrato) - min(angoloSpallaSinistroFiltrato);
fprintf('Range of Motion per angoloSpallaSinistroFiltrato: %.2f gradi\n', rom_angoloSinistro);

%% ROM XSENS
% Calcola il ROM per "right_knee_flexion_extension"
rom_left_elbow = max(left_shoulder_abduction_adduction) - min(left_shoulder_abduction_adduction);
fprintf('Range of Motion per left_shoulder_abduction_adduction: %.2f gradi\n', rom_left_elbow);

%% RMSE E ROM CALCOLATO NELLE DIVERSE FINESTRE TEMPORALI

% Disegna il grafico dell'angolo destro del ginocchio e dei dati ottici
figure;
plot(time1, angoloSpallaSinistroFiltrato, 'b-', 'LineWidth', 2);
hold on;
plot(time2, left_shoulder_abduction_adduction, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Angle (degrees)');
title('Angle vs. Time');
grid on;

% Chiedi all'utente di selezionare manualmente i punti iniziali e finali delle finestre
fprintf('Clicca sul grafico per selezionare i punti iniziali e finali delle finestre.\n');
fprintf('Premi "Enter" quando hai finito.\n');

% Inizializza un vettore per memorizzare gli RMSE e i ROM delle finestre
rmse_values = zeros(4, 1);
rom_values = zeros(4, 2); % Un array 2x4 per memorizzare il ROM per entrambi i segnali in ciascuna finestra

% Inizializza i vettori per memorizzare le deviazioni standard
dev_std_angolo = zeros(3, 1);
dev_std_opto = zeros(3, 1);
dev_std_diff = zeros(3, 1);

for window = 1:4
    % Utilizza ginput per consentire all'utente di selezionare i punti sul grafico
    [x, ~] = ginput(2); % L'utente seleziona 2 punti per ciascuna finestra
    
    % Trova gli indici dei punti nel vettore tempo
    window_start_index = find(time1 >= x(1), 1);
    window_end_index = find(time1 <= x(2), 1, 'last');
    
    % Seleziona gli angoli all'interno della finestra corrente
    angolo_window = angoloSpallaSinistroFiltrato(window_start_index:window_end_index);
    
    % Seleziona i dati ottici all'interno della finestra corrente
    opto_window = left_shoulder_abduction_adduction(window_start_index:window_end_index);
    
    % Calcola la differenza tra i due segnali
    diff = opto_window - angolo_window;
    
    % Calcola il quadrato di ciascun elemento della differenza
    diff_squared = diff.^2;
    
    % Calcola la somma degli elementi quadrati
    sum_squared = sum(diff_squared);
    
    % Calcola la media degli elementi quadrati
    mean_squared = sum_squared / numel(diff_squared);
    
    % Calcola la radice quadrata della media degli elementi quadrati per ottenere l'RMSE
    rmse = sqrt(mean_squared);
    
    % Memorizza l'RMSE nel vettore
    rmse_values(window) = rmse;
    
    % Calcola il ROM per entrambi i segnali nella finestra corrente
    rom_angolo = max(angolo_window) - min(angolo_window);
    rom_opto = max(opto_window) - min(opto_window);
    
    % Memorizza il ROM nel vettore
    rom_values(window, 1) = rom_angolo;
    rom_values(window, 2) = rom_opto;
     % Seleziona i dati delle finestre 2, 3 e 4
    if window >= 2
        angolo_finestre_selezionate = angoloSpallaSinistroFiltrato(window_start_index:window_end_index);
        opto_finestre_selezionate = left_shoulder_abduction_adduction(window_start_index:window_end_index);

        % Calcola la deviazione standard per entrambi i segnali nelle finestre 2, 3 e 4
        dev_std_angolo(window-1) = std(angolo_finestre_selezionate);
        dev_std_opto(window-1) = std(opto_finestre_selezionate);

        % Calcola la differenza tra i due segnali nelle finestre 2, 3 e 4
        diff_finestre_selezionate = opto_finestre_selezionate - angolo_finestre_selezionate;

        % Calcola la deviazione standard sulla differenza per le finestre 2, 3 e 4
        dev_std_diff(window-1) = std(diff_finestre_selezionate);
    end
end

% Visualizza i risultati
fprintf('RMSE per le finestre selezionate:\n');
for window = 1:4
    fprintf('Finestra %d: %.2f gradi\n', window, rmse_values(window));
end

fprintf('ROM per angoloSpallaSinistroFiltrato nelle finestre selezionate:\n');
for window = 1:4
    fprintf('Finestra %d: %.2f gradi\n', window, rom_values(window, 1));
end

fprintf('ROM per left_shoulder_abduction_adduction nelle finestre selezionate:\n');
for window = 1:4
    fprintf('Finestra %d: %.2f gradi\n', window, rom_values(window, 2));
end

fprintf('Deviazione standard per angoloSinistroAncaFiltrato nelle finestre 2, 3 e 4: %.2f gradi\n', mean(dev_std_angolo));
fprintf('Deviazione standard per left_hip_flexion_extension nelle finestre 2, 3 e 4: %.2f gradi\n', mean(dev_std_opto));
fprintf('Deviazione standard sulla differenza tra opto e angolo nelle finestre 2, 3 e 4: %.2f gradi\n', mean(dev_std_diff));

