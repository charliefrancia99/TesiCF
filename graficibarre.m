%prima altermnativa
% Definizione dei dati
condizioni = {'Right KFE', 'Left KFE', 'Right HFE', 'Left HFE'};
% Definizione dei dati
condizioni = {'Right KFE', 'Left KFE', 'Right HFE', 'Left HFE'};
RMSE = [15.23, 15.76, 14.87, 15.75]; % Sostituisci con i tuoi valori
SD = [2.61, 2.62, 2.81, 2.75]; % Sostituisci con i tuoi valori
percent_RMSE = [13.33, 13.57, 14.87, 15.45]; % Sostituisci con i tuoi valori

% Creazione del grafico a barre
figure;
bar(RMSE);
hold on;

% Plottaggio delle deviazioni standard come errori sull'altezza delle barre
errorbar(1:numel(RMSE), RMSE, SD, 'k.', 'LineWidth', 1);

% Aggiunta delle etichette agli assi e al titolo
set(gca, 'XTick', 1:numel(condizioni), 'XTickLabel', condizioni);
ylabel('RMSE');
title('Valori RMSE, Deviazioni standard e Percentuale RMSE per le condizioni');

% Aggiunta della legenda
legend('RMSE', 'Deviazione Standard', 'Location', 'Best');

% Aggiunta della percentuale RMSE sopra ogni barra
for i = 1:numel(RMSE)
    text(i, RMSE(i), sprintf('%.2f%%', percent_RMSE(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

hold off;

%% seconda alternatuva

% Definizione dei dati
condizioni = {'Right KFE', 'Left KFE', 'Right HFE', 'Left HFE'};
RMSE = [15.23, 15.76, 14.87, 15.75]; % Sostituisci con i tuoi valori di RMSE
percent_RMSE = [13.33, 13.57, 14.87, 15.45]; % Sostituisci con i tuoi valori di percentuale di RMSE
SD_RMSE = [2.61, 2.62, 2.81, 2.75]; % Sostituisci con i tuoi valori di deviazione standard per RMSE

% Creazione dei colori
colore_verde = [0.2, 0.6, 0.2]; % Verde non troppo acceso
colore_arancione = [1, 0.5, 0]; % Arancione non troppo acceso

% Creazione del grafico a barre
figure;
bar_width = 0.4;
x = 1:numel(condizioni);

% Barre per RMSE con asse y sinistro
yyaxis left;
h1 = bar(x - bar_width/2, RMSE, bar_width, 'FaceColor', colore_verde);
hold on;

% Plottaggio delle deviazioni standard come errori solo sulla prima barra di ogni condizione
for i = 1:numel(condizioni)
    if ~isnan(SD_RMSE(i))
        errorbar(x(i) - bar_width/2, RMSE(i), SD_RMSE(i), 'k.', 'LineWidth', 1);
    end
end

% Barre per percentuale di RMSE con asse y destro
yyaxis right;
h2 = bar(x + bar_width/2, percent_RMSE, bar_width, 'FaceColor', colore_arancione);

% Aggiunta delle etichette agli assi e al titolo
set(gca, 'XTick', 1:numel(condizioni));
% Utilizzo del formato LaTeX per il testo delle condizioni
set(gca, 'XTickLabel', {'$$\textbf{Right KFE}$$', '$$\textbf{Left KFE}$$', '$$\textbf{Right HFE}$$', '$$\textbf{Left HFE}$$'}, 'TickLabelInterpreter', 'latex');
ylabel('RMSE', 'interpreter', 'latex');
yyaxis left;
ylim([0 25]);
ylabel('Valori di RMSE', 'interpreter', 'latex');
yyaxis right;
ylabel('Percentuale di RMSE', 'interpreter', 'latex');
ylim([0 25]);

% Modifica del colore del testo degli assi y
set(gca, 'YColor', colore_verde); % Per l'asse sinistro
yyaxis right;
set(gca, 'YColor', colore_arancione); % Per l'asse destro

% Aggiunta della legenda in alto a destra
legend([h1, h2], {'RMSE', 'Percentuale RMSE'}, 'Location', 'NorthEast', 'interpreter', 'latex');

hold off;



