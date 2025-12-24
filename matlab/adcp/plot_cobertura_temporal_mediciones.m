function F1 = plot_cobertura_temporal_mediciones()
%% PLOT_COBERTURA_TEMPORAL_MEDICIONES
%
% Construye una figura de disponibilidad temporal de las mediciones
% utilizadas en el estudio, incluyendo:
%
%   1) Radar HF (WERA)
%   2) ADCP Nortek Signature en Punta Morro (PM)
%   3) ADCP Nortek Signature en Bajo San Miguel (BSM)
%   4) ADCP Aquadopp en Isla Todos Santos (ITS)
%
% La función lee archivos de texto con marcas de tiempo previamente
% extraídas de los registros instrumentales originales y genera una
% representación gráfica de los periodos de medición disponibles para
% cada sistema.
%
% IMPORTANTE:
% Esta función requiere la presencia de los siguientes archivos en el
% mismo directorio que este archivo .m:
%
%   - Fechas_radar.txt
%   - fechas_PM_signature.txt
%   - fechas_PM_signature_2.txt
%   - fechas_BSM_sig_fil.txt
%   - fechas_BSM_signature_2.txt
%   - dias_ITS_1.txt
%   - dias_ITS_2.txt
%   - dias_ITS_3.txt
%
% Salida:
%   F1 : handle de la figura generada
%
% La figura se exporta automáticamente en formatos EPS, PNG y PDF.
%
% Autor: Carlos F. Herrera Vázquez
% Año: 2025
% Institución: CICESE
%
% ======================================================================

close all; clc

%% ===================== CONFIGURACIÓN ===================== %%
% Usar el directorio donde se encuentra esta función
scriptPath = mfilename('fullpath');
[scriptDir,~,~] = fileparts(scriptPath);
cd(scriptDir);

fecha = cell(4,1);

%% ===================== RADAR HF (WERA) ===================== %%
% Archivo con columnas: YEAR  DOY
tiempo_radar = dlmread('Fechas_radar.txt');

YEAR = tiempo_radar(:,1);
DOY  = tiempo_radar(:,2);

fecha{1} = datetime(YEAR,1,1) + days(DOY-1);

%% ===================== ADCP PUNTA MORRO (PM) ===================== %%
% Archivo 1
tiempo = readtable('fechas_PM_signature.txt');
anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, 'UniformOutput', false));

t0 = datetime([anio tiempo.Var6 tiempo.Var7 tiempo.Var8/100], ...
              'Format','yyyy-MM-dd HH:mm:ss');

% Archivo 2
tiempo = readtable('fechas_PM_signature_2.txt');
anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, 'UniformOutput', false));

t1 = datetime([anio tiempo.Var6 tiempo.Var7 tiempo.Var8/100], ...
              'Format','yyyy-MM-dd HH:mm:ss');

fecha{2} = [t0; t1];

%% ===================== ADCP BAJO SAN MIGUEL (BSM) ===================== %%
% Archivo 1
tiempo = readtable('fechas_BSM_sig_fil.txt');
anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, 'UniformOutput', false));

t0 = datetime([anio tiempo.Var6 tiempo.Var7 tiempo.Var8/100], ...
              'Format','yyyy-MM-dd HH:mm:ss');

% Archivo 2
tiempo = readtable('fechas_BSM_signature_2.txt');
anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, 'UniformOutput', false));

t1 = datetime([anio tiempo.Var6 tiempo.Var7 tiempo.Var7*0], ...
              'Format','yyyy-MM-dd HH:mm:ss');

fecha{3} = [t0; t1];

%% ===================== ADCP ISLA TODOS SANTOS (ITS) ===================== %%
% Tiempos en segundos desde 1970-01-01
time0 = dlmread('dias_ITS_1.txt'); time0(time0==0) = [];
time1 = dlmread('dias_ITS_2.txt'); time1(time1==0) = [];
time2 = dlmread('dias_ITS_3.txt'); time2(time2==0) = [];

fecha{4} = datetime(datevec([time2; time0; time1]));

%% ===================== FIGURA ===================== %%
F1 = figure('Color','w');

anchoCarta = 11;
altoCarta  = 8.5;

set(F1,'Units','inches',...
       'Position',[1 1 anchoCarta altoCarta],...
       'PaperUnits','inches',...
       'PaperSize',[anchoCarta altoCarta],...
       'PaperOrientation','landscape',...
       'PaperPosition',[0 0 anchoCarta altoCarta]);

ax = subplot(4,1,[1 2]);

plot(fecha{1},ones(numel(fecha{1}),1)*1,'.k','MarkerSize',6); hold on
plot(fecha{2},ones(numel(fecha{2}),1)*2,'.r','MarkerSize',6);
plot(fecha{3},ones(numel(fecha{3}),1)*3,'.b','MarkerSize',6);
plot(fecha{4},ones(numel(fecha{4}),1)*4,'.g','MarkerSize',6);

ax.XGrid = 'on';
ax.YGrid = 'on';
ax.YLim  = [0 5];

ax.XTick = ax.XLim(1):calmonths(3):ax.XLim(end);
ax.XTickLabelRotation = 90;

ax.YTick = 1:4;
ax.YTickLabel = { ...
    '$\rm EPB_{WERA}^{HF\ Radar}$', ...
    '$\rm PM_{Signature}^{ADCP}$', ...
    '$\rm BSM_{Signature}^{ADCP}$', ...
    '$\rm ITS_{Aquadopp}^{ADCP}$'};
ax.TickLabelInterpreter = 'latex';

legend({'Radar HF','PM','BSM','ITS'},...
       'Interpreter','latex',...
       'FontSize',12,...
       'Location','southoutside',...
       'Orientation','horizontal');

%% ===================== EXPORTACIÓN ===================== %%
print(F1,'-depsc2','serie_temporal_disponibilidad')
print(F1,'-dpng','serie_temporal_disponibilidad')
print(F1,'-dpdf','serie_temporal_disponibilidad')

end

