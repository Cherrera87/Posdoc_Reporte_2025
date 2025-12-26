function F1 = plot_cobertura_temporal_mediciones(path_input, path_figure)
%% PLOT_COBERTURA_TEMPORAL_MEDICIONES
%
% Construye una figura de disponibilidad temporal de las mediciones
% de radar HF y ADCPs utilizados en el estudio.
%
% ENTRADAS:
%   path_input  : ruta donde se encuentran los archivos *.txt
%   path_figure : ruta donde se guardarán las figuras
%
% SALIDA:
%   F1 : handle de la figura generada
%
% Autor: Carlos F. Herrera Vázquez
% Año: 2025
% Institución: CICESE
% ======================================================================

%% ===================== VALIDACIÓN DE ENTRADAS ===================== %%
if nargin < 2
    error('Se requieren dos entradas: path_input y path_figure');
end

if ~isfolder(path_input)
    error('El directorio path_input no existe');
end

if ~isfolder(path_figure)
    mkdir(path_figure)
end

fecha = cell(4,1);

%% ===================== RADAR HF (WERA) ===================== %%
file_radar = fullfile(path_input,'Fechas_radar.txt');
tiempo_radar = readmatrix(file_radar);
fecha{1} = datetime(num2str(tiempo_radar), 'InputFormat','yyyyDDD');

%% ===================== ADCP PUNTA MORRO (PM) ===================== %%
% Archivo 1
file_PM1 = fullfile(path_input,'fechas_PM_signature.txt');
tiempo = readtable(file_PM1);

anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, ...
    'UniformOutput', false));

t0 = datetime([anio tiempo.Var6 tiempo.Var7 fix(tiempo.Var8/100) fix(tiempo.Var8/100)-(tiempo.Var8/100) tiempo.Var8*0], ...
    'Format','yyyy-MM-dd HH:mm:ss');

% Archivo 2
file_PM2 = fullfile(path_input,'fechas_PM_signature_2.txt');
tiempo = readtable(file_PM2);

anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, ...
    'UniformOutput', false));

t1 = datetime([anio tiempo.Var6 tiempo.Var7 fix(tiempo.Var8/100) fix(tiempo.Var8/100)-(tiempo.Var8/100) tiempo.Var8*0], ...
    'Format','yyyy-MM-dd HH:mm:ss');

fecha{2} = [t0; t1];

%% ===================== ADCP BAJO SAN MIGUEL (BSM) ===================== %%
% Archivo 1
file_BSM1 = fullfile(path_input,'fechas_BSM_sig_fil.txt');
tiempo = readtable(file_BSM1);

anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, ...
    'UniformOutput', false));

t0 = datetime([anio tiempo.Var6 tiempo.Var7 fix(tiempo.Var8/100) fix(tiempo.Var8/100)-(tiempo.Var8/100) tiempo.Var8*0], ...
    'Format','yyyy-MM-dd HH:mm:ss');

% Archivo 2
file_BSM2 = fullfile(path_input,'fechas_BSM_signature_2.txt');
tiempo = readtable(file_BSM2);

anio = str2double(cellfun(@(x) x(end-3:end), tiempo.Var5, ...
    'UniformOutput', false));

t1 = datetime([anio tiempo.Var6 tiempo.Var7 fix(tiempo.Var8/100) fix(tiempo.Var8/100)-(tiempo.Var8/100) tiempo.Var8*0], ...
    'Format','yyyy-MM-dd HH:mm:ss');

fecha{3} = [t0; t1];

%% ===================== ADCP ISLA TODOS SANTOS (ITS) ===================== %%
file_ITS1 = fullfile(path_input,'dias_ITS_1.txt');
file_ITS2 = fullfile(path_input,'dias_ITS_2.txt');
file_ITS3 = fullfile(path_input,'dias_ITS_3.txt');

time0 = dlmread(file_ITS1); time0(time0==0) = [];
time1 = dlmread(file_ITS2); time1(time1==0) = [];
time2 = dlmread(file_ITS3); time2(time2==0) = NaN;

fecha{4} = datetime(datevec([time2(:); time0(:); time1(:)]));

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

plot(fecha{1}, ones(numel(fecha{1}),1)*1,'.k','MarkerSize',6); hold on
plot(fecha{2}, ones(numel(fecha{2}),1)*2,'.r','MarkerSize',6);
plot(fecha{3}, ones(numel(fecha{3}),1)*3,'.b','MarkerSize',6);
plot(fecha{4}, ones(numel(fecha{4}),1)*4,'.g','MarkerSize',6);

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
nameFig = fullfile(path_figure,'serie_temporal_disponibilidad');

print(F1,'-depsc2',nameFig)
print(F1,'-dpng',  nameFig)
print(F1,'-dpdf',  nameFig)

end

