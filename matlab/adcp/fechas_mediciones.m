close all; clear all; clc

% Obtiene la carpeta donde se encuentra este script
scriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(scriptPath);

% Cambia el directorio actual a la carpeta del script
cd(scriptDir);

%-------- RADAR
tiempo_radar = dlmread("Fechas_radar.txt");
a = [num2str(tiempo_radar)];

YEAR = str2num(a(:,1:end-3));
DOY = str2num(a(:,end-2:end));

fecha{1} = datetime(YEAR,1,1) + DOY-1;

%---------------------- PM
tiempo = readtable("fechas_PM_signature.txt");

a = tiempo.Var5;
tiempo.Var1 = str2num(cell2mat(cellfun(@(x) x(end-3:min(end, numel(x))), a, 'UniformOutput', false)));
t0 = datetime(double([tiempo.Var1 tiempo.Var6 ...
    tiempo.Var7 tiempo.Var8/100 tiempo.Var1*0 tiempo.Var1*0]));

tiempo = readtable("fechas_PM_signature_2.txt");
a = tiempo.Var5;
tiempo.Var1 = str2num(cell2mat(cellfun(@(x) x(end-3:min(end, numel(x))), a, 'UniformOutput', false)));

t1 = datetime(double([tiempo.Var1 tiempo.Var6 ...
    tiempo.Var7 tiempo.Var8/100 tiempo.Var1*0 tiempo.Var1*0]));
fecha{2} =  [t0;t1];

%---------------------- BSM
tiempo = readtable("fechas_BSM_sig_fil.txt");
a = tiempo.Var5;
tiempo.Var1 = str2num(cell2mat(cellfun(@(x) x(end-3:min(end, numel(x))), a, 'UniformOutput', false)));

t0 = datetime(double([tiempo.Var1 tiempo.Var6 ...
    tiempo.Var7 tiempo.Var8/100 tiempo.Var1*0 tiempo.Var1*0]));

tiempo = readtable("fechas_BSM_signature_2");
a = tiempo.Var5;
tiempo.Var1 = str2num(cell2mat(cellfun(@(x) x(end-3:min(end, numel(x))), a, 'UniformOutput', false)));

t1 = datetime(double([tiempo.Var1 tiempo.Var6 ...
    tiempo.Var7 tiempo.Var8/100 tiempo.Var1*0 tiempo.Var1*0]));
fecha{3} =  [t0;t1];

%---------------------- ITS
% time:units = "since 1970-01-01 00:00:00" ;
time0 = dlmread('dias_ITS_1.txt');time0 = time0(:);time0(time0==0) = [];
time1 = dlmread('dias_ITS_2.txt');time1 = time1(:);time1(time1==0) = [];
time2 = dlmread('dias_ITS_3.txt');time2 = time2(:);time2(time2==0) = [];
fecha{4}  = datetime(datevec([time2;time0;time1]));

%%
close all; clc
% ------------- Figura
F1 = figure(1);

% Dimensiones en orientación landscape
anchoCarta = 11; % Ancho de una hoja carta en pulgadas (landscape)
altoCarta = 8.5; % Alto de una hoja carta en pulgadas (landscape)

% Configurar las propiedades de la figura
set(gcf, 'Units', 'inches', ...         % Unidades en pulgadas
         'Position', [1, 1, anchoCarta, altoCarta], ... % Posición y tamaño de la figura
         'PaperUnits', 'inches', ...    % Unidades del papel
         'PaperSize', [anchoCarta, altoCarta], ... % Tamaño del papel
         'PaperOrientation', 'landscape', ... % Orientación horizontal
         'PaperPositionMode', 'manual', ... % Control manual del papel
         'PaperPosition', [0, 0, anchoCarta, altoCarta]); % Ocupa toda la hoja


s1 = subplot(4,1,[1 2]);
P1 = plot(fecha{1},ones(numel(fecha{1}),1),'.k');
hold on
P2 = plot(fecha{2},ones(numel(fecha{2}),1)*2,'.r');
P3 = plot(fecha{3},ones(numel(fecha{3}),1)*3,'.b');
P4 = plot(fecha{4},ones(numel(fecha{4}),1)*4,'.g');

ms = 6;
P1.MarkerSize = ms;
P2.MarkerSize = ms;
P3.MarkerSize = ms;
P4.MarkerSize = ms;

s1.XGrid = 'on';
s1.YGrid = 'on';

s1.YLim = [0 5];

% datetick('x','mm/yy','keeplimits')

s1.XTick = [s1.XLim(1):calmonths(3):s1.XLim(end)];
s1.XTickLabelRotation = 90;
s1.YTick = [1:4];
s1.YTickLabel = {'$\rm EPB_{WERA}^{HF\ Radar}$', '$\rm PM_{Signature}^{ADCP}$'...
    '$\rm BSM_{Signature}^{ADCP}$' '$\rm ITS_{Aquadopp}^{ADCP}$'};
s1.TickLabelInterpreter = 'Latex'


% Crear proxies para la leyenda con marcador más grande
P01 = plot(nan, nan, '.k');
P02 = plot(nan, nan, '.r');
P03 = plot(nan, nan, '.b');
P04 = plot(nan, nan, '.g');

ms = 8;
P01.MarkerSize = ms;
P02.MarkerSize = ms;
P03.MarkerSize = ms;
P04.MarkerSize = ms;

LE = legend([P01 P02 P03 P04],{'$\rm EPB_{WERA}^{HF\ Radar}$', '$\rm PM_{Signature}^{ADCP}$'...
    '$\rm BSM_{Signature}^{ADCP}$' '$\rm ITS_{Aquadopp}^{ADCP}$'},...
    'Interpreter','latex','Position',[0.10 0.03 0.36 0.29],'FontSize',12);

LE.Position = [0.443341404358355 0.318428184281843 0.215254237288135 0.126030691496491];

print(F1,'-depsc2','serie_tempral')
print(F1,'-dpng','serie_tempral')
print(F1,'-dpdf','serie_tempral')