function Convierte_CRAD_a_MAT_WERA(path_crad, path_save)
%==========================================================================
% Convierte_CRAD_a_MAT_WERA
%
% Convierte archivos .crad generados por el radar HF WERA a archivos .mat
% utilizando el toolbox matWERA (Voulgaris et al.).
%
% Cada archivo .crad es leído individualmente y almacenado como un archivo
% .mat independiente, preservando la estructura original de las variables.
%
% ENTRADAS:
%   path_crad : ruta raíz donde se localizan los archivos .crad
%   path_save : directorio donde se guardarán los archivos .mat
%
% SALIDAS:
%   Archivos .mat con campos:
%     Time, lat, lon, x, y, u, uvar, uacc, pwr, ang, Range
%
% NOTAS:
%   - Requiere el toolbox matWERA en el path de MATLAB.
%   - El nombre del archivo se construye a partir de la marca temporal
%     contenida en el nombre del archivo .crad.
%
%==========================================================================

close all; clc;

% -------------------------------------------------------------------------
% Verificación de rutas
% -------------------------------------------------------------------------
if ~isfolder(path_crad)
    error('La ruta de archivos .crad no existe.');
end

if ~isfolder(path_save)
    mkdir(path_save)
end

% -------------------------------------------------------------------------
% Búsqueda recursiva de archivos .crad
% -------------------------------------------------------------------------
files = dir(fullfile(path_crad, '**', '*.crad'));

if isempty(files)
    warning('No se encontraron archivos .crad en la ruta especificada.');
    return
end

fprintf('Archivos .crad encontrados: %d\n', length(files));

% -------------------------------------------------------------------------
% Lectura y conversión
% -------------------------------------------------------------------------
for ii = 1:length(files)

    fprintf('Procesando archivo %d de %d\n', ii, length(files));

    try
        % Nombre completo del archivo
        fname = fullfile(files(ii).folder, files(ii).name);

        % Lectura del archivo CRAD
        [dat.Time, dat.lat, dat.lon, dat.x, dat.y, ...
         dat.u, dat.uvar, dat.uacc, dat.pwr, ...
         dat.ang, dat.Range] = read_WERA_crad(fname);

        % -------------------------------------------------------------
        % Construcción de la fecha a partir del nombre del archivo
        % Formato típico: yyyyDDDHHmm
        % -------------------------------------------------------------
        fecha = datetime(files(ii).name(1:11), ...
                         'InputFormat', 'yyyyDDDHHmm');

        % Nombre de salida
        fname_out = fullfile(path_save, ...
            ['crad_' datestr(fecha,'yyyymmddHHMM') '.mat']);

        % Guardado
        save(fname_out, '-struct', 'dat', '-v7.3');

    catch ME
        warning('Error al procesar %s\n%s', files(ii).name, ME.message);
    end

    clear dat

end

fprintf('Conversión finalizada.\n');

end

