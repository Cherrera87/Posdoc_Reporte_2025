function Extrae_burst_ADCP_signature()
% Extrae_burst_ADCP_signature
%
% DESCRIPCIÓN:
%   Segmenta archivos .mat (ADCP Signature nivel 0) en bursts individuales,
%   utilizando la variable Data.Burst_EnsembleCount para detectar reinicios
%   de conteo y así delimitar cada burst.
%
%   El script carga cada archivo, conserva metadatos (Config, Descriptions,
%   Units) y guarda una salida .mat por burst con los campos relevantes en S.
%
% ENTRADAS:
%   No recibe argumentos. La configuración se define en la sección
%   "CONFIGURACIÓN DEL USUARIO" al inicio del script:
%     - inputDir: carpeta con archivos .mat nivel 0
%     - outputDir: carpeta donde se guardarán bursts
%     - idCaso: prefijo de nombre para los archivos
%
% SALIDAS:
%   Archivos .mat en outputDir con nombre:
%      <idCaso><yyyymmddHHMM>_b###.mat
%   donde yyyymmddHHMM viene de Burst_Time(1) y ### es el número de burst.
%
% NOTAS:
%   - Requiere que cada archivo .mat contenga al menos:
%       Config, Descriptions, Units, Data.Burst_EnsembleCount, Data.Burst_Time
%   - Si algunos campos de Data tienen dimensiones diferentes (p.ej. no
%     indexables por la primera dimensión), se omiten con advertencia.
%
% REFERENCIAS:
%   - Emery, W. J., & Thomson, R. E. (2001). Data Analysis Methods in Physical
%     Oceanography. Elsevier.
%   - Manual/Guía del fabricante del ADCP Signature (Teledyne RDI / Nortek,
%     según aplique a tu conversión).
%
% AUTOR:
%   Carlos F. Herrera Vázquez
%
% FECHA:
%   2025-12
%
% -------------------------------------------------------------------------

%% CONFIGURACIÓN DEL USUARIO (edita aquí)
cfg = struct();

% Carpeta con archivos .mat (nivel 0)
cfg.inputDir  = '/ruta/a/tus/mats/level0/';

% Carpeta donde se guardarán los bursts
cfg.outputDir = '/ruta/a/salida/bursts/';

% Prefijo de nombre
cfg.idCaso    = 'Sig_BSM_';

% Patrón de búsqueda (por defecto *.mat)
cfg.filePattern = '*.mat';

% Si deseas excluir el último archivo (por ejemplo por archivo incompleto)
cfg.dropLastFile = false;

% Control de verbosidad
cfg.verbose = true;

% (Opcional) Procesamiento paralelo: normalmente NO hace falta aquí.
cfg.useParallel = false;     % true/false
cfg.maxWorkers  = [];        % [] usa feature('numcores')-1; o pon un número.

%% Validaciones básicas
assert(isfolder(cfg.inputDir),  'inputDir no existe: %s',  cfg.inputDir);
if ~isfolder(cfg.outputDir)
    mkdir(cfg.outputDir);
end

%% (Opcional) Parallel pool
if cfg.useParallel
    try
        if isempty(gcp('nocreate'))
            if isempty(cfg.maxWorkers)
                nW = max(feature('numcores')-1, 1);
            else
                nW = max(cfg.maxWorkers, 1);
            end
            parpool('Processes', nW);
        end
    catch ME
        warning('No se pudo iniciar parpool. Continuando en serial. Detalle: %s', ME.message);
        cfg.useParallel = false;
    end
end

%% Lista de archivos
files = dir(fullfile(cfg.inputDir, cfg.filePattern));
if isempty(files)
    error('No se encontraron archivos con patrón %s en %s', cfg.filePattern, cfg.inputDir);
end

% Ordenar por nombre (o por fecha si prefieres)
[~, idxSort] = sort({files.name});
files = files(idxSort);

if cfg.dropLastFile && numel(files) >= 1
    files(end) = [];
end

if cfg.verbose
    fprintf('Encontrados %d archivos en %s\n', numel(files), cfg.inputDir);
end

%% Procesamiento
for ii = 1:numel(files)
    fpath = fullfile(files(ii).folder, files(ii).name);

    if cfg.verbose
        fprintf('\n[%d/%d] Procesando: %s\n', ii, numel(files), files(ii).name);
    end

    data = load(fpath);

    % Validar estructura mínima esperada
    requiredTop = {'Config','Descriptions','Units','Data'};
    for k = 1:numel(requiredTop)
        if ~isfield(data, requiredTop{k})
            warning('Archivo %s NO contiene "%s". Se omite.', files(ii).name, requiredTop{k});
            continue;
        end
    end

    if ~isfield(data.Data,'Burst_EnsembleCount') || ~isfield(data.Data,'Burst_Time')
        warning('Archivo %s no tiene Data.Burst_EnsembleCount o Data.Burst_Time. Se omite.', files(ii).name);
        continue;
    end

    burstCount = data.Data.Burst_EnsembleCount(:);
    if numel(burstCount) < 2
        warning('Archivo %s: Burst_EnsembleCount muy corto. Se omite.', files(ii).name);
        continue;
    end

    % Determinar campos a copiar de Data:
    campos = fieldnames(data.Data);

    % Excluir campos que empiezan con 'IBurst' (como en tu versión)
    ind = ~cellfun('isempty', regexp(campos, '^IBurst', 'once'));
    campos = campos(~ind);

    % Detectar límites de burst: reinicio cuando diff < 1
    jumps = find(diff(burstCount) < 1);
    p_ini = [1; jumps + 1];
    p_end = [jumps; numel(burstCount)];

    nBursts = numel(p_ini);

    if cfg.verbose
        fprintf('  Bursts detectados: %d\n', nBursts);
    end

    % Loop por burst
    for ib = 1:nBursts
        idxBurst = p_ini(ib):p_end(ib);

        S = struct();
        S.Config       = data.Config;
        S.Descriptions = data.Descriptions;
        S.Units        = data.Units;
        S.Data         = struct();

        % Copiar campos de Data indexando la primera dimensión
        for ik = 1:numel(campos)
            fname = campos{ik};

            try
                A = data.Data.(fname);

                % Campos vacíos se copian tal cual
                if isempty(A)
                    S.Data.(fname) = A;
                    continue;
                end

                % Si es vector/array indexable por 1a dimensión, recortamos
                nd = ndims(A);
                subs = repmat({':'}, 1, nd);
                subs{1} = idxBurst;
                S.Data.(fname) = A(subs{:});

            catch ME
                % Si no se puede indexar (dimensiones incompatibles), se omite
                warning('  (burst %d) No pude recortar campo "%s": %s. Se omite ese campo.', ...
                    ib, fname, ME.message);
            end
        end

        % Timestamp desde el primer Burst_Time del burst
        try
            t0 = S.Data.Burst_Time(1);
            fecha = datestr(t0, 'yyyymmddHHMM');
        catch
            % fallback: usar nombre del archivo + burst index
            fecha = regexprep(files(ii).name, '\W', '');
        end

        outName = sprintf('%s%s_b%03d.mat', cfg.idCaso, fecha, ib);
        outPath = fullfile(cfg.outputDir, outName);

        save(outPath, '-struct', 'S');

        if cfg.verbose
            fprintf('    Guardado: %s\n', outName);
        end
    end
end

if cfg.verbose
    fprintf('\nListo. Salida en: %s\n', cfg.outputDir);
end

end

