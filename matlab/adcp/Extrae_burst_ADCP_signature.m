function Extrae_burst_ADCP_signature(path_input, path_output, idCaso)
% EXTRAE_BURST_ADCP_SIGNATURE
%
% Segmenta archivos .mat (ADCP Signature nivel 0) en bursts individuales,
% utilizando Data.Burst_EnsembleCount para detectar reinicios de conteo.
%
% ENTRADAS:
%   path_input  : carpeta con archivos .mat nivel 0
%   path_output : carpeta donde se guardarán los bursts
%   idCaso      : prefijo de nombre para los archivos de salida
%
% SALIDAS:
%   Archivos .mat en path_output con nombre:
%     <idCaso><yyyymmddHHMM>_b###.mat
%
% AUTOR:
%   Carlos F. Herrera Vázquez
%
% FECHA:
%   2025-12
% ====================================================================== %

%% ======================= VALIDACIONES ======================= %%
if nargin < 2
    error('Se requieren al menos path_input y path_output.');
end
if nargin < 3 || isempty(idCaso)
    idCaso = 'Sig_';
end

assert(isfolder(path_input), 'inputDir no existe: %s', path_input);

if ~isfolder(path_output)
    mkdir(path_output);
end

%% ======================= CONFIGURACIÓN ======================= %%
cfg = struct();

cfg.inputDir      = path_input;
cfg.outputDir     = path_output;
cfg.idCaso        = idCaso;

cfg.filePattern   = '*.mat';
cfg.dropLastFile  = false;
cfg.verbose       = true;

% Paralelización (normalmente innecesaria aquí)
cfg.useParallel   = false;
cfg.maxWorkers    = [];

%% ======================= PARALLEL POOL (OPCIONAL) ======================= %%
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

%% ======================= LISTA DE ARCHIVOS ======================= %%
files = dir(fullfile(cfg.inputDir, cfg.filePattern));
if isempty(files)
    error('No se encontraron archivos .mat en %s', cfg.inputDir);
end

% Ordenar por nombre
[~, idxSort] = sort({files.name});
files = files(idxSort);

if cfg.dropLastFile && numel(files) >= 1
    files(end) = [];
end

if cfg.verbose
    fprintf('Encontrados %d archivos en %s\n', numel(files), cfg.inputDir);
end

%% ======================= PROCESAMIENTO ======================= %%
for ii = 1:numel(files)

    fpath = fullfile(files(ii).folder, files(ii).name);

    if cfg.verbose
        fprintf('\n[%d/%d] Procesando: %s\n', ii, numel(files), files(ii).name);
    end

    data = load(fpath);

    %% Validar estructura mínima
    requiredTop = {'Config','Descriptions','Units','Data'};
    ok = true;
    for k = 1:numel(requiredTop)
        if ~isfield(data, requiredTop{k})
            warning('Archivo %s NO contiene "%s". Se omite.', ...
                    files(ii).name, requiredTop{k});
            ok = false;
        end
    end
    if ~ok, continue; end

    if ~isfield(data.Data,'Burst_EnsembleCount') || ...
       ~isfield(data.Data,'Burst_Time')
        warning('Archivo %s no tiene Burst_EnsembleCount o Burst_Time. Se omite.', ...
                files(ii).name);
        continue;
    end

    burstCount = data.Data.Burst_EnsembleCount(:);
    if numel(burstCount) < 2
        warning('Archivo %s: Burst_EnsembleCount muy corto. Se omite.', ...
                files(ii).name);
        continue;
    end

    %% Campos a copiar
    campos = fieldnames(data.Data);

    % Excluir campos que empiezan con 'IBurst'
    ind = ~cellfun('isempty', regexp(campos, '^IBurst', 'once'));
    campos = campos(~ind);

    %% Detectar límites de burst
    jumps = find(diff(burstCount) < 1);
    p_ini = [1; jumps + 1];
    p_end = [jumps; numel(burstCount)];

    nBursts = numel(p_ini);

    if cfg.verbose
        fprintf('  Bursts detectados: %d\n', nBursts);
    end

    %% Loop por burst
    for ib = 1:nBursts

        idxBurst = p_ini(ib):p_end(ib);

        S = struct();
        S.Config       = data.Config;
        S.Descriptions = data.Descriptions;
        S.Units        = data.Units;
        S.Data         = struct();

        %% Copiar campos de Data
        for ik = 1:numel(campos)

            fname = campos{ik};

            try
                A = data.Data.(fname);

                if isempty(A)
                    S.Data.(fname) = A;
                    continue;
                end

                nd = ndims(A);
                subs = repmat({':'}, 1, nd);
                subs{1} = idxBurst;

                S.Data.(fname) = A(subs{:});

            catch ME
                warning('  (burst %d) No pude recortar campo "%s": %s. Campo omitido.', ...
                        ib, fname, ME.message);
            end
        end

        %% Timestamp
        try
            t0 = S.Data.Burst_Time(1);
            fecha = datestr(t0, 'yyyymmddHHMM');
        catch
            fecha = regexprep(files(ii).name, '\W', '');
        end

        %% Guardado
        outName = sprintf('Burst_%s%s.mat', cfg.idCaso, fecha);
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

