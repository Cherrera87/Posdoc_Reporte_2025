function Corrige_declinacion_ADCP_signature()
% Corrige_declinacion_ADCP_signature
%
% DESCRIPCIÓN:
%   Aplica la corrección por declinación magnética a las velocidades
%   horizontales (ENU) medidas por ADCP Signature, convirtiéndolas desde
%   un sistema de referencia magnético a un sistema geográfico verdadero.
%
%   La corrección se realiza muestra a muestra utilizando el modelo IGRF,
%   y se almacena directamente en cada archivo .mat de burst, agregando
%   nuevos campos a la estructura Data.
%
% CAMPOS DE ENTRADA REQUERIDOS (en Data):
%   - Burst_Time            : tiempo (datenum)
%   - Burst_Velocity_ENU    : [Nt x 3 x Nz] velocidades ENU magnéticas
%
% CAMPOS GENERADOS (en Data):
%   - Burst_Velocity_ENU_true : velocidades ENU referidas al norte verdadero
%   - Declinacion_deg         : declinación magnética (° Este positivo)
%   - Heading_true (opcional) : rumbo corregido, si existe Burst_Heading
%
% CONFIGURACIÓN:
%   Se define al inicio del script (lat, lon, carpetas, paralelización).
%
% REFERENCIAS:
%   - Emery, W. J., & Thomson, R. E. (2001). Data Analysis Methods in
%     Physical Oceanography. Elsevier.
%   - International Association of Geomagnetism and Aeronomy (IAGA),
%     International Geomagnetic Reference Field (IGRF).
%
% AUTOR:
%   Carlos F. Herrera Vázquez
%
% FECHA:
%   2025-12
%
% -------------------------------------------------------------------------

%% ===================== CONFIGURACIÓN DEL USUARIO ===================== %%
cfg = struct();

% Coordenadas del instrumento (escalares o vectores)
cfg.lat   = 31.800;     % grados Norte
cfg.lon   = -116.700;   % grados Este
cfg.alt_m = 0;          % altitud sobre el nivel del mar [m]

% Carpeta que contiene subcarpetas *_Burst
cfg.parentBurstDir = '/ruta/a/carpeta_con_bursts/';

% Procesamiento paralelo
cfg.useParallel = true;

% Verbosidad
cfg.verbose = true;

%% ===================== VALIDACIONES BÁSICAS ========================== %%
assert(isfolder(cfg.parentBurstDir), ...
    'La carpeta especificada no existe: %s', cfg.parentBurstDir);

folders = dir(fullfile(cfg.parentBurstDir,'*_Burst'));
folders = folders([folders.isdir]);

if isempty(folders)
    error('No se encontraron carpetas *_Burst en %s', cfg.parentBurstDir);
end

%% ===================== PARALELIZACIÓN ================================ %%
if cfg.useParallel
    try
        if isempty(gcp('nocreate'))
            nW = max(feature('numcores')-1,1);
            parpool('Processes', nW);
        end
    catch ME
        warning('No se pudo iniciar parpool. Continuando en serial. (%s)', ME.message);
        cfg.useParallel = false;
    end
end

%% ===================== LOOP PRINCIPAL ================================ %%
loopFolders = 1:numel(folders);

if cfg.useParallel
    parfor ifol = loopFolders
        procesa_carpeta_burst(folders(ifol), cfg);
    end
else
    for ifol = loopFolders
        procesa_carpeta_burst(folders(ifol), cfg);
    end
end

if cfg.verbose
    fprintf('\nCorrección de declinación completada.\n');
end

end

%% ===================================================================== %%
%% ========================= FUNCIONES ================================= %%
%% ===================================================================== %%

function procesa_carpeta_burst(folderInfo, cfg)

path_files = fullfile(folderInfo.folder, folderInfo.name);
files = dir(fullfile(path_files,'*.mat'));

for ii = 1:numel(files)

    fname = fullfile(files(ii).folder, files(ii).name);

    if cfg.verbose
        fprintf('Procesando: %s\n', fname);
    end

    try
        S = load(fname,'Data');
        if ~isfield(S,'Data')
            warning('No existe "Data" en %s. Se omite.', fname);
            continue
        end

        Data = S.Data;

        % Campos mínimos
        req = {'Burst_Time','Burst_Velocity_ENU'};
        if ~all(isfield(Data, req))
            warning('Faltan campos requeridos (%s) en %s.', ...
                strjoin(req,', '), fname);
            continue
        end

        integrar_y_guardar(fname, Data, cfg);

        if cfg.verbose
            fprintf('  OK\n');
        end

    catch ME
        warning('Error en %s: %s', fname, ME.message);
    end
end
end

% --------------------------------------------------------------------- %

function integrar_y_guardar(ruta_mat, Data, cfg)

tiempo = datetime(Data.Burst_Time,'ConvertFrom','datenum');

[Data_corr, ~] = corregir_declinacion_ENU_igrf( ...
    Data, cfg.lat, cfg.lon, cfg.alt_m, tiempo);

% Copiar resultados
Data.Burst_Velocity_ENU_true = Data_corr.Burst_Velocity_ENU_true;
Data.Declinacion_deg         = Data_corr.Declinacion_deg;

if isfield(Data_corr,'Heading_true')
    Data.Heading_true = Data_corr.Heading_true;
end

save(ruta_mat, 'Data', '-append');
end

% --------------------------------------------------------------------- %

function [Data_out, D_vec] = corregir_declinacion_ENU_igrf(Data, lat, lon, alt_m, tiempo)

V = Data.Burst_Velocity_ENU;        % [Nt x 3 x Nz]
Nt = size(V,1);
Nc = size(V,3);

tiempo = tiempo(:);
assert(numel(tiempo)==Nt,'Burst_Time debe coincidir con Nt.');

% Expandir coordenadas
lat   = expand_vector(lat,   Nt);
lon   = expand_vector(lon,   Nt);
alt_m = expand_vector(alt_m, Nt);

% Declinación magnética (IGRF)
D_vec = declinacion_igrf_vec(lat, lon, alt_m, tiempo);   % [Nt x 1]

% Componentes horizontales
E = double(squeeze(V(:,1,:)));   % Este magnético
N = double(squeeze(V(:,2,:)));   % Norte magnético

ct = cosd(D_vec);
st = sind(D_vec);

% Rotación ENU magnético → ENU verdadero
E_true = bsxfun(@times,E,ct) + bsxfun(@times,N,st);
N_true = -bsxfun(@times,E,st) + bsxfun(@times,N,ct);

V_true = V;
V_true(:,1,:) = reshape(single(E_true), Nt, 1, Nc);
V_true(:,2,:) = reshape(single(N_true), Nt, 1, Nc);

Data_out = Data;
Data_out.Burst_Velocity_ENU_true = V_true;
Data_out.Declinacion_deg = single(D_vec);

% Rumbo (si existe)
if isfield(Data,'Burst_Heading') && numel(Data.Burst_Heading)==Nt
    Data_out.Heading_true = mod(double(Data.Burst_Heading(:)) + D_vec, 360);
end
end

% --------------------------------------------------------------------- %

function D_deg = declinacion_igrf_vec(lat, lon, alt_m, fecha)
% Calcula declinación magnética usando IGRF
% D positiva hacia el Este

y  = year(fecha);
d0 = datetime(y,1,1);
d1 = datetime(y+1,1,1);
decYear = y + days(fecha - d0) ./ days(d1 - d0);

[~, ~, D_deg, ~, ~] = igrfmagm(alt_m/1000, lat, lon, decYear);
end

% --------------------------------------------------------------------- %

function v = expand_vector(v, N)
if isscalar(v)
    v = repmat(v, N, 1);
else
    v = v(:);
end
end

