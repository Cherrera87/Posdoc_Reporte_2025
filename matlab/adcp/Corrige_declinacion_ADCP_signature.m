function Corrige_declinacion_ADCP_signature(parentBurstDir, lat, lon, alt_m, useParallel, verbose)
% Corrige_declinacion_ADCP_signature
%
% Aplica corrección por declinación magnética a velocidades ENU medidas por
% ADCP Signature (referidas a norte magnético) para obtener ENU verdadero.
% Se calcula muestra a muestra con IGRF y se guarda en cada archivo burst:
%   - Data.Burst_Velocity_ENU_true
%   - Data.Declinacion_deg
%   - Data.Heading_true (si existe Burst_Heading)
%
% ENTRADAS:
%   parentBurstDir : carpeta que contiene subcarpetas *_Burst
%   lat, lon, alt_m: coordenadas (escalares o vectores por muestra)
%   useParallel    : true/false
%   verbose        : true/false
%
% Autor: Carlos F. Herrera Vázquez
% Fecha: 2025-12
% -------------------------------------------------------------------------

%% ===================== DEFAULTS / VALIDACIÓN ===================== %%
if nargin < 1 || isempty(parentBurstDir)
    error('Se requiere parentBurstDir');
end
if nargin < 2 || isempty(lat),  lat = 31.800;   end
if nargin < 3 || isempty(lon),  lon = -116.700; end
if nargin < 4 || isempty(alt_m), alt_m = 0;     end
if nargin < 5 || isempty(useParallel), useParallel = true; end
if nargin < 6 || isempty(verbose),     verbose = true;     end

cfg = struct();
cfg.lat = lat;
cfg.lon = lon;
cfg.alt_m = alt_m;
cfg.parentBurstDir = parentBurstDir;
cfg.useParallel = useParallel;
cfg.verbose = verbose;

assert(isfolder(cfg.parentBurstDir), ...
    'La carpeta especificada no existe: %s', cfg.parentBurstDir);

files = dir(fullfile(cfg.parentBurstDir,'Burst*.mat'));

if isempty(files)
    error('No se encontraron archivos de Burst en %s', cfg.parentBurstDir);
end

%% ===================== PARALELIZACIÓN ===================== %%
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

%% ===================== LOOP PRINCIPAL ===================== %%
loopFiles = 1:numel(files);

if cfg.useParallel
    parfor ifol = loopFiles
        fid= files(ifol);
        procesa_carpeta_burst(fid, cfg);
    end
else
    for ifol = loopFiles
        fid = files(ifol);
        procesa_carpeta_burst(fid, cfg);
    end
end

if cfg.verbose
    fprintf('\nCorrección de declinación completada.\n');
end

end

%% ===================================================================== %%
%% ========================= FUNCIONES ================================= %%
%% ===================================================================== %%

function procesa_carpeta_burst(fid, cfg)
    fname = fullfile(fid.folder,fid.name);
    if cfg.verbose
        fprintf('Procesando: %s\n', fid.name);
    end


    try
        S = load(fname,'Data');
        if ~isfield(S,'Data')
            warning('No existe "Data" en %s. Se omite.', fid.name);
        end

        Data = S.Data;

        % Campos mínimos
        req = {'Burst_Time','Burst_Velocity_ENU'};
        if ~all(isfield(Data, req))
            warning('Faltan campos requeridos (%s) en %s.', ...
                strjoin(req,', '), fname);
        end

        integrar_y_guardar(fname, Data, cfg);

        if cfg.verbose
            fprintf('  OK\n');
        end

    catch ME
        warning('Error en %s: %s', fname, ME.message);
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

% Componentes horizontales (magnéticas)
E = double(squeeze(V(:,1,:)));   % [Nt x Nc]
N = double(squeeze(V(:,2,:)));   % [Nt x Nc]

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

% igrfmagm espera alt en km
[~, ~, D_deg, ~, ~] = igrfmagm(alt_m/1000, lat, lon, decYear);
end

% --------------------------------------------------------------------- %

function v = expand_vector(v, N)
if isscalar(v)
    v = repmat(v, N, 1);
else
    v = v(:);
    assert(numel(v)==N, 'Vector debe ser escalar o de longitud Nt.');
end
end
