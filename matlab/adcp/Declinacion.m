close all; clear; clc

lat = 31.800; lon = -116.700; alt_m = 0;
folders = dir('/media/root01/Pancho_SSD/ADCP-DATA/*_Burst');

parfor ifol = 1:length(folders)
    path_files = fullfile(folders(ifol).folder,folders(ifol).name);
    files = dir(fullfile(path_files,'*.mat'));

    for ii = 1:length(files)
        fname = fullfile(files(ii).folder, files(ii).name);
        fprintf('Procesando: %s\n', fname);
        try
            S = load(fname,'Data');
            if ~isfield(S,'Data')
                warning('No hay variable "Data" en %s. Saltando.', fname);
                continue
            end
            Data = S.Data;

            % Checks mínimos
            req = {'Burst_Time','Burst_Velocity_ENU'};
            if ~all(isfield(Data, req))
                warning('Faltan campos requeridos (%s) en %s. Saltando.', strjoin(req,', '), fname);
                continue
            end

            integrar_corr_en_Data_y_guardar(fname, Data, lat, lon, alt_m);
            fprintf('OK: %s\n', files(ii).name);
        catch ME
            warning('Error en %s: %s', files(ii).name, ME.message);
            % continúa con el siguiente archivo
        end
    end
end

%--------------------------------------------------------------------------
% ------------------------------- Funciones -------------------------------

function integrar_corr_en_Data_y_guardar(ruta_mat, Data, lat, lon, alt_m)
tiempo = datetime(Data.Burst_Time,'ConvertFrom','datenum');

[Data_corr, ~] = corregir_declinacion_ENU_igrf(Data, lat, lon, alt_m, tiempo);

% Copiamos campos útiles de Data_corr -> Data
Data.Burst_Velocity_ENU_true = Data_corr.Burst_Velocity_ENU_true;
Data.Declinacion_deg         = Data_corr.Declinacion_deg;
if isfield(Data_corr, 'Heading_true')
    Data.Heading_true = Data_corr.Heading_true;
end

% Guardar: crea con -v7.3 si no existe; si existe, usa -append
if isfile(ruta_mat)
    save(ruta_mat, 'Data', '-append');     % reescribe Data
else
    save(ruta_mat, 'Data', '-v7.3');       % crea nuevo grande
end
end

function [Data_out, D_vec] = corregir_declinacion_ENU_igrf(Data, lat, lon, alt_m, tiempo)
V = Data.Burst_Velocity_ENU;
Nt = size(V,1); Nc = size(V,3);
tiempo = tiempo(:);
assert(numel(tiempo)==Nt,'tiempo debe tener Nt elementos.');

% Expandir lat/lon/alt si son escalares
if isscalar(lat),   lat   = repmat(lat,   Nt,1); else, lat = lat(:);   end
if isscalar(lon),   lon   = repmat(lon,   Nt,1); else, lon = lon(:);   end
if isscalar(alt_m), alt_m = repmat(alt_m, Nt,1); else, alt_m = alt_m(:); end

% Declinación por muestra (IGRF), grados Este positivo
D_vec = declinacion_igrf_vec(lat, lon, alt_m, tiempo);   % [Nt x 1]

% Rotación horizontal magnético -> verdadero
E = double(squeeze(V(:,1,:)));   % [Nt x Nc] E
N = double(squeeze(V(:,2,:)));   % [Nt x Nc] N
ct = cosd(D_vec);  st = sind(D_vec);

% Compatibilidad amplia: bsxfun
E_true = bsxfun(@times,E,ct) + bsxfun(@times,N,st);
N_true = -bsxfun(@times,E,st) + bsxfun(@times,N,ct);

V_true = V;                          % copia
V_true(:,1,:) = reshape(single(E_true), Nt, 1, Nc);
V_true(:,2,:) = reshape(single(N_true), Nt, 1, Nc);

Data_out = Data;
Data_out.Burst_Velocity_ENU_true = V_true;
Data_out.Declinacion_deg = single(D_vec);

if isfield(Data, 'Burst_Heading') && numel(Data.Burst_Heading)==Nt
    Data_out.Heading_true = mod(double(Data.Burst_Heading(:)) + D_vec, 360);
end
end

function D_deg = declinacion_igrf_vec(lat, lon, alt_m, fecha)
% fecha puede ser vector datetime
y  = year(fecha);
d0 = datetime(y,1,1);
d1 = datetime(y+1,1,1);
decYear = y + days(fecha - d0) ./ days(d1 - d0);
% igrfmagm: [XYZ,H,D,I,F] -> D es 5ª salida
[~, ~, D_deg, ~, ~] = igrfmagm(alt_m/1000, lat, lon, decYear);
end
