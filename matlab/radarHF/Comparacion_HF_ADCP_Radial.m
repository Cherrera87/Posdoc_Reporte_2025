function U = Comparacion_HF_ADCP_Radial(path_files, path_save)
% =========================================================
% COMPARACION_HF_ADCP_RADIAL
%
% Compara mediciones de corrientes superficiales obtenidas
% con radar HF (WERA) y perfiles de velocidad ADCP,
% construyendo una velocidad "tipo radar" a partir del ADCP
% mediante un promedio vertical ponderado con la función
% de Bragg.
%
% INPUTS:
%   path_files : ruta donde se encuentran los archivos .mat
%                sincronizados HF–ADCP
%   path_save  : ruta donde se guardará la salida
%
% OUTPUT:
%   U : estructura con velocidades HF y ADCP (PM, BSM)
%
% Autor: Carlos F. Herrera Vázquez
% Fecha: Diciembre 2025
% =========================================================

close all; clc;

%% ================== CONFIGURACION ========================

files = dir(fullfile(path_files,'*.mat'));

% Eliminar archivos no útiles (ajusta si cambia el criterio)
files(1:9) = [];

% Coordenadas de sitios
xPM.Lat  = 31.857783;  xPM.Lon  = -116.682983;
xBSM.Lat = 31.853867;  xBSM.Lon = -116.804500;

% Azimut del haz radar (grados, desde Este CCW)
angle.pm  = 157.2;
angle.bsm = 124.2;

% Profundidades medias locales (m)
m_pm  = 20.1092;
m_bsm = 15.8613;

% Longitud de onda de Bragg (WERA ~25 MHz)
lambda_Bragg = 6;                % [m]
kBragg = 2*pi/lambda_Bragg;

%% ================== PREALOCACION =========================

nFiles = numel(files);

U.PM.HF   = nan(nFiles,1);
U.BSM.HF  = nan(nFiles,1);
U.PM.ADCP = nan(nFiles,1);
U.BSM.ADCP= nan(nFiles,1);

%% ================== LOOP PRINCIPAL =======================

for ii = 1:nFiles

    fprintf('Procesando %d de %d\n', ii, nFiles);

    fname = fullfile(files(ii).folder,files(ii).name);
    load(fname,'WERA','PM','BSM','pm','bsm')

    %% -------- CONTROL DE CALIDAD ADCP --------------------
    yPM  = qcUr_eta(PM.Data,  PM.Config, ...
                    pm.Z,  pm.P,  pm.U,  pm.Ur, pm.x, m_pm);

    yBSM = qcUr_eta(BSM.Data, BSM.Config, ...
                    bsm.Z, bsm.P, bsm.U, bsm.Ur, bsm.x, m_bsm);

    %% -------- INTERPOLACION DEL RADAR HF -----------------
    Ux = WERA.u .* cosd(WERA.ang);
    Vy = WERA.u .* sind(WERA.ang);

    Fu = scatteredInterpolant(WERA.lon(:),WERA.lat(:),Ux(:),...
                              'linear','nearest');
    Fv = scatteredInterpolant(WERA.lon(:),WERA.lat(:),Vy(:),...
                              'linear','nearest');

    % PM
    u_pm = Fu(xPM.Lon,xPM.Lat);
    v_pm = Fv(xPM.Lon,xPM.Lat);
    U.PM.HF(ii) = sign(atan2d(v_pm,u_pm)) * hypot(u_pm,v_pm);

    % BSM
    u_bsm = Fu(xBSM.Lon,xBSM.Lat);
    v_bsm = Fv(xBSM.Lon,xBSM.Lat);
    U.BSM.HF(ii) = sign(atan2d(v_bsm,u_bsm)) * hypot(u_bsm,v_bsm);

    %% -------- ADCP "TIPO RADAR HF" -----------------------
    U.PM.ADCP(ii)  = radar_sintetico_from_adcp(yPM,  kBragg);
    U.BSM.ADCP(ii) = radar_sintetico_from_adcp(yBSM, kBragg);

end

%% ================== GUARDADO =============================

if ~exist(path_save,'dir')
    mkdir(path_save)
end

save(fullfile(path_save,'HF_ADCP_comparacion.mat'),'U');

end


%----------------------------------------
%-----------Funciones--------------------

function Ur_mean = radar_sintetico_from_adcp(ySite,kBragg)

Ur = ySite.Ur_eta;     % [Nt x Nz]
Z  = ySite.Z;          % [Nt x Nz]

[Nt,~] = size(Ur);

zcol = mean(Z,1,'omitnan')';
[zcol,idx] = sort(zcol,'ascend');
Ur = Ur(:,idx);

h = abs(min(zcol));

eK = cosh(2*kBragg*(zcol + h)) ./ (sinh(kBragg*h).^2);
Wz = eK / trapz(zcol,eK);

Ur_eff = nan(Nt,1);
for it = 1:Nt
    prof = Ur(it,:)';
    ok = ~isnan(prof);
    if sum(ok) > 2
        Ur_eff(it) = trapz(zcol(ok), prof(ok).*Wz(ok));
    end
end

Ur_mean = nanmean(Ur_eff);

end


function y = qcUr_eta(Data,Config,Z,P,U,Ur,x,m_ref)

y.eta = P - m_ref;

Nt = size(Ur,1);
Nz = size(Ur,2);

mask = true(Nt,Nz);

for b = [1 2]
    corr = squeeze(Data.Burst_Correlation_Beam(:,b,:));
    amp  = squeeze(Data.Burst_Amplitude_Beam(:,b,:));
    mask = mask & corr > 60 & amp > 20;
end

Ur(~mask) = NaN;

z_surf = -fliplr(0:Config.Burst_CellSize:m_ref);
y.Z = z_surf;

for it = 1:Nt
    y.Ur_eta(it,:) = interp1(Z - y.eta(it),Ur(it,:),z_surf,'linear',NaN);
end

y.Ur = Ur;

end


