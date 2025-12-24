function Procesa_WERA_ADCP_Sincronizados(path_wera, path_pm, path_bsm, path_save)
%==========================================================================
% Procesa_WERA_ADCP_Sincronizados
%
% Identifica mediciones simultáneas entre radar HF (WERA) y ADCPs,
% proyecta las corrientes ADCP sobre la dirección radial del radar,
% estima oleaje direccional (PUV) y calcula deriva de Stokes direccional.
%
% ENTRADAS:
%   path_wera : carpeta con archivos WERA .mat (convertidos desde .crad)
%   path_pm   : carpeta con bursts ADCP Punta Morro
%   path_bsm  : carpeta con bursts ADCP Bajo San Miguel
%   path_save : carpeta de salida
%
% SALIDA:
%   Archivos HFADCP_YYYYMMDDHHMM.mat
%
%==========================================================================

close all; clc;

if ~isfolder(path_save)
    mkdir(path_save)
end

%------------------------------------------------------------------
% Identificación de mediciones simultáneas
%------------------------------------------------------------------
path_files = {path_wera, path_pm, path_bsm};

for ii = 1:3
    files{ii} = dir(fullfile(path_files{ii}, '*.mat'));
    names = {files{ii}.name};
    fechas_str = cellfun(@(x) x(end-15:end-4), names, 'UniformOutput', false);
    fechas{ii} = datetime(fechas_str, 'InputFormat', 'yyyyMMddHHmm');
end

fechas_comunes = fechas{1};
for k = 2:length(fechas)
    fechas_comunes = intersect(fechas_comunes, fechas{k});
end

for k = 1:length(fechas)
    [~, idx] = ismember(fechas_comunes, fechas{k});
    files{k} = files{k}(idx);
end

%------------------------------------------------------------------
% Parámetros geométricos del radar
%------------------------------------------------------------------
angle.pm  = 157.2;
angle.bsm = 124.2;

m_pm  = 20.1092;
m_bsm = 15.8613;

%------------------------------------------------------------------
% Procesamiento principal
%------------------------------------------------------------------
for ii = 1:length(files{1})

    fprintf('Procesando %d de %d\n', ii, length(files{1}));

    try
        WERA = load(fullfile(files{1}(ii).folder, files{1}(ii).name));
        PM   = load(fullfile(files{2}(ii).folder, files{2}(ii).name));
        BSM  = load(fullfile(files{3}(ii).folder, files{3}(ii).name));

        pm  = procesar_adcp(PM.Data,  PM.Config,  1, angle.pm  + 180, m_pm);
        bsm = procesar_adcp(BSM.Data, BSM.Config, 1, angle.bsm + 180, m_bsm);

        pm.mpres  = m_pm;
        bsm.mpres = m_bsm;

        fname_out = ['HFADCP_' files{1}(ii).name(end-15:end-4) '.mat'];
        save(fullfile(path_save, fname_out), ...
            'WERA','PM','BSM','pm','bsm','-v7.3');

    catch ME
        warning('Error en %s\n%s', files{1}(ii).name, ME.message);
    end
end

fprintf('Procesamiento finalizado.\n');

end

function resultado = procesar_adcp(Data, Config, ik, angle, msl)

resultado.angle = angle;

resultado.z = double(Config.Burst_BlankingDistance + ...
              (1:double(Config.Burst_NCells)) * Config.Burst_CellSize);

resultado.P = double(Data.Burst_Pressure);
resultado.x = datetime(datevec(Data.Burst_Time));
resultado.Z = resultado.z - msl;

resultado.U = double(squeeze(Data.Burst_Velocity_ENU_true(:,1,:)));
resultado.V = double(squeeze(Data.Burst_Velocity_ENU_true(:,2,:)));

resultado.u = detrend(resultado.U(:,ik));
resultado.v = detrend(resultado.V(:,ik));

resultado.fs = double(Config.Burst_SamplingRate);
resultado.dt = 1 / resultado.fs;

resultado.nF = 100;
resultado.hp = 0.6;
resultado.hv = resultado.z(ik);
parms = [0.03 200 0.03 0];
resdir = 2;

[resultado.Su, resultado.Sp, resultado.Dir, resultado.Spread, resultado.F, ...
 resultado.dF, resultado.DOF, ka] = ...
    wds(resultado.u, resultado.v, resultado.P, ...
        resultado.dt, resultado.nF, resultado.hp, resultado.hv, parms);

icod = find(~isnan(resultado.Sp));
if ~isempty(icod)
    [resultado.S2DC, resultado.Xi, resultado.Yi, resultado.ZW, resultado.theta] = ...
        make2dcd(resultado.Sp, resultado.F, resultado.Dir, resultado.Spread, resdir);
    resultado.S2DC = resultado.S2DC * pi/180;
else
    resultado.S2DC = NaN;
    resultado.theta = linspace(0, 2*pi, 36);
end

resultado.info = waveinfo(resultado.Su, resultado.Sp, resultado.Dir, ...
                          resultado.Spread, resultado.F, resultado.dF, ...
                          ka, resultado.P, resultado.hp);

[resultado.us_ztheta, resultado.uz, resultado.vz, resultado.zvec, ...
 resultado.theta, resultado.valid_count] = ...
    radial_stokes_drift_profile(resultado.S2DC, resultado.F, resultado.dF, ...
                                ka, resultado.info.h, resultado.theta, ...
                                resultado.info.zvec);

theta_rad = deg2rad(90 - angle);
resultado.Ur = resultado.U * cos(theta_rad) + resultado.V * sin(theta_rad);

resultado.time = resultado.x(1);
end

%--------------------------------------------------------------------------
% ----------------------------- FUNCIONES ---------------------------------
function resultado = procesar_adcp(Data, Config, ik, angle,msl)
% PROCESAR_ADCP Procesa un conjunto de datos de un ADCP tipo Signature
%
% INPUTS:
%   Data   : estructura con campos Burst_Velocity_ENU, Burst_Pressure, etc.
%   Config : estructura con metadatos del instrumento
%   ik     : índice de celda vertical a procesar (entero positivo)
%   angle  : ángulo en grados para obtener la componente radial (ej. 157.2)
%
% OUTPUT:
%   resultado : estructura con espectros, deriva de Stokes y componente radial

    resultado.angle = angle;

    resultado.z = double(Config.Burst_BlankingDistance + (1:double(Config.Burst_NCells)) * Config.Burst_CellSize);
    resultado.P = double(Data.Burst_Pressure);
    resultado.x = datetime(datevec(Data.Burst_Time));
    resultado.Z = resultado.z - msl;  % profundidad relativa

    resultado.U = double(squeeze(Data.Burst_Velocity_ENU_true(:,1,:)));
    resultado.V = double(squeeze(Data.Burst_Velocity_ENU_true(:,2,:)));

    resultado.u = detrend(resultado.U(:,ik));
    resultado.v = detrend(resultado.V(:,ik));

    resultado.fs = double(Config.Burst_SamplingRate);
    resultado.dt = 1 / resultado.fs;
    resultado.nF = 100;
    resultado.hp = 0.6;
    resultado.hv = resultado.z(ik);
    parms = [0.03 200 0.03 0];
    resdir = 2;

    [resultado.Su, resultado.Sp, resultado.Dir, resultado.Spread, resultado.F, ...
        resultado.dF, resultado.DOF, ka] = wds(resultado.u, resultado.v, resultado.P, ...
                                               resultado.dt, resultado.nF, resultado.hp, resultado.hv, parms);

    icod = find(~isnan(resultado.Sp));
    if ~isempty(icod)
        [resultado.S2DC, resultado.Xi, resultado.Yi, resultado.ZW, resultado.theta] = ...
            make2dcd(resultado.Sp, resultado.F, resultado.Dir, resultado.Spread, resdir);
        resultado.S2DC = resultado.S2DC * pi/180;
    else
        resultado.S2DC = NaN;
        resultado.Xi = NaN;
        resultado.Yi = NaN;
        resultado.ZW = NaN;
        resultado.theta = linspace(0, 2*pi, 36);  % evitar errores si está vacío
    end

    resultado.info = waveinfo(resultado.Su, resultado.Sp, resultado.Dir, resultado.Spread, ...
                              resultado.F, resultado.dF, ka, resultado.P, resultado.hp);

    [resultado.us_ztheta, resultado.uz, resultado.vz, resultado.zvec, resultado.theta, resultado.valid_count] = ...
        radial_stokes_drift_profile(resultado.S2DC, resultado.F, resultado.dF, ka, ...
                                    resultado.info.h, resultado.theta, resultado.info.zvec);

    theta_rad = deg2rad(90 - angle);
    resultado.Ur = resultado.U * cos(theta_rad) + resultado.V * sin(theta_rad);

    resultado.time = resultado.x(1);
end




%--------------------------------- WDS ------------------------------------
function [Su,Sp,Dir,Spread,F,dF,DOF,ka] = wds(u,v,p,dt,nF,hp,hv,parms);
% [Su,Sp,Dir,Spread,F,dF,DOF] = wds(u,v,p,dt,nF,hp,hv,[parms]);
%
% Su is the surface elevation spectra (m^2/Hz), based on the velocity data
% Sp is the surface elevation spectra (m^2/Hz), based on the pressure data
% Dir is the wave direction (deg)
% Spread is the wave spreading (deg)
%     all are in matrices (size = [nf nt], where
%     nf is the number of output frequency bands
%     nt is the number of input time series
% F is the center frequency of each band
% dF is the bandwidth of each band
% DOF is th number of degrees of freedom in each band
%
% u,v are east and north components of velocity (m/s),
% p is pressure (m)
%     all are in matrices (size = [np nt]), where
%     np = time series length (must be even; power of 2 is fastest),
%     nt = number of time series
% dt is the sample interval in s (typically 0.5 or 1 s)
% nF is the nominal number of output frequencies; the result is nf <= nF
% hp is the height of the pressure sensor above the bottom (m)
%     (this means the water depth is the mean pressure plus hp)
% hv is the height of the velocity cell above the pressure sensor (m)
%
% parms = [lf, maxfac, minspec, Ndir]
%   lf is the low-frequency cutoff. F<lf are not output
%   maxfac is the largest factor scaling pressure to surface elevation
%     spectra and directions at F above this cutoff are returned as NaNs
%   minspec is the minimum spectral level for which direction is computed.
%     directions for spectra below this level are returned as NaNs
%   Ndir is the direction of the "north" component (degrees)
%   default: parms = [0.03 200 0.03 0];
%
% Note: The routine is faster if the length of the time series
%   factors into small integers, and fastest if it is a power of 2
%
% Copyright (C) 2001, Lee Gordon, NortekUSA LLC

if nargin<8,
    parms=[0.03 200 0.1 0];
end;
lf=parms(1);
maxfac=parms(2);
minspec=parms(3);
Ndir=parms(4);

[np,nt]=size(p);            % no. of points & no. of time series
if mod(np,2)==1,          % make even number of points
    'Time series lengths must be even'    % (it's just easier)
    return;
end;

Dt=np*dt;
f=[1:(np/2)]/Dt;            % frequency array up to Nyquist

uf=(fft(u));             % compute spectrum
vf=(fft(v));             %
pf=(fft(p));

up=abs(uf(2:(np/2+1),:).^2);       % compute power spectrum & ignore zero freq
vp=abs(vf(2:(np/2+1),:).^2);       % (this uses first half of spectrum)
pp=abs(pf(2:(np/2+1),:).^2);
up=up*2/np^2/f(1);          % scale power spectrum
vp=vp*2/np^2/f(1);
pp=pp*2/np^2/f(1);

pup=real((pf.*conj(uf))*2/np^2/f(1)); % scaled cross-spectra
pup=pup(2:(np/2+1),:);    %limit to the same frequencies as power spectra
pvp=real((pf.*conj(vf))*2/np^2/f(1));
pvp=pvp(2:(np/2+1),:);
puv=real((uf.*conj(vf))*2/np^2/f(1));
puv=puv(2:(np/2+1),:);

[F, Cuu] = logavg(f, up, nF);  % average into log bands
[F, Cvv] = logavg(f, vp, nF);
[F, Cpp] = logavg(f, pp, nF);
[F, Cpu] = logavg(f, pup, nF);
[F, Cpv, dF, Ns, Ne] = logavg(f, pvp, nF);
[F, Cuv] = logavg(f, puv, nF);
DOF=2*(Ne-Ns+1);

aa=find(F>lf);            % low frequency cutoff
lF=length(aa);            % number of frequencies we keep

F=F(aa);
dF=dF(aa);
DOF=DOF(aa);

Cuu=Cuu(aa,:);
Cvv=Cvv(aa,:);
Cpp=Cpp(aa,:);
Cpu=Cpu(aa,:);
Cpv=Cpv(aa,:);
Cuv=Cuv(aa,:);

mp=ones(lF,1)*mean(p);        % vertical scaling

F2=F*ones(1,nt);
ka=wavek(F2,mp+hp);

zp = -mp;
zv = -mp+hv;
h = mp+hp;
% sinhkh=sinh(ka.*(hp+mp));        %sinh scaling for water depth
% coshkh=cosh(ka.*(hp+mp));        %cosh scaling for water depth
% coshkhz=cosh(ka.*mp);            %scaling for pressure sensor elevation
% coshkhZ=cosh(ka.*(hp+hv));       %scaling for velocity elevation
sinhkh=sinh(ka.*h);        %sinh scaling for water depth
coshkh=cosh(ka.*h);        %cosh scaling for water depth
coshkhz=cosh(ka.*zp);            %scaling for pressure sensor elevation
coshkhZ=cosh(ka.*zv);       %scaling for velocity elevation

ac=find((coshkh.^2)>maxfac);
Su=(Cuu+Cvv).*(sinhkh./(2*pi*F2)./coshkhZ).^2;
Sp=Cpp.*(coshkh./coshkhz).^2;
ad=union(find(Sp<minspec),ac);




Dir=57.296*atan2(Cpu,Cpv)+Ndir;
Dir=mod(Dir+180,360)-180;

R2=((Cuu - Cvv).^2 + 4*Cuv.^2).^.5./(Cuu+Cvv);
Spread = 57.296*((1-R2)/2).^.5;

Su(ac)=nan;
Sp(ac)=nan;
Dir(ad)=nan;
Spread(ad)=nan;
end

%------------------------- LOGAVG -----------------------------------------

%function [ F, S] = logavg( f, s, n )
function [ F, S, dF, Ns, Ne ] = logavg( f, s, n )
%
%  [ F, S, dF, Ns, Ne ] = logavg( f, s, n )
%
%  Logarithmically averages the input spectrum, s, defined over
%  frequencies f into n uniformly-spaced "log10" frequency bands.
%  The averaged spectrum and frequencies are returned in S and F,
%
%  dF holds the bandwidth of each band
%  Ns and Ne contain the start and end indices of each band
%
%  If s is a matrix, averaging is performed over each column.
%
%  Copyright (C) 2001, Lee Gordon, NortekUSA LLC

% data must be in columns
f=f(:);
[nc,nt]=size(s);
if(nc~=length(f))
    'column length must be same as length of f'
    return;
end;

% some preliminaries

lf = log( f ) ;
dlf = 1.000000001*(lf(length(f)) - lf(1)) / n ;  % log frequency increment
NDX = 1 + floor((lf-lf(1))/dlf);
AA=[find(diff(NDX)>0)' length(f)];    % array of transitions plus final f

Cs=cumsum(s);
Cf=cumsum(f);
F=[Cf(AA(1)) diff(Cf(AA)')]./[AA(1) diff(AA)];
S=[Cs(AA(1),:); diff(Cs(AA,:))]./([AA(1) diff(AA)]'*ones(1,nt));

F=F';
dF=[AA(1) diff(AA)]*(f(10)-f(9));
AA=AA';
dF=dF';

Ns=[1 AA(1:length(AA)-1)'+1]';    %starting positions of each band
Ne=AA;                          %ending positions of each band
end
%------------------------------------------------------------------MAKE2DCD
function [S2DC,Xi,Yi,ZW,theta] = make2dcd(EF,FR,MD,SPR,resdir);
%constants
twopi=2*pi;

%nd is number of directions;
nd=360/resdir;

%step of directions
dtheta=twopi/nd;

%all angles
theta = 0:dtheta:twopi-dtheta;

%spreading per frequency in radians
Srad=SPR.*(pi/180);

%calculating the s-parameters (Kuik, Van Vledder, Holthuisen: a method for the routine
%analysis of pitch-and-roll buoy wave data: journal of physical oceanography, volume 18,
%eq. (38) p. 1027)
s=2./(Srad.^2)-1;
%for all directions

for ifr=1:length(EF)  %for all frequencies

    for id=1:nd    %for all directions
        dth=theta(id)-(MD(ifr)*pi/180);
        %directional energy distribution per frequency = E(f,theta) is calculated
        %this is also called the 2D-energy spectrum
        %(Kuik, Van Vledder, Holthuisen: a method for the routine
        %analysis of pitch-and-roll buoy wave data: journal of physical oceanography, volume 18,
        %eq. (2) p. 1021)
        Dspr(id,ifr)= abs((cos(dth/2))^(2*s(ifr)));
    end
    %normalizing factor is the frequency spectrum E(f)
    %an integration per frequency over all directions of E(f,theta) is done
    Norfact(ifr)= sum(Dspr(:,ifr))*twopi/nd;

    for id=1:nd
        %normalized directional energy distribution per frequency = Df(theta) is calculated
        NDspr(id,ifr)=Dspr(id,ifr)/Norfact(ifr);
    end

    for id=1:nd
        %energy per direction and frequency is calculated
        S2DC(id,ifr)=EF(ifr)*NDspr(id,ifr);
    end

end

%polar grid definition
Xspec = sin(theta)'*FR';
Yspec = cos(theta)'*FR';    %all frequencies are projected in all directions on x- and y-axis
gridx = -0.5:0.01:0.5;      %frequency domain
gridy = -0.5:0.01:0.5;

%2 matrices (21 on 21) are generated: Xi with gridx as rows and Yi with gridy as colums
[Xi,Yi]= meshgrid(gridx,gridy);

%S2DC=double(S2DC); % ESO
ZW = griddata(Xspec,Yspec,S2DC,Xi,Yi);
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function K=wavek(F,H);
% function K=wavek(F,H);
% where K is wavenumber (rad/m)
%       F is frequency (Hz)
%       H is depth (m)
%
% Copyright (C) 2001, Lee Gordon, NortekUSA LLC

g=9.80171;

% This routine use an approximate equation, then sharpens the result with
% one interpolation. The result is good to around 1 part in 10^-6.

% The equation came out of a textbook, but I have long since forgotton
% which one. If you know, please tell me! lgordon@nortekusa.com

e1=4*pi^2*F.^2.*H/g;        %f4 = omega^2 * h1/g
e2 = 1+0.6666666*e1 + 0.355555555*e1.^2 + 0.1608465608*e1.^3 + ...
    0.0632098765*e1.^4 + 0.0217540484*e1.^5 + 0.0065407983*e1.^6;
e3 = +e1.^2 + e1./e2;
K1=sqrt(e3)./H;

%compute error as basis for interpolation

o1=sqrt(g*K1.*tanh(K1.*H));
e1=o1.^2.*H/g;
e2 = 1+0.6666666*e1 + 0.355555555*e1.^2 + 0.1608465608*e1.^3 + ...
    0.0632098765*e1.^4 + 0.0217540484*e1.^5 + 0.0065407983*e1.^6;
e3 = +e1.^2 + e1./e2;
K2=sqrt(e3)./H;

%interpolate
K=2*K1-K2;

end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function info = waveinfo(Su, Sp, Dir, Spread, F, dF, ka, p, hp)
% WAVEINFO Calcula parámetros integrales del oleaje, incluyendo deriva de Stokes

% Crear índice lógico para datos válidos (sin NaN)
valid = ~(isnan(Su) | isnan(Sp) | isnan(Dir) | isnan(Spread) | ...
    isnan(F) | isnan(dF) | isnan(ka));

% Aplicar máscara de datos válidos
Su = Su(valid);
Sp = Sp(valid);
Dir = Dir(valid);
Spread = Spread(valid);
F = F(valid);
dF = dF(valid);
ka = ka(valid);

% Remover cualquier NaN residual (por precaución)
Su(isnan(Su)) = 0;
Sp(isnan(Sp)) = 0;
Dir(isnan(Dir)) = 0;
Spread(isnan(Spread)) = 0;
F(isnan(F)) = 0;
dF(isnan(dF)) = 0;
ka(isnan(ka)) = 1e-6;  % evitar división por cero en sinh

Wi = dF(:);
omega = 2 * pi * F;

% Momentos espectrales
m0u = sum(Su .* Wi);
m0p = sum(Sp .* Wi);
m1u = sum(F .* Su .* Wi);
m1p = sum(F .* Sp .* Wi);

% Altura significativa
info.Hs_u = 4 * sqrt(m0u);
info.Hs_p = 4 * sqrt(m0p);

% Frecuencias características
[~, idx_p] = max(Sp);
[~, idx_u] = max(Su);
info.fp_p = F(idx_p);
info.fp_u = F(idx_u);
info.fm_p = m1p / m0p;
info.fm_u = m1u / m0u;

% Dirección media, dispersión media y dirección pico
info.theta_bar = sum(Dir .* Sp .* Wi) / m0p;
info.spread_bar = sum(Spread .* Sp .* Wi) / m0p;
info.theta_p = Dir(idx_p);

% Profundidad total
h = mean(p) + hp;

% Deriva de Stokes superficial
% info.us0 = 2 * sum(omega .* Su .* Wi);
% info.us0 = 2 * sum(omega .* Su .* (cosh(2*ka*h) ./ sinh(2*ka*h)) .* Wi);
kh = ka * h;
info.us0 = 2 * sum(omega .* ka .* Su .* (cosh(2*kh) ./ (2 * sinh(kh).^2)) .* Wi);% LIBRO ARDHUIN

% Perfil vertical de deriva de Stokes
ka_vec = ka(:);
zvec = linspace(-h, 0, 100);
zvec(end) = 0;
info.zvec = zvec;

us_profile = zeros(size(zvec));
for k = 1:length(zvec)
    z = zvec(k);
    kh = ka_vec * h;
    factor = cosh(2 * ka_vec .* (z + h)) ./ (2 * sinh(kh).^2);
    us_profile(k) = 2 * sum(omega .* ka_vec .* Su .* factor .* Wi);
end
info.us_profile = us_profile;
info.h = h;
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function [us_ztheta, uz, vz, zvec, theta, valid_count] = radial_stokes_drift_profile(S2DC, F, dF, ka, h, theta, zvec)
% RADIAL_STOKES_DRIFT_PROFILE - Perfil vectorial de deriva de Stokes u_s(z, theta)
%
% Calcula el perfil vertical y direccional de la deriva de Stokes, incluyendo sus
% componentes horizontales proyectadas u (este-oeste) y v (norte-sur).
%
% Entradas:
%   S2DC   - Espectro direccional [n_theta x n_freq]
%   F      - Frecuencias lineales [1 x n_freq]
%   dF     - Ancho de banda por frecuencia [1 x n_freq]
%   ka     - Número de onda por frecuencia [1 x n_freq]
%   h      - Profundidad total (m)
%   theta  - Ángulos (rad) [1 x n_theta]
%   zvec   - Vector de profundidades deseadas (negativas hacia abajo)
%
% Salidas:
%   us_ztheta   - Matriz [nz x n_theta] de deriva radial
%   uz          - Matriz [nz x n_theta] componente zonal
%   vz          - Matriz [nz x n_theta] componente meridional
%   zvec        - Vector de profundidades
%   theta       - Direcciones (radianes)
%   valid_count - Número de frecuencias válidas usadas en cada punto

omega = 2 * pi * F(:);
ka = ka(:);
dF = dF(:);

n_theta = size(S2DC, 1);
n_freq = length(F);
n_z = length(zvec);

us_ztheta = zeros(n_z, n_theta);
valid_count = zeros(n_z, n_theta);

kh = ka * h;
denominator = 2 * sinh(kh).^2;  % vector [n_freq x 1]

for k = 1:n_z
    z = zvec(k);
    numerator = cosh(2 * ka * (z + h));  % vector [n_freq x 1]
    factor = numerator ./ denominator;  % vector [n_freq x 1]
    for j = 1:n_theta
        Sj = S2DC(j, :).';          % espectro en theta_j [n_freq x 1]
        mask = ~isnan(Sj);
        valid_count(k, j) = sum(mask);
        if valid_count(k, j) > 0
            us_ztheta(k, j) = 2 * sum(omega(mask) .* ka(mask) .* Sj(mask) .* factor(mask) .* dF(mask));
        else
            us_ztheta(k, j) = NaN;
        end
    end
end

% Proyecciones horizontales
theta_row = reshape(theta, 1, []);      % [1 x n_theta]
uz = us_ztheta .* cos(theta_row);       % componente zonal (u)
vz = us_ztheta .* sin(theta_row);       % componente meridional (v)
end


