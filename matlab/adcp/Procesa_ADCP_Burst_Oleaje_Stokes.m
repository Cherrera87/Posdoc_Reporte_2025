function Procesa_ADCP_Burst_Oleaje_Stokes()
% Procesa_ADCP_Burst_Oleaje_Stokes
%
% DESCRIPCIÓN:
%   Procesa archivos ADCP en modo burst para estimar:
%     - Espectros de oleaje (u, v y presión)
%     - Espectro direccional E(f,θ)
%     - Parámetros integrales del oleaje
%     - Perfil vertical y direccional de la deriva de Stokes
%
%   El script es AUTOCONTENIDO e incluye todas las funciones necesarias.
%
% ENTRADAS:
%   El usuario debe editar la sección CONFIGURACIÓN.
%
% SALIDAS:
%   Archivo .mat con:
%     - dat  : estructuras espectrales e integrales
%     - Us3D : perfiles de deriva de Stokes
%
% REFERENCIAS:
%   Emery & Thomson (2001)
%   Gordon (2001), NortekUSA
%   Kuik et al. (1988)
%   Ardhuin et al. (2009, 2010)
%
% AUTOR:
%   Carlos F. Herrera Vázquez
%
% FECHA:
%   2025-12
% ======================================================================= %

%% ======================= CONFIGURACIÓN ======================= %%
cfg = struct();

cfg.path_burst  = 'PATH_A_CARPETA_BURST';   % <-- EDITAR
cfg.output_mat  = 'PATH_A_SALIDA/data_waves.mat';
cfg.case_id     = 'Sitio / Periodo';

cfg.min_duration_sec = 300;    % duración mínima válida
cfg.nF     = 100;              % bandas espectrales
cfg.hp     = 0.5;              % altura sensor presión (m)
cfg.resdir = 2;                % resolución direccional (deg)
cfg.cell_index = 3;            % celda de velocidad usada
cfg.verbose = true;

%% ======================= VALIDACIONES ======================== %%
assert(isfolder(cfg.path_burst),'La carpeta Burst no existe.');
files = dir(fullfile(cfg.path_burst,'*.mat'));
assert(~isempty(files),'No hay archivos .mat en la carpeta.');

m = matfile(cfg.output_mat,'Writable',true);
dat_all  = cell(numel(files),1);
Us3D_all = cell(numel(files),1);

%% ======================= LOOP PRINCIPAL ====================== %%
for ii = 1:numel(files)

    if cfg.verbose
        fprintf('[%s] %d/%d\n', cfg.case_id, ii, numel(files));
    end

    S = load(fullfile(files(ii).folder,files(ii).name));
    if ~isfield(S,'Data') || ~isfield(S,'Config'), continue; end

    Data   = S.Data;
    Config = S.Config;

    t0 = datetime(Data.Burst_Time(1),'ConvertFrom','datenum');
    t1 = datetime(Data.Burst_Time(end),'ConvertFrom','datenum');
    if seconds(t1-t0) < cfg.min_duration_sec, continue; end

    %% Variables básicas
    z = double(Config.Burst_BlankingDistance + ...
        (1:double(Config.Burst_NCells))*Config.Burst_CellSize);

    p  = double(Data.Burst_Pressure);
    fs = double(Config.Burst_SamplingRate);
    dt = 1/fs;

    U = double(Data.Burst_Velocity_ENU_true(:,1,:));
    V = double(Data.Burst_Velocity_ENU_true(:,2,:));

    ik = cfg.cell_index;
    u = detrend(U(:,ik));
    v = detrend(V(:,ik));

    if mod(length(u),2), u(end)=[]; v(end)=[]; p(end)=[]; end

    hv = z(ik);
    parms = [0.03 200 0.003 0];

    %% Espectro
    [dat.Su,dat.Sp,dat.Dir,dat.Spread,dat.F,dat.dF,dat.DOF,ka] = ...
        wds(u,v,p,dt,cfg.nF,cfg.hp,hv,parms);

    if all(isnan(dat.Sp)), continue; end

    %% Espectro direccional
    [dat.S2DC,dat.Xi,dat.Yi,dat.ZW,dat.theta] = ...
        make2dcd(dat.Sp,dat.F,dat.Dir,dat.Spread,cfg.resdir);
    dat.S2DC = dat.S2DC*pi/180;

    %% Parámetros integrales
    dat.info = waveinfo(dat.Su,dat.Sp,dat.Dir,dat.Spread,...
                        dat.F,dat.dF,ka,p,cfg.hp);

    %% Deriva de Stokes
    [Us3D.us_ztheta,Us3D.uz,Us3D.vz,...
     Us3D.zvec,Us3D.theta,Us3D.valid_count] = ...
        radial_stokes_drift_profile(dat.S2DC,dat.F,dat.dF,ka,...
                                    dat.info.h,dat.theta,dat.info.zvec);

    dat.time = t0;

    dat_all{ii}  = dat;
    Us3D_all{ii} = Us3D;
end

m.dat  = dat_all;
m.Us3D = Us3D_all;

if cfg.verbose
    disp('Procesamiento finalizado.');
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% === Espectros mediante spectf (promedio en frecuencia) ===
% Queremos ~nF bandas de salida, así que elegimos Nfa ≈ (#bins)/(nF)
Nraw = np/2;                             % nº de bins de frecuencia (positivos)
Nfa  = 2;%max(1, floor(Nraw / nF));         % nº de bins elementales por banda
a0   = 0;                                % no saltar bandas bajas

% Auto-espectros
Suu = spectf(u, dt, Nfa, a0);   % [f Hz, Sxx (m^2/s^2 * s)]
Svv = spectf(v, dt, Nfa, a0);
Spp = spectf(p, dt, Nfa, a0);

% Cross-espectros: S(:,2)=Sxx, S(:,3)=Syy, S(:,4)=Sxy, S(:,5)=phase, S(:,6)=coh
Spu = spectf(p, u, dt, Nfa, a0);   % P→x, u→y
Spv = spectf(p, v, dt, Nfa, a0);
Suv = spectf(u, v, dt, Nfa, a0);

% Extraemos frecuencias y densidades
F   = Suu(:,1);       % Hz
Cuu = Suu(:,2);       % m^2/s^2 per Hz
Cvv = Svv(:,2);
Cpp = Spp(:,2);

Cpu = real(Spu(:,4)); % Re{S_pu}
Cpv = real(Spv(:,4)); % Re{S_pv}
Cuv = real(Suv(:,4)); % Re{S_uv}

% Ancho de banda aproximado (dF) y DOF
dF  = [diff(F); F(end)-F(end-1)];       % Hz aprox
DOF = 2 * Nfa * ones(size(F));         % ~2*Nfa grados de libertad por banda


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

function S = spectf(x,y,dt,Nfa,a0)
%        S = SPECTf(x,dt,Nfa)
%        S = SPECTf(x,y,dt,Nfa)   cross spectrum
% 
% Frequency averaged power spectrum estimate,  GEOPHYSICAL NORMALIZATION
% Trend is removed, Blackman-Harris window is used. K.K.Kahma 1990-05-19
%
%     x , y  = data vectors
%     dt = sampling interval in seconds
%     Nfa = number of elementary frequency bands which are averaged
%
%     S(:,1) = f      (1/second == Hz)
%     S(:,2) = Sxx    (unit*unit*second)
%
% If cross spectrum is calculated 
%     S(:,3) = Syy
%     S(:,4) = Sxy
%     S(:,5) = phase angle = 180/pi*atan2(-imag(Sxy),real(Sxy))
%     S(:,6) = coherence   = abs(Sxy./sqrt(Sxx.*Syy))
%
%     positive phase means x leads y

%        S = SPECTF(x,y,dt,Nfa,a0)
% Elementary frequency bands 0:a0-1 (matlab index 1:a0) are ignored. 
% Default a0 is 0, i.e. all bands including zero (mean value) are incuded.
% 

x = x(:).';      	% Make sure x is a row vector
N  = max(size(x));      % Number of data points
window=blackhar(N).';

if max(size(y)) ~= N,
   if (max(size(y)) == 1) | (nargin < 5)

% ***************
   % Spectrum 
% ***************

   if (nargin < 4), Nfa = 0; end    % default a0
   if (nargin < 3), dt = 31; end    % default Nfa
   a0 = Nfa; Nfa = dt; dt = y; 

   Nfft=0; maxb=0; C=0; df=0;         % To define these variables before Xx
   Xx = fft(window.*detrend(x));
   Nfft = length(Xx);                 % Number of points in FFT
   maxb = Nfft/2+1;
   Xx(maxb+1:Nfft)=[];
   Xx(maxb) = Xx(maxb)/2;

   C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
   df = 2*pi/(dt*Nfft);

   if Nfa==1
      f = [a0:maxb-1]*df;
      Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
   else
      if Nfa > 20
%       When Nfa is large enough this is as fast as vectorized
%       averaging and it requires far less memory    
        m=0; a=a0+1; b=a0+Nfa;
        while b <= maxb
           m=m+1;
           Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
           f(m) = df*((a+b-2)/2);
           a=a+Nfa; b=b+Nfa;
        end
      else
        m=fix((maxb-a0) / Nfa);
        f=([1:m]*Nfa+(a0-0.5-Nfa/2))*df;
        b=a0+m*Nfa;

%        Old bin averaging loop
%        sx=zeros(m,Nfa);
%        for i=1:Nfa 
%           sx(:,i) = abs(Xx(a0+i:Nfa:b)).^2;
%        end
%        Pxx=(sum(sx.')*C);

        sx=zeros(Nfa,m);
        sx(:) = abs(Xx(a0+1:b)).^2;  
        Pxx=(sum(sx)*C);
      end
      a=a0+1+m*Nfa;
      if a <= maxb
         m=m+1;
         c = maxb+1-a; 
         Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
         f(m) = df*(a+maxb-2)/2;
      end
    end 
    clear Xx window
    S = [f/2/pi;2*pi*Pxx].';
 
  else

  error('x and y are not of same size'); end

else

% **********************
   % Cross spectrum
% **********************

   if (nargin < 5), a0 = 0; end    % default a0
   if (nargin < 4), Nfa = 31; end  % default Nfa

   y = y(:).';
   Nfft=0; maxb=0; C=0; df=0;
   Xx = fft(window.*detrend(x));
   Nfft = length(Xx);                 % Number of points in FFT
   maxb = Nfft/2+1;
   Xx(maxb+1:Nfft)=[];
   Xx(maxb) = Xx(maxb)/2;

   C = dt/(Nfa*pi*norm(window)^2);    % Scaling coefficient
   df = 2*pi/(dt*Nfft);

   Yy = fft(window.*detrend(y));
   Yy(maxb) = Yy(maxb)/2;
   Yy(maxb+1:Nfft)=[];

   if Nfa==1
      f = [a0:maxb-1]*df;
      Pxx = (abs(Xx(a0+1:maxb)).^2)*C;
      Pyy = (abs(Yy(a0+1:maxb)).^2)*C;
      Pxy = (conj(Xx(a0+1:maxb)).*Yy(a0+1:maxb))*C;
   else
      if Nfa > 20
         m=0; a=a0+1; b=a0+Nfa;
         while b <= maxb
            m=m+1;
            Pxx(m) = sum(abs(Xx(a:b)).^2)*C;
            Pyy(m) = sum(abs(Yy(a:b)).^2)*C;
            Pxy(m) = sum(conj(Xx(a:b)).*Yy(a:b))*C;
            f(m) = df*((a+b-2)/2);
            a=a+Nfa; b=b+Nfa;
         end
      else
         m=fix((maxb-a0) / Nfa);
         f=([1:m]*Nfa+(a0-0.5-Nfa/2))*df;
         b=a0+m*Nfa;
%         sx=zeros(m,Nfa);
%         for i=1:Nfa 
%            sx(:,i)  = abs(Xx(a0+i:Nfa:b)).^2;
%            sy(:,i)  = abs(Yy(a0+i:Nfa:b)).^2;
%            sxy(:,i) = conj(Xx(a0+i:Nfa:b)).*Yy(a0+i:Nfa:b);
%         end
         sx=zeros(Nfa,m);
         sx(:) = abs(Xx(a0+1:b)).^2;  
         Pxx=(sum(sx)*C);
         sx(:) = abs(Yy(a0+1:b)).^2;  
         Pyy=(sum(sx)*C);
         sx(:) = conj(Xx(a0+1:b)).*Yy(a0+1:b);  
         Pxy=(sum(sx)*C);
         a=a0+1+m*Nfa;
      end
   
      if a <= maxb
         m=m+1; 
         c = maxb+1-a; 
         Pxx(m) = sum(abs(Xx(a:maxb)).^2)*C*Nfa/c;
         Pyy(m) = sum(abs(Yy(a:maxb)).^2)*C*Nfa/c;
         Pxy(m) = sum(conj(Xx(a:maxb)).*Yy(a:maxb))*C*Nfa/c;
         f(m) = df*(a+maxb-2)/2;
      end
   end
   phase = 180/pi*atan2(-imag(Pxy),real(Pxy));
   coh   = abs(Pxy./sqrt(Pxx.*Pyy));
   clear Xx Yy window sx sy sxy
   S = [f/2/pi;2*pi*Pxx;2*pi*Pyy;2*pi*Pxy;phase;coh].';
end
end

function w = blackhar(n)
%
%  BLACKHAR(N) returns the N-point Blackman-Harris window as a column vector
%  K.Kahma 1989-07-20

m = (0:n-1)' * ( 2*pi/(n-1) ) ;

w = (.35875 - .48829*cos(m) + .14128*cos(2*m) - 0.01168*cos(3*m) ) ;
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
        %Dspr(id,ifr) = dir_von_mises(theta(id), MD(ifr)*pi/180, SPR(ifr));
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

function D = dir_von_mises(theta, theta0, spr_deg)
% Distribución direccional Von Mises
% theta, theta0 en radianes
% spr_deg = desviación angular σθ en grados

% Convertir spread a radianes
sigma = spr_deg * pi/180;

% Concentración
kappa = 1/(sigma^2 + eps);

% Normalización
I0 = besseli(0, kappa);

% Distribución
D = exp(kappa*cos(theta - theta0)) / (2*pi*I0);
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
