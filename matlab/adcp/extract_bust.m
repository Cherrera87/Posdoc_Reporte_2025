close all; clear all; clc;

% % Cierra cualquier pool existente
% if ~isempty(gcp('nocreate'))
%     delete(gcp('nocreate'));
% end
% 
% % Detecta número de núcleos disponibles
% numCores = feature('numcores')-1;
% 
% % Crea un nuevo parpool con el número máximo de workers
% c = parcluster('Processes');
% c.NumWorkers = numCores; % Puedes ajustar esto si deseas usar menos núcleos
% 
% % Guarda el perfil actualizado (opcional, si quieres que se mantenga en futuras sesiones)
% saveProfile(c);
% 
% % Crea el pool
% parpool('Processes', numCores);

for arch = 5
    tic
    if arch == 1 % Rutas PM 2024
        path_files = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/PM_signature/2025_04_01_PM_signature/level0/';
        path_save = '/media/root01/Pancho_SSD/ADCP-DATA/PM_2024_Burst/';
        id_caso = 'Sig_PM_';

    elseif arch == 2 % Rutas BSM 2022
        path_files = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/BSM_signature/2022_03_14_BSM_signature/level0/';
        % path_save = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/BSM_signature/2022_03_14_BSM_signature/Burst/';
        path_save = '/media/root01/Pancho_SSD/ADCP-DATA/BSM_2022_Burst/';
        id_caso = 'Sig_BSM_';
    elseif arch == 3 % Rutas PM 2022
        path_files = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/PM_signature/2022_03_24_PM_signature/level0/';
        % path_save = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/PM_signature/2022_03_24_PM_signature/Burst/';
        path_save = '/media/root01/Pancho_SSD/ADCP-DATA/PM_2022_Burst/';
        id_caso = 'Sig_PM_';
    elseif arch == 4 % Rutas PM 2025a
        path_files = '/media/root01/Pancho_SSD/ADCP-DATA/PM-abril-junio_2025/100599/converted/';
        % path_save = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/PM_signature/2022_03_24_PM_signature/Burst/';
        path_save = '/media/root01/Pancho_SSD/ADCP-DATA/PM_2025a_Burst/';
        id_caso = 'Sig_PM_';   
    elseif arch == 5 % Rutas bsm 2024
        path_files = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/BSM_signature/2025_10_21_BSM_signature(recuperado)/level0/';
        % path_save = '/run/user/1000/gvfs/smb-share:server=nas1-oleaje.cicese.mx,share=oleaje1/sensores_fijos/signature/BSM_signature/2022_03_14_BSM_signature/Burst/';
        path_save = '/media/root01/Pancho_SSD/ADCP-DATA/BSM_2024_Burst/';
        id_caso = 'Sig_BSM_';    
    end    

    files = dir([path_files '*.mat']);
    if arch == 3
        files(end) = [];
    end

    for ii = 1:length(files)
        S = struct();
       disp(['Procesando archivo: ',id_caso, ': ', files(ii).name])
        data = load([files(ii).folder '/' files(ii).name]);

        S.Config = data.Config;
        S.Descriptions = data.Descriptions;
        S.Units = data.Units;

        campos = fieldnames(data.Data);
        % Buscar los que empiezan con 'IBurst' y eliminarlos
        ind = regexp(campos, '^IBurst');  % Busca con expresión regular
        campos = campos(cellfun('isempty', ind));

        p0 = [find(diff(data.Data.Burst_EnsembleCount)<1)];
        p_ini =  [1; p0+1];
        p_end = [p0;length((data.Data.Burst_EnsembleCount))];
        for ij =1 :length(p_ini)

            p0 = [p_ini(ij):p_end(ij)];
            Data=[];
            for ik = 1:length(campos)
                nd = ndims(data.Data.(campos{ik}));
                subs = repmat({':'}, 1, nd); % Crea {:, :, :}
                subs{1} = p0;               % Solo queremos los primeros 10 de la primera dimensión

                S.Data.(campos{ik}) = data.Data.(campos{ik})(subs{:});
            end
                        
            fecha = datestr(S.Data.Burst_Time(1),'yyyymmddHHMM');
            save([path_save id_caso fecha],'-struct','S')
            fecha = [];     
        end
    end

end