# Procesamiento de oleaje y corrientes con ADCP y Radar HF

Este repositorio contiene los códigos MATLAB desarrollados en el marco del proyecto
posdoctoral para la estimación del oleaje direccional, la deriva de Stokes y la interpretación
de corrientes superficiales obtenidas con radar HF, usando mediciones de perfiladores
acústicos (ADCP) en la Bahía de Todos Santos, Ensenada, México.

Los métodos implementados están documentados en:

- Reporte Técnico 01: Procesamiento de datos de ADCP  
- Reporte Técnico 02: Análisis de corrientes con radar HF  

## Contenido

- `matlab/adcp/`  
  Segmentación de bursts, corrección de declinación magnética (IGRF),
  estimación espectral PUV, espectros direccionales y deriva de Stokes.

- `matlab/radarHF/`  
  Conversión de archivos CRAD del radar WERA, proyección radial de corrientes ADCP
  y comparación HF–ADCP.

- `data/sample/`  
  Datos mínimos de ejemplo para ejecutar los scripts.

## Flujo de trabajo

1. Organizar bursts de ADCP  
2. Corregir declinación magnética  
3. Estimar espectros y deriva de Stokes  
4. Convertir datos de radar HF  
5. Proyectar corrientes ADCP en la dirección radial  
6. Comparar radar HF vs ADCP  

## Requisitos

- MATLAB R2020 o superior
- Toolboxes:
  - Signal Processing
  - Parallel Computing (opcional)
  - IGRF / Aerospace Toolbox
  - matWERA (para radar HF)

## Ejemplo rápido

```matlab
addpath matlab/adcp
plot_cobertura_temporal_mediciones
