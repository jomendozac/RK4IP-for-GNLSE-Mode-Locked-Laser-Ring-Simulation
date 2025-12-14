% SMLL.m
% Programa Principal para Oscilador Paramétrico Óptico basado en Fibra / Láser de Bloqueo de Modos
% Algoritmo usado: RK4IP (Método de Fourier de paso dividido de 4to orden)

close all;  % Cierra todas las ventanas de gráficos que estén abiertas para empezar limpio.
clear all;  % Borra todas las variables de la memoria RAM para evitar conflictos con datos viejos.
clc;        % Limpia la consola de comandos de MATLAB (la pantalla de texto).

% --- GESTIÓN DE CARPETAS ---
% Verifica si existe una carpeta llamada 'L100x120x100_10dB' en el directorio actual.
if ~exist('outputring', 'dir')
    mkdir('outputring'); % Si no existe, la crea (mkdir = make directory). Aquí se guardarán las fotos.
end

% --- CONSTANTES FÍSICAS Y DE LA LUZ ---
c0 = 3e8;               % Velocidad de la luz en el vacío (300,000,000 m/s).
lambdaP = 1550e-9;      % Longitud de onda central del láser (1550 nanómetros, estándar en telecomunicaciones).
wP = 2*pi*c0/lambdaP;   % Frecuencia angular de la luz (omega). Es el "color" de la luz en radianes/segundo.

% --- CONFIGURACIÓN DEL TIEMPO (GRILLA TEMPORAL) ---
Tmax = 120e-12;          % Define el tamaño de la ventana de simulación (60 picosegundos hacia cada lado).
Tmin = -Tmax;           % El tiempo empieza en -60 ps.
tol = 1e-3;             % Tolerancia de error para el cálculo matemático (RK4IP). Define qué tan preciso queremos ser.

global t N f            % Define variables globales para que las otras funciones (.m) puedan usarlas.
N = 2^17;               % Número de puntos de la simulación. 2^17 = 131,072 puntos.
                        % Se usa una potencia de 2 para que la Transformada Rápida de Fourier (FFT) sea muy veloz.

deltat = 2*Tmax/(N-1);  % Resolución temporal: cuánto tiempo pasa entre un punto y el siguiente.
t = linspace(Tmin, Tmax-deltat, N); % Crea el vector de tiempo: una lista de números desde Tmin hasta Tmax.

% --- CONFIGURACIÓN DE LA FRECUENCIA (GRILLA ESPECTRAL) ---
deltaf = 1/(2*Tmax);    % Resolución de frecuencia (inversa del tiempo total).
% Crea el vector de frecuencias ordenado como lo necesita MATLAB para la FFT:
% Primero las frecuencias positivas (0 a max) y luego las negativas (-max a casi 0).
f = [0:deltaf:(N/2-1)*deltaf, -N*deltaf/2:deltaf:-deltaf]; 
w = 2*pi*f;             % Convierte frecuencia (Hz) a frecuencia angular (rad/s).
wshift = -w + wP;       % Ajusta la frecuencia para que esté centrada respecto a la frecuencia del láser.

% --- PREPARACIÓN PARA GRAFICAR ---
wl = 1e9*2*pi*c0./wshift; % Convierte las frecuencias a Longitud de Onda (nm) para los ejes de las gráficas.
tm = 1e12*t;              % Convierte el tiempo a Picosegundos (ps) para que los números en los ejes sean legibles.

% --- PARÁMETROS DEL EFECTO RAMAN (Dispersión inelástica) ---
% El efecto Raman transfiere energía de frecuencias altas a bajas dentro del pulso.
% fR = 0.18;              % Valor típico para sílice.
fR = 0;                 % Valor actual: 0.0. Significa que el efecto Raman está APAGADO en esta simulación.
t1 = 12.2e-15;            % Tiempo característico de la vibración molecular 1.
t2 = 32e-15;              % Tiempo característico de la vibración molecular 2.
t_shock = 1/wP;           % Parámetro de "auto-inclinación" (self-steepening) del pulso de la GNLSE.
%este tshock es el mismo 1/w0 de la GNLSE y que esta computado dentro de
%N_op

tr = t - t(1);            % Crea un vector de tiempo desplazado que empieza exactamente en 0.
% Fórmula de la respuesta de impulso Raman (Modelo de Blow y Wood) directamente de Agrawal 2.3.40:
hR = (t1^2+t2^2)/t1/t2^2 * exp(-tr/t2) .* sin(tr/t1);
hR = hR ./ trapz(tr, hR); % Normaliza la curva para que el área bajo la curva sea 1 (conservación de energía).
hR_f = fft(hR);           % Convierte la respuesta Raman al dominio de la frecuencia (usando FFT).

% tt = 50e-15; % (Comentado) Parámetro de dispersión de ganancia antiguo.

% --- CONFIGURACIÓN DE LA FIBRA DOPADA (EDF - Medio Activo) ---
% Esta es la fibra que amplifica la luz (le da energía).
n2 = 2.6e-20;               % Índice de refracción no lineal (Efecto Kerr: la intensidad de luz cambia la velocidad de la luz).
Esat = 50e-12;              % Energía de Saturación (50 pJ). Límite de energía que la fibra puede dar.
L_edf = 0.10;               % Longitud de la fibra dopada: 20 centímetros.
AdB_edf = 10.0;             % Pérdida básica de la fibra en dB por metro.
MFD_edf = 8.0e-6;           % Diámetro del campo modal (grosor del haz de luz en la fibra).
alfa_edf = AdB_edf/4.343;   % Convierte la pérdida de dB/m a unidades lineales (Nepers/m).
Aeff_edf = pi*(MFD_edf/2)^2;% Área efectiva del núcleo de la fibra.
gamma_edf = 2*pi*n2/lambdaP/Aeff_edf; % Coeficiente de no linealidad (Gamma). Qué tan fuerte es el efecto no lineal.

delta_wl = 40e-9;           % Ancho de banda del amplificador (40 nm). Rango de colores que puede amplificar.
Omega = (2*pi*c0/lambdaP^2) * delta_wl; % Ancho de banda en frecuencia angular.
tt = 1/Omega;               % Tiempo relacionado con la respuesta de la ganancia.

D_edf = 10.0e6;             % Dispersión de la fibra dopada. Define cuánto se "desparrama" el pulso en el tiempo.
% Calcula beta2 (dispersión de segundo orden) a partir de D:
b2_edf = lambdaP^2/(2*pi*c0*1e-12) * (D_edf); 
beta_edf = [0 b2_edf];      % Vector que guarda los coeficientes de dispersión.

% Ajuste de unidades de beta (para que coincidan con la escala temporal):
for ii = 1:length(beta_edf)
    beta_edf(ii) = (1e-12)^(ii) * beta_edf(ii); % Multiplica por 10^-12 según el orden.
end

% --- CONFIGURACIÓN DE LA FIBRA MONOMODO (SMF - Pasiva) ---
% Fibra estándar de transporte (sin ganancia).
L_smf1 = 0.4;               % Longitud tramo 1 (30 cm).
L_smf2 = 0.4;               % Longitud tramo 2 (30 cm).
L_smf3 = 0.4;               % Longitud tramo 3 (30 cm).
AdB_smf = 0.2;              % Pérdida muy baja (0.2 dB/km).
alfa_smf = (AdB_smf/1000)/4.343; % CONVIERTE pérdida a unidades lineales por metro.
MFD_smf = 10.8e-6;          % Diámetro del núcleo un poco más grande que la EDF.
Aeff_smf = pi*(MFD_smf/2)^2;% Área efectiva.
gamma_smf = 2*pi*n2/lambdaP/Aeff_smf; % Gamma para la SMF.
D_smf = 17.0e6;             % Dispersión típica de fibra estándar (positiva/anómala).
b2_smf = lambdaP^2/(2*pi*c0*1e-12) * (D_smf); % Cálculo de beta2.
beta_smf = [0 b2_smf];

for ii = 1:length(beta_smf)
    beta_smf(ii) = (1e-12)^(ii) * beta_smf(ii); % Ajuste de unidades.
end

% --- CONFIGURACIÓN DE LA FIBRA COMPENSADORA (DCF) ---
% Fibra especial para corregir la dispersión de las otras fibras.
L_dcf = 0.1;                % Longitud (30 cm).
AdB_dcf = 1.0;              % Pérdida más alta (1 dB/km).
alfa_dcf = (AdB_dcf/1000)/4.343;
MFD_dcf = 6.0e-6;           % Núcleo muy pequeño -> Alta no linealidad.
Aeff_dcf = pi*(MFD_dcf/2)^2;
gamma_dcf = 2*pi*n2/lambdaP/Aeff_dcf; %definicion del libro ecuacion antes de 4,27
D_dcf = -40e6;              % Dispersión NEGATIVA (Normal) fuerte para cancelar la positiva de las otras.
b2_dcf = lambdaP^2/(2*pi*c0*1e-12) * (D_dcf); %relaciona dispersion con el beta en 4,16
beta_dcf = [0 b2_dcf];

for ii = 1:length(beta_dcf)
    beta_dcf(ii) = (1e-12)^(ii) * beta_dcf(ii); % Ajuste de unidades.
end

% --- INICIALIZACIÓN DE LA SEÑAL (RUIDO) ---
% El láser empieza desde ruido aleatorio (emisión espontánea).
phase = pi*rand([1,N]);     % Genera fases aleatorias entre 0 y pi.
Unoise = 1e-4 .* randn(1,N) .* exp(-1i*phase); % Crea ruido blanco de muy baja potencia.
Unoise_f = fft(Unoise);     % Calcula el espectro (frecuencias) del ruido.
It0 = abs(Unoise).^2;       % Calcula la intensidad en el tiempo (Potencia inicial).
If0 = abs(Unoise_f).^2;     % Calcula la intensidad en frecuencia (Espectro inicial).
Uin0 = Unoise;              % Establece este ruido como la entrada inicial de la simulación.

w2 = wl((wl>1500) & (wl<1600)); % Define un rango de visualización entre 1500 y 1600 nm.

% --- PARÁMETROS DEL ABSORBEDOR SATURABLE (SA) ---
% El SA actúa como una puerta que se abre solo para pulsos intensos.
a0 = 0.05;      % Profundidad de modulación (cuánto absorbe como máximo).
ac = 0.55;      % Pérdidas no saturables (pérdida constante).
Ts = 750e-15;   % Tiempo de recuperación (qué tan rápido se "cierra" la puerta).
Es = 1e-12;     % Energía de saturación.
% Llama a la función SA una vez para ver el estado inicial (asegúrate de tener el archivo SA.m).
alfa0 = SA(Uin0, a0, ac, Ts, Es); 

% --- PARÁMETROS DE LA CAVIDAD (ANILLO) ---
FB = 70/100;    % Feedback: El 70% de la luz vuelve a entrar al anillo.
CL1 = 96/100;   % Pérdida por acople/conexión 1 (Eficiencia del 96%).
CL2 = 96/100;   % Pérdida por acople/conexión 2.
GdB = 20;       % Ganancia objetivo en dB (aprox 20 dB es una ganancia alta para facilitar el arranque).
gain = GdB/L_edf/4.343; % Calcula cuánta ganancia por metro necesita la fibra dopada.

indc = 0;       % Contador de vueltas de la simulación.
indx = 0;       % Contador para nombrar los archivos de imagen guardados.
loopmax = 300;   % Número máximo de vueltas (loops) que dará la simulación por el anillo

% --- BUCLE PRINCIPAL (MAIN LOOP) ---
% Aquí empieza la simulación física, vuelta tras vuelta.
fprintf('Starting Simulation...\n');

for loop = 1:loopmax
    indc = indc + 1; % Incrementa el número de vuelta actual.
    
    % --- 1. PROPAGACIÓN POR FIBRA DOPADA (EDF) ---
    %podemos usar los archivos propagation EDF o EDFL en CASO DE...
    % La señal entra en la fibra que le da energía (amplificación).
    Uin_edf = Uin0;
    % Llama a la función matemática que resuelve la ecuación de Schrödinger en la fibra dopada.
    % Se puede usar la rutina %Propagation_EDF tambien
    Uout_edf = Propagation_EDFL(Uin_edf, tol, L_edf, gain, Esat, alfa_edf, ...
                                beta_edf, gamma_edf, t_shock, fR, hR_f);
                            
    % --- 2. PROPAGACIÓN POR FIBRA SMF 1 ---
    % La luz sale de la EDF y pasa por un tramo de fibra normal.
    Uin_smf1 = sqrt(CL1) * Uout_edf; % Aplica la pérdida del conector (CL1).
    Uout_smf1 = Propagation_SMF(Uin_smf1, tol, L_smf1, alfa_smf, beta_smf, ...
                                     gamma_smf, t_shock, fR, hR_f);
                                 
    % --- 3. ABSORBEDOR SATURABLE (SA) ---
    % La luz llega al dispositivo que "limpia" el pulso, eliminando el ruido de baja intensidad.
    alfa = SA(Uout_smf1, a0, ac, Ts, Es); % Calcula cuánto debe absorber en cada instante.
    Uout_sa = sqrt(1 - alfa) .* Uout_smf1; % Aplica la absorción a la señal.
    
    % --- 4. PROPAGACIÓN POR FIBRA SMF 2 ---
    % Otro tramo de fibra normal después del SA.
    Uin_smf2 = sqrt(CL2) * Uout_sa; % Aplica pérdida de conector (CL2).
    Uout_smf2 = Propagation_SMF(Uin_smf2, tol, L_smf2, alfa_smf, beta_smf, ...
                                     gamma_smf, t_shock, fR, hR_f);
                                 
    % --- SALIDA DEL ACOPLADOR (OUTPUT COUPLER) ---
    % Aquí se divide la luz: una parte sale (output) y otra sigue en el anillo.
    Uout_loop = sqrt(1-FB) * Uout_smf2; % Esta es la parte que SALE (30%).
    Uouf_loop = fft(Uout_loop);         % Calcula el espectro de la salida.
    Iout = abs(Uout_loop).^2;           % Intensidad temporal de salida.
    Iouf = 10*log10(abs(Uouf_loop).^2); % Intensidad espectral en decibelios (dB).
    
    % --- 5. PROPAGACIÓN POR FIBRA SMF 3 (Feedback path) ---
    % La luz que SE QUEDA en el anillo (70%) viaja por otro tramo de fibra.
    Uin_smf3 = sqrt(FB) * Uout_smf2;    % Aplica el factor de Feedback (FB).
    Uout_smf3 = Propagation_SMF(Uin_smf3, tol, L_smf3, alfa_smf, beta_smf, ...
                                     gamma_smf, t_shock, fR, hR_f);
                                 
    % --- 6. PROPAGACIÓN POR FIBRA DCF (Compensación) ---
    % Fibra final para corregir la dispersión acumulada antes de empezar de nuevo.
    Uin_dcf = sqrt(CL1) * Uout_smf3;
    Uout_dcf = Propagation_SMF(Uin_dcf, tol, L_dcf, alfa_dcf, beta_dcf, ...
                                    gamma_dcf, t_shock, fR, hR_f);
    
    % La salida de la DCF se convierte en la entrada de la vuelta siguiente (cierre del bucle).
    Uin0 = sqrt(CL1)*Uout_dcf;
    
    % --- GRÁFICOS EN TIEMPO REAL ---
    figure(1); % Selecciona la ventana de gráfico 1.
    
    % Gráfica superior: Dominio del Tiempo (Pulso)
    subplot(3, 1, [1]);
    plot(tm, It0, 'r', 'LineWidth', 2);  % Dibuja el ruido inicial en rojo.
    hold on;
    plot(tm, Iout, 'b', 'LineWidth', 1); % Dibuja el pulso actual en azul.
    xlabel('Time (ps)', 'FontSize',14); 
    xlim([-120 120]);                      % Fija el eje X entre -60 y 60 ps.
    set(gca,'XTick',-120:10:120);
    ylim([0 10]);                        % Fija el eje Y (potencia) entre 0 y 10.
    set(gca,'YTick',0:2:10);
    %legend('input', 'out');
    ylabel ('Power (a.u)','FontSize',12);
    grid on;                             % Pone cuadrícula de fondo.
    set(gca,'FontSize',12);
    hold off;

    % Gráfica central: Absorción en el tiempo
    subplot(3, 1, [2]);                    
    plot(tm,alfa0*100,'r','LineWidth',2); % Absorción inicial (rojo).
    hold on;
    plot(tm,alfa*100,'b','LineWidth',1);  % Absorción actual (azul).
    xlabel ('Time (ps)','FontSize', 14);
    xlim([-120 120]);
    set(gca,'XTick',-120:10:120);
    ylim([54 62]);                        % Zoom en el eje Y para ver detalles de absorción.
    % present time signal
    set(gca,'YTick',54:2:62);
    %legend('\alpha(t)');
    ylabel ('Absorption (%)','FontSize',12);
    grid on;
    set(gca,'FontSize',12);
    hold off;

    % Gráfica inferior: Espectro (Longitud de onda)
    subplot(3, 1, [3]);
    plot(wl,10*log10(If0),'LineWidth',2,'Color','r'); % Espectro inicial (rojo).
    hold off;
    plot(wl,Iouf,'b','LineWidth',1);      % Espectro actual (azul).
    text(1510,10, num2str(indc),'FontSize',20); % Escribe el número de vuelta en la pantalla.
    xlabel ('Wavelength (nm)');
    ylabel ('Power (dB)');
    xlim([1500 1600]);                    % Muestra solo de 1500 a 1600 nm.
    set(gca,'XTick',1500:10:1600);
    ylim([-80 100]);                      % Rango dinámico en dB.
    set(gca,'YTick',-80:40:100);
    grid on;
    set(gca,'FontSize',12);
    hold off;
    
    drawnow; % Fuerza a MATLAB a pintar la gráfica en este instante.

    % --- GUARDADO DE IMÁGENES ---
    % Comprueba si el número de vuelta es par (2, 4, 6...).
    if mod(indc,2)==0
        indx = indx +1;    % Aumenta el contador de fotos guardadas.
        cd outputring % Entra en la carpeta de resultados.
        saveas(gcf, strcat('file', num2str(indx), '.png')); % Guarda lo que se ve en pantalla como PNG.
        cd ../               % Vuelve a salir a la carpeta principal (CRUCIAL).
    end
end