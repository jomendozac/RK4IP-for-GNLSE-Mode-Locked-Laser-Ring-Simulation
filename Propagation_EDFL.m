function Uout = Propagation_EDFL(Uin, tol, L, g0, Es, alpha, beta, gam, ts, fR, hR)
    % Function Propagation_EDFL: Propagación en Fibra Dopada con Erbio (Activa)
    % Simula la amplificación de la señal considerando saturación y ancho de banda finito.
    
    global f t % Importa vectores globales de frecuencia y tiempo.
    
    % Calcula el paso temporal 'dt' para realizar integrales de energía.
    delta_t = abs(t(2)-t(1));
    
    % Vector de frecuencias angulares.
    omega = 2*pi*f;
    
    % --- CONSTANTES LOCALES ---
    c0 = 3e8;               % Velocidad de la luz (m/s).
    lambdaP = 1550e-9;      % Longitud de onda de bombeo/central.
    wP = 2*pi*c0/lambdaP;   % Frecuencia angular central.
    
    % Ancho de banda de ganancia (Gain Bandwidth).
    delta_wl = 20e-9;       % 20 nm de ancho de banda.
    
    % Ancho de banda en frecuencia angular (OmegaP).
    % Convierte el ancho de banda de nm a rad/s.
    OmegaP = (2*pi*c0/lambdaP^2)*delta_wl;
    
    % Inicializa el campo eléctrico.
    Efield = Uin;
    
    % Mensaje de estado.
    fprintf(1, '\nCaculation in process...');
 
    % Paso inicial de propagación (pequeño para precisión al arranque).
    dz = L / 10000;
    
    % Guardar parámetros físicos iniciales.
    alpha0 = alpha;
    beta0 = beta;
    gamma0 = gam;
    
    % Inicializar contadores de distancia y pasos.
    Z_prop = 0;
    ii = 0;
   
    % --- BUCLE DE PROPAGACIÓN ---
    while Z_prop < L
        ii = ii+1;
        
        % Ajuste para no exceder la longitud total L en el último paso.
        if Z_prop + dz > L
            dz = L - Z_prop;
        end
        Z_prop = Z_prop + dz;
        
        % --- CÁLCULO DE PARÁMETROS DINÁMICOS DEL AMPLIFICADOR ---
        
        % Constante de tiempo asociada al ancho de banda de ganancia (Parámetro T2).
        tt= 1/OmegaP;
        
        % Cálculo de la Ganancia Saturada (g).
        % La ganancia baja si la energía del pulso (integral de |E|^2) es alta comparada con Es.
        g = g0./(1 + sum(abs(Efield).^2)*delta_t/Es);
        
        % Actualización de Alpha (Pérdida Neta).
        % En la ecuación, alpha positivo es pérdida. Restamos la ganancia 'g'.
        alpha = alpha0 - g;
        
        % Dispersión de Ganancia (Efecto Filtro Espectral).
        % Se añade un término imaginario a la dispersión beta2.
        % Físicamente, esto actúa como un filtro gaussiano que limita el ancho de banda del pulso.
        beta(2) = beta0(2) + 1i*tt*tt*g;
        
        % Perfil de no linealidad (constante en este caso, m=0).
        m = 0;                              
        gam = gamma0*(1 + m*Z_prop);
        
        % --- MÉTODO RK4IP ADAPTATIVO ---
        
        % 1. Paso Fino (Uf): Avanza 'dz' dando dos medios pasos (dz/2).
        Uf = rk4ip(rk4ip(Efield, omega, dz/2, beta, alpha, gam, ts, fR, hR), omega, dz/2, beta, alpha, gam, ts, fR, hR);
        
        % 2. Paso Tosco (Uc): Avanza 'dz' dando un solo paso completo.
        Uc = rk4ip(Efield, omega, dz, beta, alpha, gam, ts, fR, hR);
        
        % Cálculo del error relativo local entre ambos métodos.
        error = sqrt(sum(abs(Uf-Uc).^2)) / sqrt(sum(abs(Uf).^2));
        
        % Factor de seguridad para ajustar el siguiente paso 'dz'.
        factor = tol/error;
        
        % --- CONTROL DE PASO ---
        if error > 2*tol
            % Si el error es inaceptable, descartamos este avance.
            % Retrocedemos Z_prop y repetimos el bucle (dz será menor en la sig. vuelta).
            Z_prop = Z_prop - dz;
        else
            % Si el error es aceptable, actualizamos el campo eléctrico.
            % Usamos extrapolación de Richardson (combinación lineal) para mayor precisión.
            Efield = 16/15 * Uf - 1/15 * Uc;
        end
        
        % Calculamos el nuevo tamaño de paso 'dz' para la siguiente iteración.
        % Se usa potencia 1/5 porque RK4 es un método de orden 4.
        dz = dz * factor^(1/5);
    end
    
    % Devuelve el campo al final de la fibra.
    Uout = Efield;
end