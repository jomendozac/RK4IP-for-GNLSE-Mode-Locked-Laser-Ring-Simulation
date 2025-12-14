function [Uout, numFFT] = Propagation_SMF(imp, tol, long, alpha, beta, gamma, ts, fR, hR)
    % Propagation_SMF: Simula la propagación en Fibra Monomodo (Pasiva).
    % Utiliza el algoritmo RK4IP con tamaño de paso adaptativo.
    
    % --- ENTRADAS ---
    % imp: Impulso/Campo óptico de entrada.
    % tol: Tolerancia de error local.
    % long: Longitud de la fibra.
    % alpha, beta, gamma: Parámetros de la fibra (Pérdidas, Dispersión, No linealidad).
    % ts, fR, hR: Parámetros de auto-inclinación y efecto Raman.

    global f t % Variables globales de frecuencia y tiempo.
    
    % Frecuencia angular.
    omega = 2*pi*f;
    
    % Carga el campo eléctrico inicial.
    E_z = imp;
    
    % Mensaje inicial en consola.
    fprintf(1, '\nCaculation in process...')
    
    % --- INICIALIZACIÓN ---
    % Define un paso inicial muy pequeño para arrancar con precisión.
    dz = long / 10000; 
    
    % Guarda los valores base de los parámetros físicos.
    alpha0 = alpha;
    beta0 = beta;
    gamma0 = gamma;
    
    Z_prop = 0; % Posición actual en la fibra (metros).
    ii = 0;     % Contador de iteraciones.
    
    % --- BUCLE PRINCIPAL DE PROPAGACIÓN ---
    while Z_prop < long
        ii = ii + 1;
        
        % Ajuste de borde: Si el paso 'dz' nos hace salir de la fibra, lo recortamos.
        if Z_prop + dz > long
            dz = long - Z_prop;
        end
        
        % Avanza la posición acumulada.
        Z_prop = Z_prop + dz;
        
        % --- VARIACIÓN LONGITUDINAL DE PARÁMETROS ---
        % Este bloque permite simular fibras que cambian sus propiedades a lo largo (ej. fibras conicas).
        % Como m=0, en esta simulación los parámetros se mantienen constantes.
        m = 0;
        alpha = alpha0*(1 + m*Z_prop);
        beta = beta0*(1 + m*Z_prop);
        gamma = gamma0*(1 + m*Z_prop);
        
        % --- MÉTODO RK4IP (Runge-Kutta 4th Order Interaction Picture) ---
        % Calculamos la evolución del pulso de dos maneras para comparar errores.
        
        % 1. Paso Fino (Uf): Avanzamos la distancia 'dz' dando DOS medios pasos (dz/2 + dz/2).
        % Esto es computacionalmente más costoso pero más preciso.
        Uf = rk4ip(rk4ip(E_z, omega, dz/2, beta, alpha, gamma, ts, fR, hR), omega, dz/2, beta, alpha, gamma, ts, fR, hR);
        
        % 2. Paso Tosco (Uc): Avanzamos la distancia 'dz' en UN solo paso completo.
        Uc = rk4ip(E_z, omega, dz, beta, alpha, gamma, ts, fR, hR);
        
        % --- CÁLCULO DE ERROR Y ADAPTACIÓN ---
        
        % Error relativo local: Diferencia entre el cálculo fino y el tosco.
        error = sqrt(sum(abs(Uf-Uc).^2)) / sqrt(sum(abs(Uf).^2));
        
        % Factor de corrección: Cuánto debemos cambiar el paso 'dz' para mantenernos en la tolerancia.
        factor = tol/error;
        
        % Decisión: ¿Aceptamos el paso?
        if error > 2*tol
            % RECHAZADO: El error es muy alto.
            % Deshacemos el avance (restamos dz a Z_prop) y repetiremos el bucle con un paso más pequeño.
            Z_prop = Z_prop - dz;
        else
            % ACEPTADO: El error está dentro de los límites.
            % Extrapolación de Richardson: Combinamos Uf y Uc para obtener una solución de 5to orden (muy precisa).
            E_z = 16/15 * Uf - 1/15 * Uc;
        end
        
        % Actualizamos el tamaño del paso 'dz' para la SIGUIENTE iteración.
        % Se usa la raíz quinta (^(1/5)) porque el método base es de orden 4.
        dz = dz * factor^(1/5);
        
        % Barra de progreso en consola:
        % '\b' es backspace (borrar carácter). Borra el porcentaje anterior y escribe el nuevo.
        fprintf(1, '\b\b\b\b\b\b%5.2f%%', Z_prop* 100.0 /long);
    end
    
    % Salida final: Campo eléctrico tras recorrer toda la fibra.
    Uout = E_z;
    
    % Estimación del coste computacional (número de FFTs realizadas).
    % 3 FFTs por paso de RK4IP * 16 pasos internos * ii iteraciones (Aproximado).
    numFFT = 3*16*ii; 
end