function alfa = SA(Ei, a0, ac, Ts, Es)
    % Function SA: Modelo de Absorbedor Saturable (Respuesta Lenta / Integradora)
    % Simula un dispositivo (como nanotubos de carbono o SESAM) cuya absorción cambia con el tiempo.
    
    % --- ENTRADAS ---
    % Ei: Campo eléctrico de entrada (el pulso que llega al absorbedor).
    % a0: Profundidad de modulación (Delta Alpha). Es la cantidad de absorción que se puede "borrar".
    % ac: Pérdidas no saturables (Alpha continuo). Absorción que siempre está ahí (impurezas, etc.).
    % Ts: Tiempo de recuperación. Cuánto tarda el material en volver a ser opaco después del pulso.
    % Es: Energía de saturación. Cuánta energía se necesita para volver transparente el material.
    
    global N t % Importa el número de puntos y el vector de tiempo globales.
    
    % Calcula el paso temporal 'dt' para las integrales.
    deltat = abs(t(2)-t(1));
    
    % Inicializa vectores con ceros para velocidad.
    Et = zeros(size(t));    % Para guardar el campo.
    alfa1 = zeros(size(t)); % Para guardar el perfil de absorción resultante.
    intga = zeros(size(t)); % Variable auxiliar para la integral.
    
    % Acumuladores para resolver la ecuación diferencial.
    intga0 = 0;
    intga1 = 0;
    
    % --- BUCLE TEMPORAL (Dinámica del Absorbedor) ---
    % A diferencia de la fibra (que se resuelve en frecuencia con FFT), 
    % el SA tiene "memoria", por lo que se resuelve punto a punto en el tiempo.
    for k = 1:N
        Et(k) = Ei(k); % Toma el valor del campo en el instante k.
        

        % --- CÁLCULO DE LA TASA DE RELAJACIÓN Y SATURACIÓN ---
        %DIRECTAMENTE DE LA ECUACION DIFERENCIAL DEL PAPER, despejando
        %dalpha

        % intga(k) representa la velocidad a la que cambia el estado del absorbedor.
        % 1/Ts: Tasa de recuperación natural (los electrones vuelven a su sitio).
        % abs(Et)^2/Es: Tasa de saturación forzada por la luz intensa.
        intga(k) = deltat * (1/Ts + abs(Et(k))^2 / Es); 
        
        % --- SOLUCIÓN DE LA ECUACIÓN DIFERENCIAL (Método de Integración) ---
        % El código resuelve la ecuación: dq/dt = (q0 - q)/Ts - q*P(t)/Es
        % Donde 'q' es la absorción saturable.
        
        % Acumulador 1 (Factor Integrante): Suma la historia de relajación
        % y saturación del alpha en cada punto de la simulacion de 1 a N
        intga0 = intga0 + intga(k);
        
        % Acumulador 2 (Término fuente):
        % Integra la recuperación natural ponderada por el estado acumulado.
        intga1 = intga1 + a0/Ts*deltat*exp(+intga0);  
        
        % --- CÁLCULO FINAL DE LA ABSORCIÓN EN EL INSTANTE k ---
        % Fórmula analítica de la solución de la ecuación diferencial.
        % ac: Pérdida base fija.
        % El resto es la parte variable que baja cuando pasa el pulso y sube después.
        alfa1(k) = ac + exp(-intga0)*(intga1 + a0);
    end
    
    % Devuelve el perfil de absorción temporal completo.
    alfa = alfa1;
end