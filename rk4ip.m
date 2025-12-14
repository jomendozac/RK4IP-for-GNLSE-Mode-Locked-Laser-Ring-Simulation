function Eout = rk4ip(A0, omega, dz, beta, alpha, gamma, ts, fR, hR)
    % Function rk4ip: Método Runge-Kutta de 4to Orden en el Marco de Interacción (Interaction Picture)
    % Este es el "corazón" numérico que resuelve la Ecuación No Lineal de Schrödinger (GNLSE).
    
    % A0: Campo inicial.
    % omega: Frecuencia.
    % dz: Tamaño del paso (distancia a avanzar).
    % beta, alpha: Parámetros lineales (dispersión y pérdidas).
    % gamma, ts, fR, hR: Parámetros no lineales.

    Efield = A0; % Carga el campo eléctrico actual.
    
    % --- CONCEPTO DEL MARCO DE INTERACCIÓN ---
    % Imaginamos que nos movemos con el pulso a la velocidad de la luz y
    % resolvemos la parte lineal (dispersión/pérdida) de forma exacta.
    % Las variables 'k' calculan las correcciones no lineales (distorsiones por intensidad).
    
    % --- Paso 1: Inicio del intervalo (z) ---
    
    % Aip: Campo auxiliar movido "medio paso" (dz/2) SOLO con efectos lineales.
    % Esto nos sitúa en el centro del intervalo para el marco de referencia.
    Aip = D_op(Efield, dz/2, beta, alpha, omega);
    
    % k1: Calcula la no linealidad al principio (Efield), y luego propaga esa
    % distorsión medio paso hacia adelante linealmente.
    % Es la "pendiente" inicial.
    k1 = D_op(N_op(Efield, dz, gamma, ts, fR, hR), dz/2, beta, alpha, omega);
    
    % --- Paso 2: Punto medio (z + dz/2) - Primera estimación ---
    
    % k2: Estima la no linealidad en el punto medio.
    % Usa 'Aip + k1/2' como una predicción de cómo estará el campo en la mitad.
    % Nota: Aquí no llamamos a D_op porque ya estamos en el marco del punto medio.
    k2 = N_op(Aip + k1/2, dz, gamma, ts, fR, hR);
    
    % --- Paso 3: Punto medio (z + dz/2) - Segunda estimación ---
    
    % k3: Refina la estimación del punto medio usando k2.
    % Esto corrige el error de la predicción anterior.
    k3 = N_op(Aip + k2/2, dz, gamma, ts, fR, hR);
    
    % --- Paso 4: Final del intervalo (z + dz) ---
    
    % k4: Estima la no linealidad al final del paso.
    % Primero tomamos la estimación actual en el medio (Aip + k3) y la avanzamos
    % linealmente el otro medio paso (dz/2) para llegar al final.
    % Luego calculamos N_op allí.
    k4 = N_op(D_op(Aip + k3, dz/2, beta, alpha, omega), dz, gamma, ts, fR, hR);
    
    % --- COMBINACIÓN FINAL (Promedio Ponderado) ---
    % La fórmula mágica de Runge-Kutta: k1 + 2*k2 + 2*k3 + k4
    % Se suman las correcciones no lineales con esos pesos específicos (Simpson's rule).
    % Luego, todo el conjunto se propaga linealmente el medio paso restante (D_op)
    % para completar la distancia total dz.
    
    Eout = D_op(Aip + k1/6 + k2/3 + k3/3, dz/2, beta, alpha, omega) + k4/6;
end