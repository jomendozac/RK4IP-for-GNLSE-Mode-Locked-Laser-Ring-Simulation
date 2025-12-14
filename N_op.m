function Aout = N_op(Ain, dz, gamma, t_shock, fR, hR_w)
    % Function N_op: Operador No Lineal N para la Ecuación GNLSE
    % Calcula la parte no lineal de la ecuación: hN(A)
    % Es decir, calcula cómo cambia el pulso debido a su propia intensidad.
    
    % Ain: Campo óptico de entrada (vector complejo en el tiempo).
    % dz: Paso de propagación (distancia).
    % gamma: Coeficiente de no linealidad (qué tan fuerte responde el material a la luz).
    % t_shock: Tiempo de choque óptico (1/omega0), relacionado con el auto-inclinamiento.
    % fR: Fracción Raman (qué porcentaje de la no linealidad es retardada/molecular).
    % hR_w: Respuesta Raman en el dominio de la frecuencia (calculada previamente).

    global t % Accede al vector de tiempo global.
    
    % Calcula el paso de tiempo 'dt' restando el segundo punto menos el primero.
    % Es necesario para calcular integrales y derivadas numéricas correctamente.
    dt = abs(t(2)-t(1));
    
    % --- CÁLCULO DEL EFECTO RAMAN ---
    % El efecto Raman es una respuesta "lenta" del material. Depende de la historia del pulso.
    % Matemáticamente es una convolución: Integral(hR(t') * |A(t-t')|^2)
    % que es el ultimo termino de la GNLSE
    % Por el teorema de convolución, en frecuencia es una simple multiplicación.
    
    % 1. abs(Ain).^2: Calcula la potencia instantánea del pulso.
    % 2. fft(...): Pasa la potencia al dominio de frecuencia.
    % 3. .* hR_w: Multiplica por la respuesta Raman del material.
    % 4. ifft(...): Vuelve al dominio del tiempo.
    % 5. * dt: Factor de escala necesario para la integral discreta.
    hR_A2 = ifft(hR_w .* fft(abs(Ain).^2)) * dt;
    %asi obtenemos el termino de aporte del efecto RAMAN 
    
    % --- TÉRMINO 1: AUTOMODULACIÓN DE FASE (SPM) + RAMAN ---
    % Computamos la parte de N que no tiene derivadas, teniendo ya la convolucion hecha del paso anterior.
    %Aqui se tienen en cuenta los dos efectos no lineales Kerr y raman asi
    % (1-fR) .* Ain .* abs(Ain).^2 : Efecto Kerr puro (instantáneo). La fase cambia con la intensidad.
    % fR .* Ain .* hR_A2 : Efecto Raman. La fase cambia con la "memoria" de
    % la intensidad, por eso Raman es "lento"
    N1 = 1i * gamma * ((1-fR) .* Ain .* abs(Ain).^2 + fR .* Ain .* hR_A2);
    
    % --- TÉRMINO 2: AUTO-INCLINACIÓN (SELF-STEEPENING / OPTICAL SHOCK) ---
    % Este efecto ocurre en pulsos muy cortos. La velocidad de grupo depende de la intensidad.
    % Hace que el pico del pulso viaje más lento que las colas, creando un frente de choque.
    % Matemáticamente es proporcional a la derivada temporal de la polarización no lineal.
    
    % Repite el cálculo de la polarización no lineal (lo que está dentro del paréntesis de N1):
    term_inside = ((1-fR) .* Ain .* abs(Ain).^2 + fR .* Ain .* hR_A2);
    
    % Calcula la derivada temporal usando 'gradient'.
    % gradient(vector, paso): Calcula derivada numérica central d/dt.
    % Multiplica por -gamma * t_shock (donde t_shock es aprox 1/frecuencia_central).
    N2 = -gamma * t_shock .* gradient(term_inside, dt);
    
    % --- SALIDA FINAL ---
    % Suma ambos efectos (N1 y N2) y multiplica por la distancia 'dz'.
    % Esto es equivalente al operador N total
    % Esto devuelve el cambio total que sufrirá el campo 'Ain' debido a efectos no lineales
    % al recorrer la distancia 'dz'.
    Aout = dz * (N1 + N2);
end