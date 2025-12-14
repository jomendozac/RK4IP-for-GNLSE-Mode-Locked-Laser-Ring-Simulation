function Aout = D_op(Ain, dz, beta, alpha, omega)
    % Function D_op: Operador Lineal D (Dispersión y Pérdidas) para la GNLSE
    % Calcula la parte lineal de la ecuación: exp(hD * z) * A
    
    % Ain: Campo óptico de entrada (en el dominio del tiempo).
    % dz: Paso de propagación (distancia pequeña que avanza la luz).
    % beta: Vector con los coeficientes de dispersión [beta1, beta2, beta3...].
    % alpha: Coeficiente de pérdidas lineales de la fibra.
    % omega: Vector de frecuencias angulares.

    % Inicializa el vector de dispersión Dw con ceros, del mismo tamaño que el vector de frecuencias.
    Dw = zeros(1, length(omega));
    
    % --- Bucle para sumar los términos de dispersión de orden superior ---
    % Suma de la Serie de Taylor para la dispersión: Sum(i * beta(j) * (-w)^j / j!)
    for jl = 1:length(beta)
        % En cada iteración (jl = 1, 2, 3...) se añade un término de dispersión.
        % 1i: Es la unidad imaginaria 'i' en MATLAB.
        % beta(jl): Es el coeficiente de dispersión de ese orden (ej. beta2, beta3).
        % factorial(jl): Divide por el factorial del orden (1!, 2!, etc.).
        % (-omega).^(jl): Eleva la frecuencia negativa a la potencia del orden actual.
        %OJO = omega es 2*pi*f y esta reemplazando al operador d/dt pues es
        %equivalente derivar A que multiplicar por omega.

        
        % NOTA: Aquí el código ejecuta una resta del término imaginario.
        Dw = Dw - 1i.* beta(jl)/factorial(jl) .* (-omega).^(jl);
    end 
 
    % --- Aplicación del Operador en el Dominio de la Frecuencia ---
    % El método Split-Step Fourier resuelve la parte lineal (D) en frecuencias:
    
    % 1. fft(Ain): Convierte la señal de entrada al dominio de la frecuencia.
    % 2. exp(dz * (-Dw - alpha/2)): Crea el "propagador".
    %    -Dw: Aplica la fase correspondiente a la dispersión acumulada.
    %    -alpha/2: Aplica la atenuación (pérdida de potencia) en la distancia dz.
    % 3. .* : Multiplica la señal en frecuencia por el propagador.
    % 4. ifft(...): Devuelve la señal al dominio del tiempo usando la Transformada Inversa de Fourier.
    
    Aout = ifft(exp(dz * (-Dw - alpha/2)) .* fft(Ain));
end