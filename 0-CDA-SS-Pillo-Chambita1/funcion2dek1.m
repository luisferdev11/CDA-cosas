%%Jorge Caballero Rosales. Ing. Robótica Industrial. 20/03/2022

%function=SolveU_Iterative_1(K1)
clc
close all
S = mfilename('fullpath'); f = filesep; ind = strfind(S, f); S1 = S(1:ind(end) - 1); cd(S1);

%Primera K1
load('Ws_60x60x20_3D_SoftK_S1SM1FDD0D2Sy0R0_GSBESO3D_V4_SoftK_se.mat')

%Mayor dificultad
%load('Ws_120x40x30__10_3_0.6_0.8_3_3_0.5_17_35_16_1HardK_S1SM2FDD2D2Sy0V0_0.03_0.15_GSBESO3D_V4_HardK_1.mat')

K1(freedofs, freedofs);
it1 = 0;
it2 = 0;
it3 = 0;
rr1 = 0;
rr2 = 0;
rr3 = 0;
X1 = 0;
X2 = 0;
X3 = 0;
%Graficamos K1 y damos pausa para que se pueda analizar la grafica.
spy(K1);
saveas(figure(1), 'Grafica de K1.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%Equilibrar y reordenar.%%%%%%%%%%%%
n = input('\n 1:Equilibrar y reordenar K1 \n 2:No realizarle ninguna operacion a K1 \n');

switch n
    case 1
        [P, R, C] = equilibrate (K1); %equilibramos
        K1new = R * P * K1 * C;
        Fnew = R * P * F(freedofs);
        %condest (K1new)
        q = dissect(K1new); %Reordenamos
        K1 = K1new(q, q);
        F = Fnew(q);

        spy(K1new, 'ro')
    otherwise
        disp('No se equilibro ni reordeno la matriz K1')
end

fprintf('\n')
disp('Presione cualquier tecla para continuar');

pause;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Seleccionamos un preaciondicionador%%%%

% Y=input('1: Seleccionar un preacondicionador. \n 2:Terminar')
%
for i = 1:100;
    x = input(' \n1:Realizar PCG con distintos métodos \n2: Tabla de resultados \n3:Terminar.\n')

    switch x
            %PCG con factorizacion incompleta LU
        case 1

            SolveU1 = tic;
            size(K1)
            [L, I] = ilu(K1(freedofs, freedofs));
            [x1, fl1, rr1, it1] = pcg(K1(freedofs, freedofs), F(freedofs), 1e-8, 3000, L, I);
            fprintf('%e - resiuduo relativo \n%d - numero de iteraciones \n', rr1, it1);

            fprintf('Se calculo el método pcg con factorizacion incompleta LU \n')
            Tiempodesolucion1 = toc(SolveU1);
            X1 = Tiempodesolucion1
            disp('Presione cualquier tecla para continuar');

            pause;
            
            
            %PCG con factorizacion incompleta ichol.
            SolveU2 = tic;
            size(K1)
            
            alpha = .01;
            L = ichol(K1(freedofs, freedofs), struct('type','ict','droptol',1e-1,'diagcomp', alpha));
            [x2, fl2, rr2, it2] = pcg(K1(freedofs, freedofs), F(freedofs), 1e-8, 3000, L, L');
            
            
            fprintf('%e - resiuduo relativo \n%d - numero de iteraciones \n', rr2, it2);
            fprintf('Se calculo el método pcg con factorizacion incompleta ichol \n')
            Tiempodesolucion2 = toc(SolveU2);
            X2 = Tiempodesolucion2
            disp('Presione cualquier tecla para continuar');
            pause;
            
            
            %PCG con matriz identidad.
            SolveU3 = tic;
            size(K1)
            L = speye(max(size(freedofs)));
            [x3, fl3, rr3, it3] = pcg(K1(freedofs, freedofs), F(freedofs), 1e-8, 3000, L, L');
            fprintf('%e - resiuduo relativo \n%d - numero de iteraciones \n', rr3, it3);
            fprintf('Se realizó con exito el PCG con matriz identidad \n')
            Tiempodesolucion3 = toc(SolveU3);

            X3 = Tiempodesolucion3
            disp('Presione cualquier tecla para continuar');
            pause;

        case 2

            fprintf('Metodo \t\t\t No de iteracion \t\t Residuo relativo \t\t Tiempo de solucion');
            fprintf('\nPCG con ILU \t\t\t %d \t\t\t %e \t\t\t %d', it1, rr1, X1);
            fprintf('\nPCG con ICHOL \t\t\t %d \t\t\t %e \t\t\t %d', it2, rr2, X2);
            fprintf('\nPCG con Identidad \t\t %d \t\t\t %e \t\t\t %d \n', it3, rr3, X3);

        case 3
            break

    end

    pause;
end



