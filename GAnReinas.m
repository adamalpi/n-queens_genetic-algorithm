function [ mejorTablero ] = GAnReinas( n ,nm, maxi)
%GAn Resulve el problema de las n-reinas con un genético
% n : número de reinas
% nm : número de muestras de cada generación
% maxi : número máximo de generaciones
% mejorValor: vector con el mejor Valor de cada generación

valorObjetivo = sum(1:n-1);
fprintf('Ataques posibles: %d\n',valorObjetivo);

if (mod(nm,2)==1)
    nm=nm+1;
end
% poblacion inicial
poblacion=NaN(nm,n);
for i = 1:nm
    poblacion(i,:) = randperm(n);
end
ObjV = EvaluaPoblacion(poblacion);



mejorValor = zeros(1,1);
mejorIndividuo = zeros(maxi,n);
gen = 0;                     
while ((gen < maxi) && (max(mejorValor) ~= valorObjetivo)),
    
    poblacion = siguienteGeneracion(poblacion,ObjV);
    
    ObjV = EvaluaPoblacion(poblacion);
    
    %Almacenar media del Valor y mejor individuo
    mejor = find(ObjV==max(ObjV));
    mejorValor(gen+1) = ObjV(mejor(1));
    mejorIndividuo(gen+1,:) = poblacion(mejor(1),:);
    
    %incrementar el numero de generaciones
    gen = gen + 1;

end
[mejorVal, mejorIndPos]= max(mejorValor);
mejorTablero=mejorIndividuo(mejorIndPos,:);
fprintf('Ataques en la solucion: %d\n', ...
    valorObjetivo-mejorVal);

%pintamos el trablero
for fila = 1:n
    for col = 1:n
        if(mejorTablero(col)==fila)
            fprintf('r ');
        else
            fprintf('- ');
        end
    end
    fprintf('\n');
end
return;

function SG = siguienteGeneracion(AG, ObjV)
%    Crea una nueva generacion
%    AP: matriz con la población actual
%    ObjV: vector con la valoración de la población actual
%    porcNewInd: porcentaje de nuevos individuos en la nueva generación. 
%    funcSelec: nombre de la función de selección
%    funcCruce: nombre de la función de cruce
%    funcMutacion: nombre de la función de mutacion
%    porcMutacion: porcetaje de individuos a los que se les aplica la
%                  mutación

[nM, nR] = size(AG);
SG = zeros(nM,nR);
pc=1/2;
pm=0.1;
n = 0; %nuevos individuos
parIndex=[NaN NaN];
parMuestras = NaN(2,nR);
aux = NaN(1,nR);
suma=sum(ObjV);
while n < nM
    %REPRODUCCION
    % seleccionamos la primera muestra
    aleatorio = rand() * suma;
    parIndex(1) = 1;
    while (parIndex(1) <= length(ObjV) && sum(ObjV(1:parIndex(1))) < aleatorio)
        parIndex(1) = parIndex(1) + 1;
    end
    parIndex(2) = parIndex(1);
    % seleccionamos la segunda muestra distinta de la primera
    while parIndex(2) == parIndex(1)
        parIndex(2) = 1;
        aleatorio = rand() * suma;
        while (parIndex(2) <= length(ObjV) && sum(ObjV(1:parIndex(2))) < aleatorio)
           parIndex(2) = parIndex(2) + 1;
        end
    end
    parMuestras(1,:)=AG(parIndex(1),:);
    parMuestras(2,:)=AG(parIndex(2),:);
    %CRUCE
    if(rand<pc)
        % se selecciona el punto de cruce al azar
        pto_cruce = round((nR-2) * rand())+1;
        %primero la muestra 1
        aux1=NaN;
        aux1(1:pto_cruce) = parMuestras(1,1:pto_cruce);
        for i = parMuestras(2,:)
            if(~any(aux1==i))
                aux1 = [aux1 i];
            end
        end
        %despues la muestra 2
        aux2=NaN;
        aux2(1:pto_cruce) = parMuestras(2,1:pto_cruce);
        for i = parMuestras(1,:)
            if(~any(aux2==i))
                aux2 = [aux2 i];
            end
        end
        parMuestras(1,:)=aux1;
        parMuestras(2,:)=aux2;
    end
    
    %MUTACION
    if (rand<pm)
        aux=parMuestras(1,:);
        c = [1 1];
        while c(1) == c(2) 
            c = round(rand(2,1) * (nR-1))+1;
        end
        parMuestras(1,c(1)) = aux(c(2));
        parMuestras(1,c(2)) = aux(c(1));
    end
    
    if (rand<pm)
        aux=parMuestras(2,:);
        c = [1 1];
        while c(1) == c(2) 
            c = round(rand(2,1) * (nR-1))+1;
        end
        parMuestras(2,c(1)) = aux(c(2));
        parMuestras(2,c(2)) = aux(c(1));
    end
    % se inserta el primer hijo en la poblacion
    SG(n+1:n+2,:)=parMuestras;
    n=n+2;
end

return;

function [ ObjV ] = EvaluaPoblacion( poblacion)
% Calcula el valor de los individuos de la población
%    poblacion: matriz con los individuos de la población
tam = size(poblacion);
ObjV=zeros(tam(1),1);
for i = 1:tam(1)
    tablero=poblacion(i,:);
    ataques = 0;
    for j = 1:(tam(2)-1)
        for k = (j+1):tam(2)
            if (k-j) == abs(tablero(j)-tablero(k)) % misma diagonal 
                ataques = ataques + 1;
            end
        end
    end
    %se maximiza el problema suponiendo inc como el numero de ataques
    %maximo, el cual se resta al numero de ataques actual.
    inc = sum(1:tam(2)-1);
    ObjV(i) = -1 * ataques + inc;
    %el inc sale de que en el peor de los casos hay inc ataques
end
return
