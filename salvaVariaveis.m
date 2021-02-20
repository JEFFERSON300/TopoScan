delete 'Abertura_Distancia.xlsx';
filename = 'Abertura_Distancia.xlsx';
xlswrite(filename,Abertura_Distancia',1,'A1:A230') % variavel que esta dentro de salva.txt dividido por 2

NumeroFraturas = 230; % precisa definir, variavel que esta dentro de salva.txt
tempAber_Dist= Abertura_Distancia(1:NumeroFraturas);

acrescenta = 0;
vetor = zeros(1,NumeroFraturas)';
Ab = zeros(1,NumeroFraturas)';

for jk=2:2:NumeroFraturas     % variavel que esta dentro de salva.txt
    acrescenta = acrescenta + 1;
    vetor(jk-1) = acrescenta;
    vetor(jk) = acrescenta;
    Ab(jk-1) = tempAber_Dist(acrescenta);
    Ab(jk) = tempAber_Dist(acrescenta);
end

VariavelSalva = [vetor coordinates(1:NumeroFraturas,:) Ab];

delete 'coordenadas.xlsx';
filename = 'coordenadas.xlsx';
xlswrite(filename,VariavelSalva,1,'A1:D230') % variavel que esta dentro de salva.txt

delete 'familias.xlsx';
filename = 'familias.xlsx';
xlswrite(filename,familias,1,'A1:B64') 

fileID = fopen('fatorReal.txt','w');
fprintf(fileID,'%f\n',FatorReal);
fclose(fileID);