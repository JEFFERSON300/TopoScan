%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Autor: Jefferson P. S. SIlva
% Versao: TopoScanv1.0
% 
% Rotina para caracterização geologica:analise topologica e scanline   
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc

%%----------------------------------------------------------------------------------
% variaveis a serem definidas

%AngFamilias = 45;       % angulo que delimita cada familia (número par ou multiplo de 9)
n = 1000;                % número maximo de fraturas
fator = 60;              % fator de multiplicação para abertura da imagem utilizada para a scanline

%%----------------------------------------------------------------------------------
% carrega imagem e imprime na tela

message = sprintf('A para todas as etapas \nB apenas scanlines e analise topologica');
uiwait(helpdlg(message));
opcao = GUI_41();
i = 0;

if opcao == 'A'

    filename = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
              '*.*','All Files' },'mytitle',pwd);
    I = imread(filename); % imagem que será carregada
    I2 = rgb2gray(I);
        
    K = medfilt2(I2);
    J = imgaussfilt(I2);
    L = stdfilt(I2);
            
    figure;
    subplot(1,3,1);
    imshow(K); 
    title('1 - Filtro de Média'); 
            
    subplot(1,3,2);
    imshow(J); 
    title('2 - Filtro Gaussiano');
            
    subplot(1,3,3);
    imshow(L); 
    title('3 - Filtro Desvio Padrão');
    
    ordem = str2num(GUI_43());
    pause(15);
            
    if ordem==1
       imwrite(K,'filtro.png');         
    end
            
    if ordem==2
       imwrite(J,'filtro.png')           
    end
            
    if ordem==3
       imwrite(K,'filtro.png')           
    end
    
    close all
    I4 = imread('filtro.png');
    
    pout_imadjust = imadjust(I4);
    pout_histeq = histeq(I4);
    pout_adapthisteq = adapthisteq(I4);

    figure ;
    subplot(1,3,1);
    imshow(pout_imadjust); 
    title('1 - Aumento de Constraste');

    subplot(1,3,2);
    imshow(pout_histeq); 
    title('2 - Histograma');
            
    subplot(1,3,3);
    imshow(pout_adapthisteq); 
    title('3 - Histograma Adaptativo');

    ordem = str2num(GUI_43());
    pause(15);

    if ordem==1
        imwrite(pout_imadjust,'contraste.png')
    end

    if ordem==2
        imwrite(pout_histeq,'contraste.png')
    end

    if ordem==3
        imwrite(pout_adapthisteq,'contraste.png')
    end
    
    close all
    
    I5 = imread('contraste.png');
    figure, imshow(I5)
    
    %%----------------------------------------------------------------------------------
    % marcar os pontos do objeto de referencia e receber o tamanho na escala
    % real

    dimensaoObjeto = 0;
    i = 0;
    coordinatesObjeto = zeros(2,2);
    button0 = 0;
    
    AngFamilias = str2double(GUI_42());
    
    message = sprintf('Marque os pontos de incio e fim do objeto de referencia com o mouse (Linhas)\n D para apagar linha \n');
    uiwait(helpdlg(message));
    
    while (dimensaoObjeto==0)
    
        [a,b,button0] = ginput(1);
        if button0~=100
            i = i + 1;
            coordinatesObjeto(i,:) = [a, b];
        end
    
        if (mod(i,2)==0)
            dimensaoObjeto = GUI_40();
            dimensaoObjeto = str2double(dimensaoObjeto);
            H = line([coordinatesObjeto(i-1,1) coordinatesObjeto(i,1)],[coordinatesObjeto(i-1,2) coordinatesObjeto(i,2)],'Color','red','LineWidth',1);
            [a,b,button0] = ginput(1);
            
            if button0 == 100
                delete(H);
                coordinatesObjeto(i,:) = [];   %
                coordinatesObjeto(i-1,:) = []; % apaga linha e os pontos pertencentes a linha
                i = i - 2;                     % 
                dimensaoObjeto = 0;
                continue
            end
        
            if button0 ~= 100 && i==2
                X =  [coordinatesObjeto(i-1,1),coordinatesObjeto(i-1,2);coordinatesObjeto(i,1),coordinatesObjeto(i,2)];
                distanciaObjeto = pdist(X,'euclidean');
                continue
            end      
        end   
    end

    FatorReal = dimensaoObjeto/distanciaObjeto;  % densidade do objeto de referencia cm/pixels

    %%----------------------------------------------------------------------------------
    % marcar os pontos das fraturas 

    message = sprintf('Marque os pontos de incio e fim das Fraturas com o mouse (Linhas)\n D para apagar a fratura anterior\n E para terminar');
    uiwait(helpdlg(message));

    ContfamiliasVet = ones(1,90/AngFamilias);  %%%%%%%%%%%
    T = size(I);
    blankimage = zeros(T(1),T(2),3);
    cont = 0;
    coordinates = zeros(n,2);
    coeficientes = zeros(n,2);
    Abertura_Distancia = zeros(n,2);
    familias = zeros(n,90/AngFamilias);  %%%%%%%%%%%%%%
    i = 0;
    button = 0;
    soma = 1;
    troca = 0;
    ij = 0;
    abertura = 0;
    contAbertura = 1;
    contDistancia = 1;
    distance = 0;
    pointsDistance = zeros(n,1);
    Area_Fracture = 0;
    Area_Figure = 0;
    
    hold on

    while (button~=101)

        bb = 0;                           % variavel que irá receber as letras para controle da marcação das aberturas
        coordinatesAB = zeros(2,2);       % variavel que ira receber os pontos das aberturas
        [x,y,button] = ginput(1);         % recebe pos pontos ou as letras definidas pelo usuario para marcação das fraturas
        if (button~=100 && button~=101)   % 100 representa a letra D (deletar linha) maiuscula e 101 a letra E (end da inserção de fraturas) maiuscula
            i = i + 1;
            coordinates(i,:) = [x, y];    % atribui os pontos recebidos pelo usuario para a variavel cordendadas
        end
        if button == 101                  % se o botão 101 letra E maiuscula for digitada pulo o loop para a proxima rodada
           continue                       %
        end                               %

        if mod(i,2)==0 && button~=101    % se dois pares de pontos já foram definidos se tem uma reta

                %------------ marcação das aberturas -------------
                message = sprintf('Marque dois pontos para definir a abertura\n x para confirmar \n');
                uiwait(helpdlg(message));

                while (bb~=120)             % 120 representa o botão x minusculo
                    [A,B,bb] = ginput(1);   % recebe os pontos ou as letras digitadas
                    if bb==120              
                        continue
                    end
                    ij = ij + 1;
                    if ij>2
                        coordinatesAB(1,:) = coordinatesAB(2,:);    % salva apenas as duas ultimas coordenadas da abertura
                        coordinatesAB(2,:) = [A, B];
                    else
                        coordinatesAB(ij,:) = [A, B];               % define coordenadas da abertura
                    end
                end
                ij = 0;
                XXX = [coordinatesAB(2,1),coordinatesAB(2,2);coordinatesAB(1,1),coordinatesAB(1,2)];    % calcula distancia entre os dois pontos
                abertura = pdist(XXX,'euclidean');                                                      % 
                %--------------------------------------------------

                H = line([coordinates(i-1,1) coordinates(i,1)],[coordinates(i-1,2) coordinates(i,2)],... % imrpime a fratura na imagem
                    'Color','red','LineWidth',11*abertura*FatorReal);                                    % 
                [x,y,button] = ginput(1);    % solicita pontos ou algo a ser feito

                if button == 100             % 100 representa a letra D maiuscula (deletar linha)
                    delete(H);
                    coordinates(i,:) = [];   %
                    coordinates(i-1,:) = []; % apaga linha e os pontos pertencentes a linha
                    coordinatesAB(:) = [];
                    i = i - 2;               % 
                    continue
                end

%                 if button == 101             % 101 representa a letra E maiuscula (terminar a inseção de fraturas)
%                    continue; 
%                 end

                if button ~= 100 && mod(i,2)==0
                    if cont == 0
                        
                        RGB1 = insertShape(blankimage,'Line',[coordinates(i-1,1), coordinates(i-1,2), ...
                            coordinates(i,1), coordinates(i,2)],'LineWidth',3,'Color', {'white'});     % desenha linhas abertua 3 pixel
                        RGB2 = insertShape(blankimage,'Line',[coordinates(i-1,1), coordinates(i-1,2), ...
                            coordinates(i,1), coordinates(i,2)],'LineWidth',1,'Color', {'white'});     % desenha linhas abertura 1 pixel
                        RGB3 = insertShape(blankimage,'Line',[coordinates(i-1,1), coordinates(i-1,2), ...
                            coordinates(i,1), coordinates(i,2)],'LineWidth',round(fator*abertura*FatorReal),'Color', {'white'});    % desenha linha com largura de pixel de acordo com a abertura
                        imwrite(RGB1,'endPoints.jpg');                  % 
                        imwrite(RGB2,'branchPoints.jpg');               % salva images
                        imwrite(RGB3,'ScanLine.jpg');                   %
                        
                        cont = cont + 1;
                    else
                    
                        blankimage1 = imread('endPoints.jpg');          %
                        blankimage2 = imread('branchPoints.jpg');       % cria imagens   
                        blankimage3 = imread('ScanLine.jpg');           %
                        RGB1 = insertShape(blankimage1,'Line',[coordinates(i-1,1), coordinates(i-1,2), ...
                            coordinates(i,1), coordinates(i,2)],'LineWidth',3,'Color', {'white'});     % desenha linhas abertua 3 pixel
                        RGB2 = insertShape(blankimage2,'Line',[coordinates(i-1,1), coordinates(i-1,2), ...
                            coordinates(i,1), coordinates(i,2)],'LineWidth',1,'Color', {'white'});     % desenha linhas abertura 1 pixel
                        RGB3 = insertShape(blankimage3,'Line',[coordinates(i-1,1), coordinates(i-1,2), ...
                            coordinates(i,1), coordinates(i,2)],'LineWidth',round(fator*abertura*FatorReal),'Color', {'white'});    % desenha linha com largura de pixel de acordo com a abertura
                        imwrite(RGB1,'endPoints.jpg');                  % 
                        imwrite(RGB2,'branchPoints.jpg');               % salva images
                        imwrite(RGB3,'ScanLine.jpg');                   %
                        
                    end
                    %------------ salva comprimento -------------
                    
%                     if (coordinates(i-1,1)>0 && coordinates(i-1,1)<=T(2) && coordinates(i,1)>0 && coordinates(i,1)<=T(2)) && ...
%                         (coordinates(i-1,2)>0 && coordinates(i-1,2)<=T(1) && coordinates(i,2)>0 && coordinates(i,2)<=T(1))
                        
                    	XX = [coordinates(i-1,1),coordinates(i-1,2);coordinates(i,1),coordinates(i,2)];
                        distance = pdist(XX,'euclidean');
                        Abertura_Distancia(contDistancia,2) = distance;
                        pointsDistance(contDistancia) = contAbertura;
                        contDistancia = contDistancia + 1;
                            
%                     end
                    %------------ salva abertura -------------
                    
                    Abertura_Distancia(contAbertura,1) = abertura;
                    contAbertura = contAbertura + 1;
                    
                    %------------------------------------------

                    if button ~= 101 && button ~= 100
                        i = i + 1;
                        coordinates(i,:) = [x, y];
                        continue
                    end

                end
        end
    end
    
    Abertura_Distancia = Abertura_Distancia*FatorReal;
    Abertura_Distancia(Abertura_Distancia<=0)=[];
    
    for ijk=1:(length(Abertura_Distancia)/2)
        Area_Fracture = Area_Fracture + Abertura_Distancia(ijk)*Abertura_Distancia(ijk+(length(Abertura_Distancia))/2);
    end
    
    Area_Figure = T(1)*T(2)*FatorReal*FatorReal;
    
    P22 = Area_Fracture/Area_Figure;
    
    %%----------------------------------------------------------------------------------
    % salva variavel i para reutilização futura
    
    delete 'salva.txt';
    fileID = fopen('salva.txt','w');
    fprintf(fileID,'%i\n',i);
    fclose(fileID);

    %%----------------------------------------------------------------------------------
    % calcula angulo das fraturas e separa as familias

    k = i;
    Contfamilias = 1;
    Angle1 = atan((coordinates(2,2)-coordinates(1,2))/(coordinates(2,1)-coordinates(1,1)));

    for i=2:2:k
        Angle2 = atan((coordinates(i,2)-coordinates(i-1,2))/(coordinates(i,1)-coordinates(i-1,1)));                
        diff = (Angle2 - Angle1)*180/pi;

        if abs(diff)>90
            diff = 180-abs(diff); % coloca angulo no primeiro quadrante
        end
        
        if (coordinates(i-1,1)>0 && coordinates(i-1,1)<=T(2) && coordinates(i,1)>0 && coordinates(i,1)<=T(2))
        	if (coordinates(i-1,2)>0 && coordinates(i-1,2)<=T(1) && coordinates(i,2)>0 && coordinates(i,2)<=T(1))
        
                XX = [coordinates(i-1,1),coordinates(i-1,2);coordinates(i,1),coordinates(i,2)];
                distance = pdist(XX,'euclidean');
                
                for j=AngFamilias:AngFamilias:90
                    if diff>=-j && diff<-(j-AngFamilias) || diff>(j-AngFamilias) && diff<=j
                        familias(ContfamiliasVet(1,Contfamilias),Contfamilias) = distance;
                        ContfamiliasVet(1,Contfamilias) = ContfamiliasVet(1,Contfamilias) + 1;
                    end
                    Contfamilias = Contfamilias + 1;
                end  
                
            end
        end
              
        Contfamilias = 1;    
    end

    XX = [coordinates(1,1),coordinates(1,2);coordinates(2,1),coordinates(2,2)];
    distance = pdist(XX,'euclidean');
    familias(ContfamiliasVet(1,1),1) = distance;

    familias = familias*FatorReal;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------------------------------------------------------------------------
%   ScanLine Linear e Analise topologica

if opcao == 'B' || i~=0
    
    %%----------------------------------------------------------------------------------
    % cria uma nova imagem com apenas o contorno das fraturas

    blankimage3 = imread('ScanLine.jpg');
    I = rgb2gray(blankimage3);
    BW1 = edge(I,'sobel');
    B = medfilt2(BW1,[2 2]);
    imwrite(B,'SobelScan.jpg');

    %%----------------------------------------------------------------------------------
    % 
    
    blankimage3 = imread('ScanLine.jpg');
    figure, imshow(blankimage3);
    message = sprintf('L para scanline linear');
    uiwait(helpdlg(message));
    [x,y,button] = ginput(1);

    %%----------------------------------------------------------------------------------
    %   ScanLine Linear

    if button == 76  % 76 representa a letra L maiuscula (seleção scanline)
        Z = GUI_37();
        NumScanLines = str2num(Z);
        Contador = 1;

        %%----------------------------------------------------------------------------------
        % recebe os pontos de inicio e fim das scanlines

        if NumScanLines~=0
            hold on
            message = sprintf('Marque as ScanLines com o mouse (Linhas)');
            uiwait(helpdlg(message));
            for k=1:2*NumScanLines
                [x, y] = ginputc(1, 'Color', 'y', 'LineWidth', 1);
                coordinatesScanLine(k,:) = [x, y];
                if (mod(k,2)==0)
                    BW1 = imread('SobelScan.jpg');
                    H = line([coordinatesScanLine(k-1,1) coordinatesScanLine(k,1)],[coordinatesScanLine(k-1,2) coordinatesScanLine(k,2)],'Color','red','LineWidth',1);
                    Image = insertShape(BW1,'Line',[coordinatesScanLine(k-1,1), coordinatesScanLine(k-1,2), coordinatesScanLine(k,1), coordinatesScanLine(k,2)],...
                        'LineWidth',1,'Color', {'white'});
                    imwrite(Image,strcat(num2str(Contador),'.jpg'));
                    Contador = Contador + 1;
                end
            end

            ScanL = zeros(n,NumScanLines);
            preencheVetor = 1;
            ContLine = 1;
            ContScan = 1;

            pause(10);
            saveas(gcf,'Scan.png')
            close all
        %%----------------------------------------------------------------------------------
        % encontra pontos cruzados scanline e fraturas

            for m=1:NumScanLines
                folder = pwd;
                baseFileName = strcat(num2str(m),'.jpg');
                fullFileName = fullfile(folder, baseFileName);
                if ~exist(fullFileName, 'file')
                    fullFileNameOnSearchPath = baseFileName; 
                    if ~exist(fullFileNameOnSearchPath, 'file')
                        errorMessage = sprintf('Erro: %s não existe o arquivo na pasta.', fullFileName);
                        uiwait(warndlg(errorMessage));
                        return;
                    end
                end
                grayImage = imread(fullFileName);
                [rows, columns, numberOfColorChannels] = size(grayImage);
                if numberOfColorChannels > 1
                    grayImage = grayImage(:, :, 2); 
                end

                binaryImage = grayImage >= 128;
                binaryImage = binaryImage(4:end-4, 4:end-4);
                binaryImage = bwmorph(binaryImage, 'skel', inf);
                branchpoints = bwmorph(binaryImage, 'branchpoints');
                bridge = bwmorph(branchpoints, 'bridge');
                [rowsSL, columnsSL] = find(bridge);

            %%----------------------------------------------------------------------------------
            % tira duplicata de pontos cruzados

                CopyRowsSL = rowsSL;
                CopyColumnsSL = columnsSL;
                Conta = 0;
                ContaB = 0;
                Xvazio = zeros(1,length(rowsSL));
                Yvazio = zeros(1,length(rowsSL));

                for i = 1:(length(CopyRowsSL))
                    for k = i+1:(length(CopyRowsSL))
                        X =  [CopyRowsSL(i),CopyColumnsSL(i);CopyRowsSL(k),CopyColumnsSL(k)];
                        dist = pdist(X,'euclidean');
                        if dist<=5
                            Conta = Conta + 1;
                            vazio(Conta) = k;   
                        end
                    end
                    for j=1:Conta
                        CopyRowsSL(vazio(j)-ContaB) = [];
                        CopyColumnsSL(vazio(j)-ContaB) = [];
                        ContaB = ContaB + 1;
                    end
                    Conta = 0;
                    ContaB = 0;
                end
            %%----------------------------------------------------------------------------------
            % calcula os coeficientes da reta da scanline e retira pontos fora desta reta

            Conta = 0;
            ContaB = 0;

            coeficiente_A = (coordinatesScanLine(2*m-1,2) - coordinatesScanLine(2*m,2))/(coordinatesScanLine(2*m-1,1)-coordinatesScanLine(2*m,1));
            coeficiente_B = coordinatesScanLine(2*m-1,2) - coeficiente_A*coordinatesScanLine(2*m-1,1);

            for kl=1:length(CopyColumnsSL)
               yy = coeficiente_A*CopyColumnsSL(kl) + coeficiente_B;
               diferenca = yy-CopyRowsSL(kl);
               if abs(diferenca) <= 30
                  Conta = Conta + 1;
                  Xvazio(Conta) = CopyColumnsSL(kl);
                  Yvazio(Conta) = CopyRowsSL(kl);
               end
            end

            CopyColumnsSL = Xvazio;
            CopyRowsSL = Yvazio;

            CopyColumnsSL(CopyColumnsSL==0)=[];
            CopyRowsSL(CopyRowsSL==0)=[];

            if mod(length(CopyColumnsSL),2)~=0 || mod(length(CopyRowsSL),2)~=0
               fprintf('Problema na identificação dos pontos de cruzamento\n');
               continue 
            end

            %%----------------------------------------------------------------------------------
            % imprime os pontos cruzados encontrados na imagem para das
            % scanlines

            workspace;
            format long g;
            format compact;
            fontSize = 20;
            imshow(imcomplement(grayImage), []);
            
            axis on;
            hold on
            title('Padrao X', 'FontSize', fontSize, 'Interpreter', 'None');
            for k = 1 : length(CopyColumnsSL);
                plot(CopyColumnsSL(k), CopyRowsSL(k), 'bd', 'MarkerSize', 12, 'LineWidth', 1);
            end
            pause(10)
            close all

            %%----------------------------------------------------------------------------------
            % calcula a abertura e a distancia entre as fraturas

            ContBP = 1;
            if coordinatesScanLine(m,1)>=0 || coordinatesScanLine(m,2)>=0
                for n=1:(length(CopyColumnsSL)+1)
                    if n==1 || n==(length(CopyColumnsSL)+1)
                        X =  [coordinatesScanLine(ContScan,1),coordinatesScanLine(ContScan,2);CopyColumnsSL(ContBP),CopyRowsSL(ContBP)];
                        dS = pdist(X,'euclidean');
                        ScanL(preencheVetor,m) = dS;
                        preencheVetor = preencheVetor + 1;
                        ContScan = ContScan + 1;
                    else
                        X =  [CopyColumnsSL(ContBP),CopyRowsSL(ContBP);CopyColumnsSL(ContBP+1),CopyRowsSL(ContBP+1)];
                        dA = pdist(X,'euclidean');
                        ScanL(preencheVetor,m) = dA;
                        preencheVetor = preencheVetor + 1;
                        ContBP = ContBP + 1;
                    end
                end            
            end

            preencheVetor = 1;

            end
        end
    pause(5);    
    end

    %%----------------------------------------------------------------------------------
    % salva arquivo de texto da scanline (aberturas,distancias)

    ScanL = ScanL*FatorReal;
    SizeScanL = size(ScanL);

    for t = 1:SizeScanL(2)

       ScanLTemp = ScanL(:,t);

       if ScanLTemp(1) == 0
          continue  
       else
          ScanLTemp = ScanLTemp(ScanLTemp>0.1);
          ScanLTemp(1) = [];

          z=[];
          x=[];
          indicex = 1;
          indicez = 1;
          contadorxy = 1;

          for n = 1:length(ScanLTemp)
              if mod(n,2)==1
                x(indicex) = ScanLTemp(contadorxy);
                indicex = indicex + 1;
              else
                z(indicez) = ScanLTemp(contadorxy);
                indicez = indicez + 1;
              end
              contadorxy = contadorxy + 1;
          end
          
          if length(x)==length(z)
              A = [z;x];
              fileID = fopen(strcat(num2str(t),'.txt'),'w');
              fprintf(fileID,'%4.4f %4.4f\r\n',A);
              fclose(fileID);
              
              %%----------------------------------------------------------------------------------
              % Função que plota lei de potência
              scanlinev31(num2str(t));
          end
          
       end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%----------------------------------------------------------------------------------
    % Analise topologica

    clc;        
    close all;  
    workspace;  
    format long g;
    format compact;
    fontSize = 20;
    ContX = 0;
    ContI = 0;

    fileID = fopen('salva.txt','r');
    i = fscanf(fileID,'%i');
    ContLines = i;

    %%----------------------------------------------------------------------------------
    % carrega a imagem e encontra o cruzamento das fraturas

    folder = pwd;
    baseFileName = 'branchPoints.jpg';
    fullFileName = fullfile(folder, baseFileName);

    if ~exist(fullFileName, 'file')
        fullFileNameOnSearchPath = baseFileName; 
        if ~exist(fullFileNameOnSearchPath, 'file')
            errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
            uiwait(warndlg(errorMessage));
            return;
        end
    end
    grayImage = imread(fullFileName);

    [rows, columns, numberOfColorChannels] = size(grayImage);
    if numberOfColorChannels > 1
        grayImage = grayImage(:, :, 2); 
    end

    binaryImage = grayImage >= 128;
    binaryImage = binaryImage(4:end-4, 4:end-4);
    binaryImage = bwmorph(binaryImage, 'skel', inf);

    branchpoints = bwmorph(binaryImage, 'branchpoints');
    [rowsBP, columnsBP] = find(branchpoints);

    %%----------------------------------------------------------------------------------
    % tira duplicata de pontos cruzados 

    CopyRowsBP = rowsBP;
    CopyColumnsBP = columnsBP;
    Conta = 0;
    ContaB = 0;
    vazio = zeros(1,length(rowsBP));

    for i = 1:(length(CopyRowsBP))
       for k = i+1:(length(CopyRowsBP))
           modulo = sqrt((CopyRowsBP(i)-CopyRowsBP(k))^2+(CopyColumnsBP(i)-CopyColumnsBP(k))^2);
           if modulo<=10
               Conta = Conta + 1;
               vazio(Conta) = k;   
           end
       end
       for j=1:Conta
           CopyRowsBP(vazio(j)-ContaB) = [];
           CopyColumnsBP(vazio(j)-ContaB) = [];
           ContaB = ContaB + 1;
       end
       Conta = 0;
       ContaB = 0;
    end

    %%----------------------------------------------------------------------------------
    % carrega a imagem e encontra o fim das fraturas

    folder = pwd;
    baseFileName = 'endPoints.jpg';
    fullFileName = fullfile(folder, baseFileName);
    if ~exist(fullFileName, 'file')
        fullFileNameOnSearchPath = baseFileName; % No path this time.
        if ~exist(fullFileNameOnSearchPath, 'file')
            errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
            uiwait(warndlg(errorMessage));
            return;
        end
    end
    grayImage = imread(fullFileName);
    [rows, columns, numberOfColorChannels] = size(grayImage);
    if numberOfColorChannels > 1
        grayImage = grayImage(:, :, 2); 
    end

    binaryImage = grayImage >= 128;
    binaryImage = binaryImage(4:end-4, 4:end-4);
    binaryImage = bwmorph(binaryImage, 'skel', inf);

    bridge = bwmorph(binaryImage, 'bridge');
    endpoints = bwmorph(bridge, 'endpoints');

    [rowsEP, columnsEP] = find(endpoints);

    %%----------------------------------------------------------------------------------
    % tira duplicata de fim de pontos

    CopyRowsEP = rowsEP;
    CopyColumnsEP = columnsEP;
    Conta = 0;
    ContaB = 0;
    vazio = zeros(1,length(rowsEP));

    for i = 1:(length(CopyRowsEP))
       for k = i+1:(length(CopyRowsEP))
           modulo = sqrt((CopyRowsEP(i)-CopyRowsEP(k))^2+(CopyColumnsEP(i)-CopyColumnsEP(k))^2);
           if modulo<=10
               Conta = Conta + 1;
               vazio(Conta) = k;   
           end
       end
       for j=1:Conta
           CopyRowsEP(vazio(j)-ContaB) = [];
           CopyColumnsEP(vazio(j)-ContaB) = [];
           ContaB = ContaB + 1;
       end
       Conta = 0;
       ContaB = 0;
    end

    %%----------------------------------------------------------------------------------
    % conta pontos fora da imagem 

    DotOut = 0;

    for k=1:length(CopyRowsEP)
        if (CopyRowsEP(k) < 10 || CopyRowsEP(k) > T(1)-10 || CopyColumnsEP(k) < 10 || CopyColumnsEP(k) > T(2)-10)
            DotOut = DotOut + 1;
        end
    end
    
    %%----------------------------------------------------------------------------------
    % tira pontos em comum I e X

    Conta = 0;
    ContaB = 0;
    vazio = zeros(1,length(CopyRowsEP));

    for i = 1:(length(CopyRowsBP))
       for k = 1:(length(CopyRowsEP))
           modulo = sqrt((CopyRowsBP(i)-CopyRowsEP(k))^2+(CopyColumnsBP(i)-CopyColumnsEP(k))^2);
           if modulo<=50
               Conta = Conta + 1;
               vazio(Conta) = k;   
           end
       end
       for j=1:Conta
           CopyRowsEP(vazio(j)-ContaB) = [];
           CopyColumnsEP(vazio(j)-ContaB) = [];
           ContaB = ContaB + 1;
       end
       Conta = 0;
       ContaB = 0;
    end

    [sharedVals,idxsIntoCopyColumnsEP] = intersect(CopyColumnsEP,CopyColumnsBP);
    for i=1:length(idxsIntoCopyColumnsEP)
        CopyColumnsEP(idxsIntoCopyColumnsEP(i))=[];
        CopyRowsEP(idxsIntoCopyColumnsEP(i))=[];
    end

    %%------------------------------------------------------------------------------------
    % conta os elementos de X e Y

    ContY = 0;
    ContX = 0;
    indiciesY = zeros(n,1);
    indiciesX = zeros(n,1);

    for i = 1:length(CopyColumnsBP)
        teste = 0;
        for k = 1:ContLines
          X =  [CopyColumnsBP(i),CopyRowsBP(i);coordinates(k,1),coordinates(k,2)];
          d = pdist(X,'euclidean');
          if d<=50
             ContY = ContY + 1; 
             indiciesY(ContY) = i;
             teste = 1;
          end
       end
       if teste == 0
           ContX = ContX + 1;
           indiciesX(ContX) = i;
       end
    end

    indiciesX = unique(indiciesX);
    if indiciesX(1)==0
        indiciesX(1) = [];
    end

    %%----------------------------------------------------------------------------------
    % imprime imagem com os padroes topologicos

    Image = imread('endPoints.jpg');
    bw2 = imcomplement(Image);
    imshow(bw2, []);
    axis on;
    title('Padrões Topológicos', 'FontSize', fontSize, 'Interpreter', 'None');
    hold on;

    for k = 1 : ContY
        plot(CopyColumnsBP(indiciesY(k)), CopyRowsBP(indiciesY(k)), 'r^', 'MarkerSize', 12, 'LineWidth', 2);
    end

    for k = 1 : ContX
        plot(CopyColumnsBP(indiciesX(k)), CopyRowsBP(indiciesX(k)), 'bd', 'MarkerSize', 12, 'LineWidth', 2);
    end

    for k = 1 : length(CopyColumnsEP)  % imprime apenas pontos de fim dentro da imagem

        if (CopyRowsEP(k) <= 8 || CopyRowsEP(k) >= T(1)-8 || CopyColumnsEP(k) <= 8 || CopyColumnsEP(k) >= T(2)-8)
        else
            plot(CopyColumnsEP(k), CopyRowsEP(k), 'go', 'MarkerSize', 12, 'LineWidth', 2);
        end
    end
    
    saveas(gcf,'Topologico.png');
    % ImageTop = imread('Topologico.png');
    % bw2 = imcomplement(ImageTop);
    % imwrite(bw2,'Topologico1.png');
    pause(5)
    
    %%----------------------------------------------------------------------------------
    % conta elementos de I dentro da imagem

    ContI = length(CopyRowsEP);
    ContI = ContI - DotOut;
    if ContI<0
        ContI = 0;
    end

    %%----------------------------------------------------------------------------------
    % diagrama ternário

    figure
    [h,hg,htick]=terplot;
    title('Diagrama Ternário');
    hter=ternaryc(ContY/(ContX+ContY+ContI),ContX/(ContX+ContY+ContI),ContI/(ContX+ContY+ContI)); % Y, X, I
    set(hter,'marker','diamond','markerfacecolor','blue','markersize',8)
    hlabels=terlabel('Y','X','I');
    
    saveas(gcf,'Ternario.png')
    pause(5)
end

%%----------------------------------------------------------------------------------
% gráficos 

% clc
% close all
% 
% xreal = [0.1 0.2 0.3 0.4 0.5]; % cm
% yreal = [3.1 6.3 12.5 25 50]; % cm
% 
% Dist = Abertura_Distancia(:,2);
% Abert = Abertura_Distancia(:,1);
% 
% pointsDistance = pointsDistance(pointsDistance>0);
% xaprox = Abert(pointsDistance);
% yaprox = Dist(Dist>0);
% 
% Erro_Abertura1 = zeros(1,length(xaprox));
% Erro_Comprimento1 = zeros(1,length(yaprox));
% 
% Erro_Abertura2 = zeros(1,length(xaprox));
% Erro_Comprimento2 = zeros(1,length(yaprox));
% 
% Erro_Abertura3 = zeros(1,length(xaprox));
% Erro_Comprimento3 = zeros(1,length(yaprox));
% 
% Erro_Abertura4 = zeros(1,length(xaprox));
% Erro_Comprimento4 = zeros(1,length(yaprox));
% 
% Erro_Abertura5 = zeros(1,length(xaprox));
% Erro_Comprimento5 = zeros(1,length(yaprox));
% 
% Cont_A1 = 1;
% Cont_A2 = 1;
% Cont_A3 = 1;
% Cont_A4 = 1;
% Cont_A5 = 1;
% Cont_L1 = 1;
% Cont_L2 = 1;
% Cont_L3 = 1;
% Cont_L4 = 1;
% Cont_L5 = 1;
% 
% for r=1:length(xaprox)
%     
%     % Erro das aberturas
%     
%     if xaprox(r)>0 && xaprox(r)<=0.15
%     	Erro_Abertura1(Cont_A1) = (xreal(1)-xaprox(r))/xreal(1);
%         Cont_A1 = Cont_A1 + 1;
%     end
%     
%     if xaprox(r)>0.15 && xaprox(r)<=0.25 
%     	Erro_Abertura2(Cont_A2) = (xreal(2)-xaprox(r))/xreal(2);
%         Cont_A2 = Cont_A2 + 1;
%     end
%     
%     if xaprox(r)>0.25 && xaprox(r)<=0.35 
%     	Erro_Abertura3(Cont_A3) = (xreal(3)-xaprox(r))/xreal(3);
%         Cont_A3 = Cont_A3 + 1;
%     end
%     
%     if xaprox(r)>0.35 && xaprox(r)<=0.45 
%     	Erro_Abertura4(Cont_A4) = (xreal(4)-xaprox(r))/xreal(4);
%         Cont_A4 = Cont_A4 + 1;
%     end
%     
%     if xaprox(r)>0.45 && xaprox(r)<=0.7
%     	Erro_Abertura5(Cont_A5) = (xreal(5)-xaprox(r))/xreal(5);
%         Cont_A5 = Cont_A5 + 1;
%     end
%     
% end
% 
% for r=1:length(yaprox)
%     
%    % Erros dos comprimentos 
%     
%     if yaprox(r)>0 && yaprox(r)<=4.7
%         Erro_Comprimento1(Cont_L1) = (yreal(1)-yaprox(r))/yreal(1);
%         Cont_L1 = Cont_L1 + 1;
%     end
%     
%     if yaprox(r)>4.7 && yaprox(r)<=9.4
%         Erro_Comprimento2(Cont_L1) = (yreal(2)-yaprox(r))/yreal(2);
%         Cont_L2 = Cont_L2 + 1;
%     end
%     
%     if yaprox(r)>9.4 && yaprox(r)<=18.75
%         Erro_Comprimento3(Cont_L1) = (yreal(3)-yaprox(r))/yreal(3);
%         Cont_L3 = Cont_L3 + 1;
%     end
%     
%     if yaprox(r)>18.75 && yaprox(r)<=37.5
%         Erro_Comprimento4(Cont_L4) = (yreal(4)-yaprox(r))/yreal(4);
%         Cont_L4 = Cont_L4 + 1;
%     end
%     
%     if yaprox(r)>37.5 && yaprox(r)<=50 
%         Erro_Comprimento5(Cont_L5) = (yreal(5)-yaprox(r))/yreal(5);
%         Cont_L5 = Cont_L5 + 1;
%     end
%     
% end
% 
% Erro_Abertura1(Erro_Abertura1==0) = [];
% Erro_Abertura2(Erro_Abertura2==0) = [];
% Erro_Abertura3(Erro_Abertura3==0) = [];
% Erro_Abertura4(Erro_Abertura4==0) = [];
% Erro_Abertura5(Erro_Abertura5==0) = [];
% 
% Erro_Comprimento1(Erro_Comprimento1==0) = [];
% Erro_Comprimento2(Erro_Comprimento2==0) = [];
% Erro_Comprimento3(Erro_Comprimento3==0) = [];
% Erro_Comprimento4(Erro_Comprimento4==0) = [];
% Erro_Comprimento5(Erro_Comprimento5==0) = [];
% 
% %%----------------------------------------------------------------------------------
% % boxplot
% 
% Error_AB = [Erro_Abertura1'; Erro_Abertura2'; Erro_Abertura3'; Erro_Abertura4'; Erro_Abertura5'];
% g1 = [ones(size(Erro_Abertura1')); 2*ones(size(Erro_Abertura2')); 3*ones(size(Erro_Abertura3'));4*ones(size(Erro_Abertura4'));5*ones(size(Erro_Abertura5'))];
% boxplot(Error_AB,g1);
% xlabel('Aberturas')
% ylabel('Erro Relativo (%)')
% title('Erro Relativo das Aberturas')
% saveas(gcf,'BoxplotAbertura.png')
% pause(5)
% 
% Error_L = [Erro_Comprimento1'; Erro_Comprimento2'; Erro_Comprimento3'; Erro_Comprimento4'; Erro_Comprimento5'];
% g2 = [ones(size(Erro_Comprimento1')); 2*ones(size(Erro_Comprimento2')); 3*ones(size(Erro_Comprimento3'));4*ones(size(Erro_Comprimento4'));5*ones(size(Erro_Comprimento5'))];
% boxplot(Error_L,g2);
% xlabel('Comprimentos')
% ylabel('Erro Relativo (%)')
% title('Erro Relativo dos Comprimentos')
% saveas(gcf,'BoxplotComprimento.png')
% pause(5)
% 
% %%----------------------------------------------------------------------------------
% % gráficos Comprimeto x abertura
% 
% plot(xreal,yreal,'r*');
% grid on;
% title('Curva Teorica');
% ylabel('Length, L (cm)')
% xlabel('Aperture, b (cm)')
% p1 = polyfit(xreal,yreal,1);
% y1 = polyval(p1,xreal);
% hold on; plot(xreal,y1);
% 
% angular = (y1(2) - y1(1))/(xreal(2)-xreal(1));
% linear = y1(1) - angular*xreal(1);
% 
% N = length(xreal); 
% r2 = 1 - N*sum((yreal - y1).^2)/(N*sum(yreal.^2) - sum(yreal)^2);
% hold on; plot(0.5,0,'.w','LineWidth',0.1);
% w2 = legend('Data',sprintf('L={%.3f}*b - %.3f{} ',angular,linear*(-1)),sprintf('R^2=%.3f{}',r2));
% set(w2,'FontSize',15);
% 
% pause(10);
% saveas(gcf,'grafico1.png')
% close all
% %%
% 
% plot(xaprox,yaprox,'r*');
% grid on;
% title('Curva Dados');
% ylabel('Length, L (cm)')
% xlabel('Aperture, b (cm)')
% p1 = polyfit(xaprox,yaprox,1);
% y1 = polyval(p1,xaprox);
% hold on; plot(xaprox,y1);
% 
% ang = (y1(2) - y1(1))/(xaprox(2)-xaprox(1));
% lin = y1(1) - ang*xaprox(1);
% 
% N1 = length(xaprox); 
% r2 = 1 - N1*sum((yaprox - y1).^2)/(N1*sum(yaprox.^2) - sum(yaprox)^2);
% hold on; plot(0.5,0,'.w','LineWidth',0.1);
% if lin>0
%     w3 = legend('Data',sprintf('L={%.3f}*b + %.3f{} ',ang,lin),sprintf('R^2=%.3f{}',r2));
%     set(w3,'FontSize',15);
% else
%     w3 = legend('Data',sprintf('L={%.3f}*b - %.3f{} ',ang,lin*(-1)),sprintf('R^2=%.3f{}',r2));
%     set(w3,'FontSize',15);
% end
% 
% pause(10);
% saveas(gcf,'grafico2.png')
% close all
