% plota pontos experimentais

clc
close all

xreal = [0.1 0.2 0.3 0.4 0.5]; % cm
yreal = [4.2 6.1 9.5 15.6 29.7]; % cm

Dist = Abertura_Distancia(:,2);
Abert = Abertura_Distancia(:,1);

xaprox = Abert(Abert>0);
yaprox = Dist(Dist>0);

Erro_Abertura1 = zeros(1,length(xaprox));
Erro_Comprimento1 = zeros(1,length(xaprox));

Erro_Abertura2 = zeros(1,length(xaprox));
Erro_Comprimento2 = zeros(1,length(xaprox));

Erro_Abertura3 = zeros(1,length(xaprox));
Erro_Comprimento3 = zeros(1,length(xaprox));

Erro_Abertura4 = zeros(1,length(xaprox));
Erro_Comprimento4 = zeros(1,length(xaprox));

Erro_Abertura5 = zeros(1,length(xaprox));
Erro_Comprimento5 = zeros(1,length(xaprox));

Cont_1 = 1;
Cont_2 = 1;
Cont_3 = 1;
Cont_4 = 1;
Cont_5 = 1;

for r=1:length(xaprox)
    
    if yaprox(r)>0 && yaprox(r)<=5.15
    	Erro_Abertura1(Cont_1) = (xreal(1)-xaprox(r))/xreal(1);
        Erro_Comprimento1(Cont_1) = (yreal(1)-yaprox(r))/yreal(1);
        Cont_1 = Cont_1 + 1;
    end
    if yaprox(r)>5.15 && yaprox(r)<=7.8
    	Erro_Abertura2(Cont_2) = (xreal(2)-xaprox(r))/xreal(2);
        Erro_Comprimento2(Cont_2) = (yreal(2)-yaprox(r))/yreal(2);
        Cont_2 = Cont_2 + 1;
    end
    if yaprox(r)>7.8 && yaprox(r)<=12.55
    	Erro_Abertura3(Cont_3) = (xreal(3)-xaprox(r))/xreal(3);
        Erro_Comprimento3(Cont_3) = (yreal(3)-yaprox(r))/yreal(3);
        Cont_3 = Cont_3 + 1;
    end
    if yaprox(r)>12.55 && yaprox(r)<=22.65
    	Erro_Abertura4(Cont_4) = (xreal(4)-xaprox(r))/xreal(4);
        Erro_Comprimento4(Cont_4) = (yreal(4)-yaprox(r))/yreal(4);
        Cont_4 = Cont_4 + 1;
    end
    if yaprox(r)>22.65 && yaprox(r)<=1000
    	Erro_Abertura5(Cont_5) = (xreal(5)-xaprox(r))/xreal(5);
        Erro_Comprimento5(Cont_5) = (yreal(5)-yaprox(r))/yreal(5);
        Cont_5 = Cont_5 + 1;
    end
    
end

Erro_Abertura1(Erro_Abertura1==0) = [];
Erro_Abertura2(Erro_Abertura2==0) = [];
Erro_Abertura3(Erro_Abertura3==0) = [];
Erro_Abertura4(Erro_Abertura4==0) = [];
Erro_Abertura5(Erro_Abertura5==0) = [];

Erro_Comprimento1(Erro_Comprimento1==0) = [];
Erro_Comprimento2(Erro_Comprimento2==0) = [];
Erro_Comprimento3(Erro_Comprimento3==0) = [];
Erro_Comprimento4(Erro_Comprimento4==0) = [];
Erro_Comprimento5(Erro_Comprimento5==0) = [];

Error_AB = [Erro_Abertura1'; Erro_Abertura2'; Erro_Abertura3'; Erro_Abertura4'; Erro_Abertura5'];
g1 = [ones(size(Erro_Abertura1')); 2*ones(size(Erro_Abertura2')); 3*ones(size(Erro_Abertura3'));4*ones(size(Erro_Abertura4'));5*ones(size(Erro_Abertura5'))];
boxplot(Error_AB,g1);
xlabel('Aberturas')
ylabel('Erro Relativo (%)')
title('Erro Relativo das Aberturas')
saveas(gcf,'BoxplotAbertura.png')
pause(5)

Error_L = [Erro_Comprimento1'; Erro_Comprimento2'; Erro_Comprimento3'; Erro_Comprimento4'; Erro_Comprimento5'];
g2 = [ones(size(Erro_Comprimento1')); 2*ones(size(Erro_Comprimento2')); 3*ones(size(Erro_Comprimento3'));4*ones(size(Erro_Comprimento4'));5*ones(size(Erro_Comprimento5'))];
boxplot(Error_L,g2);
xlabel('Comprimentos')
ylabel('Erro Relativo (%)')
title('Erro Relativo dos Comprimentos')
saveas(gcf,'BoxplotComprimento.png')
pause(5)

plot(xreal,yreal,'r*');
grid on;
title('Curva Teorica');
ylabel('Length, L (cm)')
xlabel('Aperture, b (cm)')
p1 = polyfit(xreal,yreal,1);
y1 = polyval(p1,xreal);
hold on; plot(xreal,y1);

angular = (y1(2) - y1(1))/(xreal(2)-xreal(1));
linear = y1(1) - angular*xreal(1);

N = length(xreal); 
r2 = 1 - N*sum((yreal - y1).^2)/(N*sum(yreal.^2) - sum(yreal)^2);
hold on; plot(0.5,0,'.w','LineWidth',0.1);
w2 = legend('Data',sprintf('L={%.3f}*b - %.3f{} ',angular,linear*(-1)),sprintf('R^2=%.3f{}',r2));
set(w2,'FontSize',15);

saveas(gcf,'grafico1.png')
close all

plot(xaprox,yaprox,'r*');
grid on;
title('Curva Dados');
ylabel('Length, L (cm)')
xlabel('Aperture, b (cm)')
p1 = polyfit(xaprox,yaprox,1);
y1 = polyval(p1,xaprox);
hold on; plot(xaprox,y1);

ang = (y1(2) - y1(1))/(xaprox(2)-xaprox(1));
lin = y1(1) - ang*xaprox(1);

N1 = length(xaprox); 
r2 = 1 - N1*sum((yaprox - y1).^2)/(N1*sum(yaprox.^2) - sum(yaprox)^2);
hold on; plot(0.5,0,'.w','LineWidth',0.1);
if lin>0
    w3 = legend('Data',sprintf('L={%.3f}*b + %.3f{} ',ang,lin),sprintf('R^2=%.3f{}',r2));
    set(w3,'FontSize',15);
else
    w3 = legend('Data',sprintf('L={%.3f}*b - %.3f{} ',ang,lin*(-1)),sprintf('R^2=%.3f{}',r2));
    set(w3,'FontSize',15);
end

saveas(gcf,'grafico2.png')
close all

