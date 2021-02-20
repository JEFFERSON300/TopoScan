%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obten as Scanline calculando os coeficientes da lei de Potencia 
%                         Projeto - Turing 2012
% Autor: Rafael F.V.C. Santos
% Atualizado: 03/02/13 - adicao de grafico com escala log-log
% ULTIMA ATUALIZAÇÃO: 10/03/14   - adicao de formato e resolucao de figs;
% Versao: 3.1
% alteracao do tamanho da scanline;
% 
% colocar o numero de dados
% O valor maximo e minimo de dados
% criar rotina para salvar dados em arquivo
% Fazer graficos e salvar em formato adequado
% colocar quantos dados foram excluidos - pq se repetiram
% Add mod tamanho das legendas
% Jah plota o r^2 no grafico

% OBS: para testar os dados das scanlines artificiais do matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = scanlinev31(name)

% clear all; clc;
% close all;

%--------------------------------------------------------------------------
s_font = 15;
s_leg = 15;
cor_data = 'bo' ;
cor_fit = 'g' ;
tik = 1;

resFig = '-r150';  
formFig = '-dtiff';

tag = 'SCL'; % indicativo para mudança de nome de figuras salvas

%--------------------------------------------------------------------------

% Carrega dados do diretorio atual
DD = load (strcat(name,'.txt')); % <---------------------------------------------|*AQUI*|

dados = DD; % <---------------------------------------------|*AQUI*|
[m, n]= size(dados); % <-----------------------------------------|*AQUI*|

dist = dados(:,1); % distancia entre aberturas
Ab = dados(:,2); % abertura da fratura

% Dados em milimetros
unid =1; % 1 para dados em (mm) e 0 para dados em (m)
% colocar termos cruzados. ex. dist(mm) e Ab(m) ou vice versa


nome='Artigo';  % <---------------------------------------------|*AQUI*|

disp('===================================================================')
disp('------------CONSTRUÇÃO AUTOMÁTICA DA LEI DE POTÊNCIA---------------')
disp('===================================================================')
disp('                                                                   '); 

% ---- -----------------------UNIDADES DE MEDIDAS--------------------------
disp('------------INFORMAÇÕES SOBRE A SCANLINE EXPERIMENTAL--------------')

if max(dados(:,1)) > 100
    disp('Unidade dos dados: ');
    disp('Os dados fornecidos estao em milimetros')
end

% tamanho da scanline teorica
% scl_T=input('Digite o tamanho da scanline em metros = ');
scl_T=58; % <----------------------------------------------------|*AQUI*|

if unid == 1
disp('-------------------Dados estao em milimetros-----------------------');
    disp('O tamanho da scanline (m) é: ');
    scl_E = (sum(dist)+sum(Ab))/1000 % OBS: MODIFICADO PQ dist(m) e Ab(mm)
end

if unid == 0
disp('------------------------Dados estao em metros----------------------');
    disp('O tamanho da scanline (m) é: ');
    scl_E = (sum(dist)+sum(Ab))*1000 
end

if scl_T ~= scl_E  
disp('O erro do comprimento da scanline em % é: '); 
err = (abs(scl_E-scl_T)/max([scl_E scl_T]))*100
end

disp('-------------------------------------------------------------------')

% vertor de Classes
classe=sort(randperm(m))'; % classe=1:m;

% abertura ordenada
Ab_ord = sort(Ab,'descend');

% Frequencia Acumulada
Fr_ac = classe/scl_E;

v=[classe Fr_ac Ab_ord];

% Retira valores repetidos segundo Ortega et al., 2006
[nao ii]=unique(v(:,3));
[n_repete aux]=size(ii);
disp('Numero de dados retirados :');
n_retirado = m - n_repete
ii=sort(ii,'ascend');

disp('-------------------------------------------------------------------')
disp(' Dados antes da retirada de classes repetidas : ')
v
disp('-------------------------------------------------------------------')

disp('-------------------------------------------------------------------')
disp(' Dados apos a retirada de classes repetidas : ')
d=[classe(ii) Fr_ac(ii) sort(nao,'descend')]
disp('-------------------------------------------------------------------')

% Salva dados efetivos da plotagem da scanline
%--------------------------------------------------------------------------
rf=['scanl_' nome '.dat'];
fid1=fopen(rf,'w');
for i=1:max(size(d))
fprintf(fid1,'%7.5g  %7.5g  %7.5g\n',d(i,1),d(i,2),d(i,3));
end
fclose(fid1);
%--------------------------------------------------------------------------

x = log10(d(:,3));  % aberturas
y = log10(d(:,2));  % frequencia acumulada

% % [b,bint,R,rint,satars]=regress(x,y)

% Faz ajuste dos falores a uma funcao potencia

[p s]=polyfit(x,y,1);

[z delta]=polyval(p,x,s);

plot(x,y,cor_data);

axis auto
grid 
xlabel('Aperture, b (mm)','fontsize',s_font); 
ylabel('Cumulative Frequency, F (m^{-1})','fontsize',s_font);
title('Log-Log');
%--------------------------------------------------------------------------
% Calculo do R^2
disp(' Coeficiente de Regressão Linear: ');
y1 = polyval(p,x);
N=length(x); 
r2 = 1 - N*sum((y - y1).^2)/(N*sum(y.^2) - sum(y)^2)
disp('-------------------------------------------------------------------')
disp('----------------------------FIM DA ANÁLISE-------------------------')
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------


k=p(1);
loga=p(2);
a=10^(loga);
% hold on; 
% plot(x,k*x+loga,'g','LineWidth',1);hold on; plot(0,0,'.w','LineWidth',0.1);
% w1=legend('Data',sprintf('F=%.3f{}log(b)+log(%.3f)',k,a),sprintf('R^2=%.3f{}',r2));
% set(w1,'FontSize',s_leg);  %<------------------------- MODIFICA Tamanho da legenda!!
% 
% print(formFig,resFig,[strcat('Artigo_log_',name) tag]);
% % saveas(gcf, ['Artigo_log_' tag], 'fig');


close

%--------------------------------------------------------------------------

% Grafico LogLog
figure

loglog(d(:,3),d(:,2),'bo');
xlabel('Aperture, b (mm)','fontsize',s_font); 
ylabel('Cumulative Frequency, F (m^{-1})','fontsize',s_font);
title('Log-Log');
hold on; plot(d(:,3),a*d(:,3).^k,'k','LineWidth',1.5);
hold on; plot(0,0,'.w','LineWidth',0.1);
w2=legend('Data',sprintf('F=%.3f{} b^{%.3f}',a,k),sprintf('R^2=%.3f{}',r2));
set(w2,'FontSize',s_leg);
grid on; 


% print('-dtiff','-r300','Artigo_loglog_gridOn');
print(formFig,resFig,[strcat('Artigo_loglog_gridOn_',name) tag]);
% saveas(gcf, ['Artigo_loglog_gridOn_' tag], 'fig');

% figure
% % gráfico sem grid
% loglog(d(:,3),d(:,2),'bo');
% xlabel('Aperture, b (mm)','fontsize',s_font);
% ylabel('Cumulative Frequency, F (m^{-1})','fontsize',s_font);
% % title('Escala Micro');
% hold on; plot(d(:,3),a*d(:,3).^k,'k','LineWidth',1.5);
% hold on; plot(0,0,'.w','LineWidth',0.1);
% w3=legend('Data',sprintf('F=%.3f{} b^{%.3f}',a,k),sprintf('R^2=%.3f{}',r2));
% set(w3,'FontSize',s_leg);
% 
% % print('-dtiff','-r300','Artigo_loglog_gridOff');
% print(formFig,resFig,[strcat('Artigo_loglog_gridOff_',name) tag]);
% % saveas(gcf, ['Artigo_loglog_gridOff_' tag], 'fig');
%--------------------------------------------------------------------------

figure 
plot(d(:,3),d(:,2),cor_data);
xlabel('Aperture, b (mm)','fontsize',s_font); 
ylabel('Cumulative Frequency, F (m^{-1})','fontsize',s_font);
title('Linear Axes');
axis equal square
grid
hold on; plot(d(:,3),a*d(:,3).^k,'g');hold on; plot(0,0,'.w','LineWidth',0.1);
w4=legend('Data',sprintf('F=%.3f{} b^{%.3f}',a,k),sprintf('R^2=%.3f{}',r2));
set(w4,'FontSize',s_leg);

% print('-dtiff','-r300','Artigo_morm');
print(formFig,resFig,[strcat('Artigo_morm_',name) tag]);
% saveas(gcf, ['Artigo_morm_' tag], 'fig');


% Plotagem da Roseta
% figure
% theta = 0;
% r =(2*pi);
% polar(theta,r);
% imagesc(d) % coloca uma image colorida correspondentes a valores
% colobar;

%--------------------------------------------------------------------------
%-------------------SALVAR INFORMAÇOES EM ARQUIVO DE TEXTO-----------------

nf=num2str(m,5);
rf=['INFO_' nome nf  '.dat'];

%  Salvar arquivo .dat 
fid=fopen(rf,'w+');
fprintf(fid,'%c','Informacoes dos dados das scanlines:  ');
fprintf(fid,'%c',date);
fprintf(fid,'\n'); fprintf(fid,'\n');

fprintf(fid,'%c','-----------------------------------------------------------');
fprintf(fid,'\n');
fprintf(fid,'%c','Numero de dados de campo = ');
fprintf(fid,'%d',m);fprintf(fid,'\n');fprintf(fid,'\n');

fprintf(fid,'%c','Numero de dados retirados = ');
fprintf(fid,'%d',n_retirado); fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%c','-----------------------------------------------------------');
fprintf(fid,'\n');
fprintf(fid,'%c','Tamanho da Scanline (m) = ');
fprintf(fid,'%d',scl_E);
fprintf(fid,'\n');fprintf(fid,'\n');
fprintf(fid,'%c','Somatório das distancias [total] (mm) = ');
fprintf(fid,'%d',sum(dist));
fprintf(fid,'\n');fprintf(fid,'\n');
fprintf(fid,'%c','Somatório das Aberturas [total] (mm) = ');
fprintf(fid,'%d',sum(Ab));
fprintf(fid,'\n');fprintf(fid,'\n');
fprintf(fid,'Somatório das Aberturas [Lei de Potencia] (mm) = ');
fprintf(fid,'%d',sum(d(:,3)));
fprintf(fid,'\n');
fprintf(fid,'%c','-----------------------------------------------------------');

fprintf(fid,'\n'); fprintf(fid,'\n');
fprintf(fid,'%c',' a =        k='); fprintf(fid,'\n');
fprintf(fid,'%5.3f     %5.3f',a, k);
fprintf(fid,'\n'); fprintf(fid,'\n');
fprintf(fid,'%c','y = a * X^(k)');
fprintf(fid,'\n');fprintf(fid,'\n');
fprintf(fid,'%c','Coeficiente de Correlacao:  ');
fprintf(fid,'%c','R^2 = ');
fprintf(fid,'%5.3f',r2);
fprintf(fid,'\n');fprintf(fid,'\n');
fprintf(fid,'%c','-----------------------------------------------------------');

fclose(fid);

end






