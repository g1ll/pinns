function [tmax,m,g] = dof2y( HdL,H0dL0,H2dL2,H1dL1,S1dH0,a,b,fi1,fi2,fi3,psi,nr,dirname,d,half,title)
% warning ('off','all');

if sum(isnan([ HdL H0dL0 H2dL2 H1dL1 S1dH0 a b ]))>=1
  %fprintf('NaN arguments:\n[ HdL H0dL0 H2dL2 H1dL1 S1dH0 a b ]\n1');
  throw(MException('Erro:pde','NaN arguments'));
end
 
if sum(isinf([ HdL H0dL0 H2dL2 H1dL1 S1dH0 a b ]))>=1
  %fprintf('Infinity arguments:\n[ HdL H0dL0 H2dL2 H1dL1 S1dH0 a b ]\n1');
  throw(MException('Erro:pde','Infinity arguments'));
end

s=0; %%save figures 1 to save
if d >=1
  dirname = strcat(dirname,'_',datestr(datetime('now'),'yyyymmdd_HHMMSS'));
end
%TESTAR GEOMETRIAS:
%dof2y(1,16.483148,0.067797,0.067797,0.32350359,0,0,0.015,0.015,0.04,1,1,'testGeo.png',1,0)

%Duplo-T
%dof2y(1,6,0.4,.4,0.5,0,0,0.015,0.015,0.04,1,1/10,'mesh-dt.png',2,0,'mesh-dt')

%Mesh Validate
%dof2y(1,20,0.14,14,0.45,0,0,0.02,0.02,0.02,1,1/10,'mesh.png',2,0,'mesh')

format long;
% Restricoes da cavidade
     %ficav = fi3 + 2*fi1 +2*fi2; % área total da cavidade H0*L0 + 2*phi1 + 2*phi2
       
    %Definindo os valores para geometria
    L = (1/HdL)^0.5;
    H = HdL*L;

    l0 = sqrt(fi3/H0dL0);
    h0 = H0dL0*l0;

    %l3 = sqrt(psi/H3dL3);
    %h3 = H3dL3*l3; 
    
    l1 = sqrt(fi1/(H1dL1));
    h1 = H1dL1*l1;

    l2 = sqrt(fi2/(H2dL2));
    h2 = H2dL2*l2;
    
    s1=S1dH0*h0;
    tmax = inf;
    try
      %if d ==1
        %fprintf('\nGeometria:\n\:');
        %[H,L,h0,l0,s1,h1,l1,h2,l2,a,b]
        %log_gls = fopen(strcat(pwd,'/','LOG_',dirname,'.txt'),'w');
      %end
      time =tic();
      %if dom == 2
      % tmax = pde_cav_2y_h(H,L,h0,l0,s1,h1,l1,h2,l2,nr,0);
      %end
      if d==1
        fprintf("\ndof2y:%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",[ HdL H0dL0 H2dL2 H1dL1 S1dH0 a b ]);
      end
      if half == 1
        tmax = pde_cav_2y_h(H,L,h0,l0,s1,h1,l1,h2,l2,a,b,nr,d,s,dirname,title);
      else
        tmax = pde_cav_2y(H,L,h0,l0,s1,h1,l1,h2,l2,a,b,nr,d,s,dirname,title);
      end
      tp = toc(time);
      %if d ==1
        %fprintf(log_gls,'Características do Problema:\n');
        %fprintf(log_gls,'HL: %.15f\n;h0l0: %.15f\n;h2l2: %.15f\n;h1l1: %.15f\n;s1h0: %.15f\n;a: %.3f\n;b: %.3f\n;fi1: %.15f\n;fi2: %.15f\n;fi3: %.15f\n;psi: %.15f\n;malha: %.15f\n;tmim: %.15f\n;time(sec): %.15f;',...
         % HdL, H0dL0,H2dL2,H1dL1,S1dH0,a,b,fi1,fi2,fi3,psi,1/nr,tmax,tp);
        %fclose(log_gls);
      %end
      %fprintf('Resultado: %.5f',tmax);
      %fprintf('\nTempo %.5f',tp);
      % fprintf('Resultado: %.5f',tmax);
      % fprintf('\nTempo %.5f',tp);
      % h0
      % l2
      % l0
      % psi
      % aux = 2*l2+l0
      % aux = aux * h0
      % H*L
      % aux/H*L
      if d==1
        fprintf("\t%0.15f",tmax);
      end
    catch error_dof
      if d>1
        fprintf("\ndof2y:%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%0.15f",[ HdL H0dL0 H2dL2 H1dL1 S1dH0 a b tmax]);
      end
      throw(MException('Erro:pde',error_dof.message));
    end
end




%TESTAR GEOMETRIAS:
% VALIDATION
% AppSci
% dof2y(1,20,0.14,14,0.45,0,0,0.02,0.02,0.02,1,100,'test1.png',1,1,'teste'): 0.0756 
% dof2y(1,6.5,0.16,1,0.45,0,0,0.00002,0.02,0.06,1,100,'test1.png',1,0,'teste'):0.0753

% TESTE MALHA

% dof2y(1,6,0.4,0.5,0.4,0,0,0.015,0.015,0.04,0.5,100,'test1.png',1,0,'teste')


% dof2y(1,16.483148,0,0,0.067797,0.32350359,0.067797,0.015,0.015,0.04,1,100,'test1.png',1,0,'teste')



%OPT GEO 7DOF
%dof2y(30,747.2638872125,-80,-58,23.1194,0.4,37.65696344365,0.015,0.015,0.04,1,100,'opt_geo_30-00',1,0,'')
%dof2y(10.333,254.2735,73,-80,2.2292,0.3,20.4322365141,0.015,0.015,0.04,1,100,'opt_geo_10-33',1,0,'')
%dof2y(3.778,18.13029806511,67,80,0.01242345789209,0.1,0.04669647280519,0.015,0.015,0.04,1,100,'opt_geo_3-77',1,0,'')
%dof2y(0.5,7.845086077602,1,-2,0.0343,0.3,0.0343,0.015,0.015,0.04,1,100,'opt_geo_0-50',1,0,'')
%dof2y(0.262,3.780127074406,1,0,0.01770685593803,0.3,0.0176,0.015,0.015,0.04,1,100,'opt_geo_0-26',1,0,'')
%dof2y(0.03,0.4318025403784,0,0,0.00201322015686,0.3,0.002012197658287,0.015,0.015,0.04,1,100,'opt_geo_0-03',1,0,'')

%OPT FIS
%fic=0.01;fi1=(fic*0.6)/4;fi2=fi1;fi3=fic*0.4;
%dof2y(0.03,1.2,1,4,0.00082740079914,0.7,0.000185099508115,fi1,fi2,fi3,1,100,'opt_geo_001',1,0,'')
%fic=0.05;fi1=(fic*0.6)/4;fi2=fi1;fi3=fic*0.4;
%dof2y(0.03,0.8436847577293,0,0,0.000953372937058,0.3,0.000950365296496,fi1,fi2,fi3,1,100,'opt_geo_005',1,0,'')
%fic=0.1;fi1=(fic*0.6)/4;fi2=fi1;fi3=fic*0.4;
%dof2y(0.03,0.4318025403784,0,0,0.00201322015686,0.3,0.002012197658287,fi1,fi2,fi3,1,100,'opt_geo_010',1,0,'')
%fic=0.15;fi1=(fic*0.6)/4;fi2=fi1;fi3=fic*0.4;
%dof2y(0.03,0.2954265589908,0,0,0.003188909237503,0.3,0.003178630955169,fi1,fi2,fi3,1,100,'opt_geo_015',1,0,'')
%fic=0.2;fi1=(fic*0.6)/4;fi2=fi1;fi3=fic*0.4;
%dof2y(0.03,0.2266949192431,0,0,0.004489782207918,0.3,0.00447666667088,fi1,fi2,fi3,1,100,'opt_geo_020',1,0,'')

%FI 0.01 OPT
%fic=0.01;fi1=(fic*0.6)/4;fi2=fi1;fi3=fic*0.4;
%dof2y(0.03,4.09,0,0,0.000183,0.3,0.000183,fi1,fi2,fi3,1,100,'opt_geo_001',1,0,'')

%FIS DIST FIC=0.1
%>>fic=0.1;fi3=fic*0.4;
%>>HdL=.03;L=(1/HdL)^0.5;H = HdL*L;
%>>h0=0.4318025403784*sqrt(fi3/0.4318025403784)

%FIS 0.8
%>>pfis=.8;fi1=(fic*pfis)/4;fi2=fi1;fi3=fic*(1-pfis);
%>>H0dL0=h0/(fi3/h0);
%>>disp([fic fi1 fi2 fi3 H0dL0]');
%>>dof_2y_opt(0.03,H0dL0,0.3,0,0,fi1,fi2,fi3,1,100,'opt_geo_01080',1,0,'') 

%FIS 0.9
%>>pfis=.9;fi1=(fic*pfis)/4;fi2=fi1;fi3=fic*(1-pfis);
%>>H0dL0=h0/(fi3/h0);
%>>disp([fic fi1 fi2 fi3 H0dL0]');
%>>dof_2y_opt(0.03,H0dL0,0.3,0,0,fi1,fi2,fi3,1,100,'opt_geo_01090',1,0,'') 


%FIS 0.95
%>>pfis=.95;fi1=(fic*pfis)/4;fi2=fi1;fi3=fic*(1-pfis);
%>>H0dL0=h0/(fi3/h0);
%>>disp([fic fi1 fi2 fi3 H0dL0]');
%>>dof_2y_opt(0.03,H0dL0,0.3,0,0,fi1,fi2,fi3,1,100,'opt_geo_01095',1,0,'') 

%FIS DIST FIC=0.2
%FIS 0.8
% >> fic=0.2;fi1=(fic*0.8)/4;fi2=fi1;fi3=fic*0.2;                                   
% >> dof_2y_opt(0.03,0.4566949192431,0.3,0,0,fi1,fi2,fi3,1,100,'opt_geo_020',1,0,'')                              

% ans =

%    0.040000000000000
%    0.040000000000000
%    0.040000000000000


% P2 =

%      0


% ans =

%      1     1     1


% RESULTADO:
% Fic: 0.20000
% HL: 0.03000
% H0L0: 0.45669
% beta: 0.00000
% H2L2: 0.005371850921554
% S1H0: 0.30000
% a: 0.00000
% H1L1: 0.005371850921554
% Tmax: 0.000759012394308
% Tempo 1.23492

% ans =

%      7.590123943075660e-04


%FIS 0.9
%fic=0.2;fi1=(fic*0.9)/4;fi2=fi1;fi3=fic*0.1;
% dof_2y_opt(0.03,0.9566949192431,0.30,0,0,fi1,fi2,fi3,1,100,'opt_geo_020',1,0,'')

% ans =

%    0.045000000000000
%    0.045000000000000
%    0.020000000000000


% P2 =

%      0


% ans =

%      1     1     1


% RESULTADO:
% Fic: 0.20000
% HL: 0.03000
% H0L0: 0.95669
% beta: 0.00000
% H2L2: 0.005721562286617
% S1H0: 0.30000
% a: 0.00000
% H1L1: 0.005721562286617
% Tmax: 0.000696500957096
% Tempo 46.11842

% ans =

%      6.965009570958525e-04

%FIS 0.95
%fic=0.2;fi1=(fic*0.95)/4;fi2=fi1;fi3=fic-2*fi1-2*fi2;
%[fi3 fi1 fi2]'
% dof_2y_opt(0.03,1.9,0.3,0,0,fi1,fi2,fi3,1,100,'opt_geo_020',1,0,'')

% fi1 =

%    0.047500000000000


% fi2 =

%    0.047500000000000


% fi3 =

%    0.010000000000000


% P2 =

%      0


% a =

%      0


% H1L1min =

%    0.005887228501047


% ans =

%      1     1     1


% RESULTADO:
% Fic: 0.20000
% Fi1: 0.04750
% Fi2: 0.04750
% Fi3: 0.01000
% HL: 0.03000
% H0L0: 1.90000
% beta: 0.00000
% H2L2: 0.005887228501047
% S1H0: 0.30000
% a: 0.00000
% H1L1: 0.005887228501047
% Tmax: 0.000671630099888
% Tempo 70.95865

% ans =

%      6.716300998875185e-04