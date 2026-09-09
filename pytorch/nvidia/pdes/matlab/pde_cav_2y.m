%Objective Function Double-Y Shaped Cavity
%%	H = Solid Heigth
%%  L = Solid Width
%%	h0 = Cavity Trunk heigth
%%	l0 = Cavity Trunk length/width
%%	....
%%	a = alpha angle
%%	b = beta angle
%%	nr = size of triangle FEM
%%	d = display plot graph
%%	s = save figures
%%  dirname = file names to save
function [saida] = pde_cav_2y(H,L,h0,l0,s1,h1,l1,h2,l2,a,b,nr,d,s,dirname,txt)
	% warning ('off','all');
	
	saida = inf;
	model = inf;
	g=0;
	%termalmodel = createpde('thermal','steadystate');
	model = createpde(1);

	%Geometria:
	%coordenadas de cada ponto (18 pontos) Desenho horário cavidade y
	x1 = -l0/2;	 				y1 = 0;
	x2 = x1;	 				y2 = s1-(h1/sind(90-a))/2;
	x3 = -(cosd(a)*l1+l0/2);	y3 = l1*sind(a)+y2;
	x4 = x3;	 				y4 = y3+h1/sind(90-a);
	x5 = x1;	 				y5 = y2+h1/sind(90-a);
	x6 = x1;	 				y6 = h0 - h2/sind(90-b);
	x7 = -(cosd(b)*l2+l0/2);	y7 = (h0+l2*sind(b))- h2/sind(90-b);
	x8 = x7;					y8 = h0+l2*sind(b);
	x9 = x1;					y9 = h0;
	x10 = -x1;	 				y10 = h0;
	x11 = -x8;					y11 = y8;
	x12 = x11;	 				y12 = y7;
	x13 = x10;					y13 = y6;
	x14 = x10;	 				y14 = y5;
	x15 = -x3;	 				y15 = y4;
	x16 = x15;					y16 = y3;
	x17 = x10;	 				y17 = y2;
	x18 = x10;	 				y18 = 0;
	
	%coordenados do sólido
	xr1 = -L/2;			yr1 = 0;
	xr2 = xr1;			yr2 = H;
	xr3 = L/2;			yr3 = yr2;
	xr4 = xr3;			yr4 = 0;

	P1 =  [2 4 xr1 xr2 xr3 xr4,...
			   yr1 yr2 yr3 yr4 ]';

	P2 =  [2 18 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18,...
				y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18]';

	
	%validade geometry
	px = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 xr1 xr2 xr3 xr4]';
	py = [y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 yr1 yr2 yr3 yr4]';

    polysh = [px py];
	regions = polyshape(polysh).NumRegions;
	% regions
	% plot(polyshape(polysh));
	 % For a polygon, gstat = 0 indicates that the polygon is closed and does not 
	 % intersect itself, i.e., it has a well-defined, unique interior region.
	 % 1 indicates an open and non-self-intersecting polygon,
	 % 2 indicates a closed and self-intersecting polygon,
	 % and 3 indicates an open and self-intersecting polygon.
	 P1 = [P1;zeros(length(P2)-length(P1),1)];
	 gstat = csgchk([P2,P1],L/2,H);
	%  [gstat regions]
	 if sum([gstat regions]==1) <= 3
		 ns = (char('P2','P1'))';
		 sf = 'P1-P2';
		 g = decsg([P2,P1],sf,ns);
	 else
		 errormsg = sprintf('Erro na geometria! gstat %d %d regions %d\n', gstat,regions);
		 throw(MException('Erro:geo',errormsg));
	 end
	ns = (char('P2','P1'))';
	sf = 'P1-P2';
	g = decsg([P2,P1],sf,ns);
	geometryFromEdges(model,g);
	%dlmwrite(strcat(pwd,'/_log_descg.csv'),g,'precision','%.15f','-append');
	if d == 2
		% f1  = figure(1);
		% axis([-0.6 0.6 0 1.1]);
		% title 'Block Geometry With Edge Labels Displayed'
		% pdegplot(model,'EdgeLabels','on')
		% set(gca,'TickLength',[0 0])	
		% caxis manual
		 ax = gca;
	     ax.DataAspectRatio = [1 1 1];
		 ax.DataAspectRatioMode = 'manual';
		% pbaspect([1 1 1])
		% daspect([1 0.5 2])
		% axis manual
		% ax = gca;
	      %ax.DataAspectRatio = [0.01 1 1];
		% ax.DataAspectRatioMode = 'manual';
		
		if s==2
			try
				saveas(f1,strcat('GEO_',dirname,'.png'));
			catch
				fprintf('\nErro save figure! geo\n');
			end
		end		
	end

	%%Condições de Contorno
	%Descobrir quais seguimentos estao nas laterais
	sg = length(g);
	s_edge = [];
	% sg
	%restrições geométricas
	if sg ~= 22
		throw(MException('Erro:geo','Erro na geometria! sg ~= 22\n'));
		if d == 2
			f1  = figure(1);
			axis([-0.6 0.6 0 1.1]);
			title 'Block Geometry With Edge Labels Displayed'
			pdegplot(model,'EdgeLabels','on');
		end
	end
	et = 1e-15; %erro truncamento
	%procurando as laterais
	for i=1:1:sg
		c = g(:,i);
		if or((abs(c(3)+l0/2) < et && abs(c(2)+L/2) < et),(abs(c(2)+l0/2) < et && abs(c(3)+L/2) < et)) %seg inf negativo
		 	s_edge = [s_edge i];
		elseif abs(c(2) + L/2) < et && abs(c(3)+L/2)< et %seg lat negativo
			s_edge = [s_edge i];
		% elseif or((c(2) == -L/2 && c(3) == L/2),(c(3) == -L/2 && c(2) == L/2)) %seg superior
		elseif or(abs(c(2) - L/2) < et && abs(c(3)+L/2) < et,abs(c(2) +L/2) < et && abs(c(3)-L/2) < et)  %seg superior
			s_edge = [s_edge i];
		% elseif c(2) == L/2 && c(3) == L/2 %seg lateral positivo
		 elseif abs(c(2) - L/2) < et && abs(c(3)-L/2) < et %seg lateral positivo
			s_edge = [s_edge i];
		elseif or((abs(c(3)-L/2) < et && abs(c(2)-l0/2) < et ),(abs(c(2)-L/2) < et && abs(c(3)-l0/2) < et)) %seg inferior positivo
			s_edge = [s_edge i];
		end
		% 'c'
		% [c(2) c(3)]
		% 'L'
		% L/2
		% 'd'
		% c(2)+L/2
 	end
	%TESTAR ERRO DE GEOMETRIA
	%dof2y(0.5,12.148946609406723,0.033076793491827,0.033076793491827,0.9,0,0,0.015,0.015,0.04,1,10,'ex_5dof_cav2y',1)
	%  s_edge
%  throw(MException('Erro:geo','TESTE\n'));
	if length(s_edge) ~= 5
		%sg
		%g
		% s_edge
		throw(MException('Erro:geo','Erro na geometria! s_edge ~= 5\n'));
		if d == 2
			f1  = figure(1);
			axis([-0.6 0.6 0 1.1]);
			title 'Block Geometry With Edge Labels Displayed'
			pdegplot(model,'EdgeLabels','on');
		end
	end
	%linhas da cavidade
	c_edge  = [];

	for i=1:sg
		if i ~= s_edge(1) && i ~=  s_edge(2) && i ~=  s_edge(3) && i ~= s_edge(4) &&  i ~= s_edge(5)
			c_edge = [c_edge i];
		end
	end
	%mostrando linhas do solido
	%s_edge
	%mostrando linhas da cavidade
	%c_edge
	
	%aplicando condições de contorno aos vertices
	applyBoundaryCondition(model,'neumann','Edge',s_edge,'q',0,'g',0);
	applyBoundaryCondition(model,'dirichlet','Edge',c_edge,'h',1,'r',0);
	% applyBoundaryCondition(model,'neumann','Edge',c_edge,'q',100,'g',0);//equivalente
	% applyBoundaryCondition(model,'neumann','Edge',c_edge,'q',0.05,'g',0)

	specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);

	mesh = generateMesh(model,'Hmax',1/nr);
	% mesh
	if d==2
		% f2  = figure(2);
		% pdeplot(model);
		if s==2
			try	
				saveas(f2,strcat('MESH_',dirname,'.png'));
			catch
				fprintf('\nErro save figure!\n');
			end
		end
	end

	results = solvepde(model);
	u = results.NodalSolution;
	saida = max(u);
	% fprintf('\nNós: %d\n',length(u));
	% fprintf('Elementos: %d\n',length(results.Mesh.Elements));
	
	if d>=1
		%set(gcf,'Renderer','Painters');
		figure('visible','off');
		f3  = figure(3);
		axis([-0.6 0.6 0 1.1]);
		title 'GEOMETRIA ÓTIMA E CAMPO DE TEMPERATURAS'
	%	com gradiente ou fluxo de calor
		
		ux = results.XGradients;
		uy = results.YGradients;
		%pdeplot(model,'XYData',u,'colormap','hot','contour','on','FlowData',[ux,uy]);	
		pdeplot(model,'XYData',u,'colormap','jet','contour','off');
		% caxis([32.5 34.5]);
		caxis([0 0.1]);
		ax = gca;
	    % ax.DataAspectRatio = [5 1 1];
		ax.DataAspectRatioMode = 'manual';
		ax.FontSize = 12;
		% ax.XAxis.Color = 'w';
		ax.XColor = 'none';
		ax.YColor = 'none';
		if strcmp(txt,'')~=1
			title(strcat(txt,sprintf(" Tmax: %.6f",saida)));
		end
		set(gca,'TickLength',[0 0]); 
		caxis manual
		if s==1
			try
				% savefig(f1,strcat(dirname,'_TEMP'));
				print(f3,strcat(dirname,'_TEMP'),'-dpng','-r200') %resolução screen png
				
				%saveas(f3,strcat(dirname,'_TEMP','.jpg'));
				%close
			catch erro
				fprintf('\nErro save figure TEMP!\n');
				fprintf('\nErro: %s',erro.message);
				%close
			end
		end
	end
end