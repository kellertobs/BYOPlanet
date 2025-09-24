% plot and print output

% set plotting scales
SCL = 200;
scl = [1/20,2];


% Figure 1: top-down disk view
figure(1); clf; set(gcf,'Visible','on','PaperUnits','centimeters','PaperSize',[24,18]);
rectangle('Position',[-1 -1  2  2],'Curvature',[1,1],'EdgeColor','k','LineStyle',':','LineWidth',2); hold on; box on; axis equal;
rectangle('Position',[-5 -5 10 10],'Curvature',[1,1],'EdgeColor','k','LineStyle',':','LineWidth',2);
scatter(X(3:N,1),X(3:N,2),Rtot(3:N).^2.*SCL,[0,0,1],'filled');
scatter(X(3:N,1),X(3:N,2),Rrck(3:N).^2.*SCL,[0,1,0],'filled');
scatter(X(3:N,1),X(3:N,2),Rmtl(3:N).^2.*SCL,[1,0,0],'filled');
scatter(X(1:2,1),X(1:2,2),[Rsun,Rggt].^2.*scl,zeros(2,3),'filled');

[~,I] = sort(M,'descend');
[~,J] = sort(r(I(3:7)),'ascend');
for n = 1:5
    ind = I(J(n)+2);
    scatter(X(ind,1),X(ind,2),(Rtot(ind)+0.01).^2.*SCL,[1,0,1],'LineWidth',1.5)
end

scatter(-9,-8,Rear^2.*SCL,[0,0,0]);
scatter(-9,-9,Rggt^2.*scl(2),[0,0,0],'filled');
text(-8,-8,'Earth','FontSize',13);
text(-8,-9,'star, gas giant not to scale','FontSize',13);

text(0.025,0.95,[num2str(time/yr,4),' years'],'Units','normalized','FontSize',13);
text(0.025,0.90,[int2str(N  ),' bodies'    ],'Units','normalized','FontSize',13);
text(0.025,0.85,[int2str(CLS),' collisions'],'Units','normalized','FontSize',13);

axis([X(1,1)-10,X(1,1)+10,X(1,2)-10,X(1,2)+10]);
title('Accretionary Disk – Top View','FontSize',15);
xlabel('Distance [AU]','FontSize',15)
ylabel('Distance [AU]','FontSize',15)

drawnow;
print(gcf,['./out/',runID,'/',runID,'_top_',int2str(k/nop)],'-dpdf','-vector','-fillpage');


% Figure 2: radial disk view
figure(2); clf; set(gcf,'Visible','on','PaperUnits','centimeters','PaperSize',[24,18]);
line([1,1],[-1,1],'Color','k','LineStyle',':','LineWidth',2); hold on; box on; axis equal;
line([5,5],[-1,1],'Color','k','LineStyle',':','LineWidth',2); hold on; box on; axis equal;
scatter(r(3:N),X(3:N,3),Rtot(3:N).^2.*SCL,[0,0,1],'filled'); hold on; box on; axis equal;
scatter(r(3:N),X(3:N,3),Rrck(3:N).^2.*SCL,[0,1,0],'filled');
scatter(r(3:N),X(3:N,3),Rmtl(3:N).^2.*SCL,[1,0,0],'filled');
scatter(r(1:2),X(1:2,3),[Rsun,Rggt].^2.*scl,zeros(2,3),'filled');

ra = 1;
for n = 1:5
    ind = I(J(n)+2);
    scatter(r(ind),X(ind,3),(Rtot(ind)+0.01).^2.*SCL,[1,0,1],'LineWidth',1.5)
    a  = 0.25; r0 = 2.5;
    ra = max(mean(r(I(J(1:5)+2)))+(r(ind)-mean(r(I(J(1:5)+2))))./1.2,ra+1);
    Xa = 2.5+X(ind,3).*2;
    line([r(ind),ra],[X(ind,3),Xa],'Color','k','LineWidth',1.5);
    scatter(ra,Xa,(Rtot(ind)*3).^2.*SCL,[0,0,1],'filled'); hold on; box on; axis equal;
    scatter(ra,Xa,(Rrck(ind)*3).^2.*SCL,[0,1,0],'filled');
    scatter(ra,Xa,(Rmtl(ind)*3).^2.*SCL,[1,0,0],'filled');
end
text(1.35,4,'most massive bodies (3x)','FontSize',13);

scatter(1.0,-1.3,Rear^2.*SCL,[0,0,0]);
scatter(1.0,-2.0,Rggt^2.*scl(2),[0,0,0],'filled');
text(2,-1.3,'Earth','FontSize',13);
text(2,-2.0,'star, gas giant not to scale','FontSize',13);

text(0.5,6.0,[num2str(time/yr,4),' years'],'FontSize',13);
text(0.5,5.4,[int2str(N  ),' bodies'      ],'FontSize',13);
text(0.5,4.8,[int2str(CLS),' collisions'  ],'FontSize',13);

plot(5,6.0,'r.','MarkerSize',20);
plot(5,5.4,'g.','MarkerSize',20);
plot(5,4.8,'b.','MarkerSize',20);

text(5.5,6.0,'metal','FontSize',13);
text(5.5,5.4,'rock' ,'FontSize',13);
text(5.5,4.8,'ice'  ,'FontSize',13);

axis([0,X(1,1)+20,X(1,1)-2.5,X(1,1)+6.5]);
title('Accretionary Disk – Radial View','FontSize',15);
xlabel('Radial Distance [AU]','FontSize',15)

drawnow;
print(gcf,['./out/',runID,'/',runID,'_rad_',int2str(k/nop)],'-dpdf','-vector','-fillpage');


% Figure 3: Planetary Zoo
figure(3); clf; set(gcf,'Visible','on','PaperUnits','centimeters','PaperSize',[24,18]);
scatter(log10(M(3:N)),log10(Rtot(3:N)),SCL/2,r(3:N),'filled'); hold on; box on;
cb = colorbar('FontSize',13); cb.Label.String = 'Distance from Star [AU]'; colormap(flipud(copper)); clim([0,15]);

for n = 1:5
    ind = I(J(n)+2);
    scatter(log10(M(ind)),log10(Rtot(ind)),(1+0.01).*SCL/2,[1,0,1],'LineWidth',1.5)
end

scatter(0,log10((1/1.5).^(1/3)),SCL/2,[0,0,0]);
text(0.1,-0.1,'Earth','FontSize',13);

axis([-3,1,-1.5,1.0]);
title('Planetary Zoo','FontSize',15);
xlabel('log_{10} Planetary Mass [Earth Mass]','FontSize',15)
ylabel('log_{10} Planetary Radius [Earth Radius]','FontSize',15)

text(-2.9,0.9,[num2str(time/yr,4),' years'],'FontSize',13);
text(-2.9,0.75,[int2str(N  ),' bodies'      ],'FontSize',13);
text(-2.9,0.6,[int2str(CLS),' collisions'  ],'FontSize',13);

drawnow;

print(gcf,['./out/',runID,'/',runID,'_zoo_',int2str(k/nop)],'-dpdf','-vector','-fillpage');
