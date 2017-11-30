function PlotLatticeState(LattState)

Nx = size(LattState,1);
Ny = size(LattState,2);


iOccup = [];
jOccup = [];
iUnocc = [];
jUnocc = [];
for i = 1:Nx
    for j = 1:Ny

%         % 4-fold lattice
%         if LattState(i,j) == 1;
%             iOccup = [iOccup i];
%             jOccup = [jOccup j];
%         else
%             iUnocc = [iUnocc i];
%             jUnocc = [jUnocc j];
%         end
        
        % 6-fold lattice
        if LattState(i,j) == 1;
            iOccup = [iOccup i-mod(j,2)*1/2];
            jOccup = [jOccup j];
        else
            iUnocc = [iUnocc i-mod(j,2)*1/2];
            jUnocc = [jUnocc j];
        end
        
    end
end


MSiz = 28;
clf
hold on
plot(iOccup,jOccup,'o','MarkerFaceColor',[0.4 0 1]*1.0,'Color',[0.4 0 1]*1.0,'MarkerSize',MSiz)
plot(iUnocc,jUnocc,'o','MarkerFaceColor',[1 1 1]*0.8,'Color',[1 1 1]*0.8,'MarkerSize',MSiz)
for i = 1:Nx
    for j = 1:Ny
        text(i-mod(j,2)*1/2,j,[num2str(i) ',' num2str(j)],'HorizontalAlignment','Center','Color','w')
    end
end
hold off

axis equal
box on
set(gcf,'Color','w')

return