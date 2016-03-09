function morphGui

% Copyright 2016 David Ochsenbein
% 
% This file is part of CHIMMP.
% 
% CHIMMP is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% 
% CHIMMP is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% close all
global hp  Ival crysnames hpu milstr

hf = figure;

set(hf,'position',[300 150 1200 700],'color','w','menubar','none')
ha1 = axes('position',[0.075 0.175 0.5 0.8]);

ha2 = axes('position',[0.6 0.1 0.35 0.4]);

hslider = [];
try
    rvec = evalin('base','r1vec');
    U = evalin('base','U');
    crystal = evalin('base','crystal');
    Cmat = evalin('base','Cmat');
catch
    updateMap
end

% Crystals = load('CrystalData');
% Crystals = Crystals.Crystals;
% Mc = struct2cell(Crystals);
% Mc = Mc(6,:);
% I = false(length(Mc),1);
% for i = 1:length(Mc)
%     if length(Mc{i}(1,:))~=3
%         I(i) = false;
%     else
%         I(i) = true;
%     end
% end
% Crystals = Crystals(I,:);
% 
% for i = 1:length(Crystals)
 crysnames= {'blub'};
% end

% 
% [~,Ival] = getCrystalInfo(crystal.heteronyms{1});
% Ival = find(Ival);

hpu = uicontrol(hf,'style','popupmenu','units','normalized','position',[0.65    0.6    0.2500    0.2000],'string',crysnames,'value',Ival,...
    'fontsize',14,'callback',@popupcallback,'backgroundcolor','w');

    function posChange(pos)

    axes(ha2);
    cla
    
    if ndims(U) == 2
        L = round([1 exp(pos(1)) exp(pos(2))]'*1e3)/1e3;
    elseif ndims(U) == 3
        sliderVal = get(hslider,'value');
        sliceI = round(sliderVal*(length(rvec{3})-1))+1;
        L = round([1 exp(pos(1)) exp(pos(2)) rvec{3}(sliceI)*exp(pos(2))]'*1e3)/1e3;
    end

    Pbase = Polyhedron(crystal.H,ones(size(crystal.H,1),1));
    P = Polyhedron(crystal.H,crystal.M*L);

    plot(P)
    axis equal
    axis off
    hold on
    
    [pFace] = colorCrystal(Pbase,P,crystal.Mfull,[]);
    
    AA = get(findall(gcf,'type','patch'),'vertices');
    BB = get(findall(gcf,'type','patch'),'faces');
    if length(AA(:,1)) == length(BB(:,1))
        set(findall(gcf,'type','patch'),'vertices',[AA; 0 0 0]);
    end
    
    set(findall(hf,'type','patch'),'facevertexcdata',pFace,'facecolor','flat','facealpha',1)

    
    htit = title(['h = [1, ',num2str(L(2),'%4.2f'),', ',num2str(L(3),'%4.2f'),', ',num2str(L(4),'%4.2f'),']'],'fontsize',14);
    
    set(htit,'units','normalized','position',[0.5 1 0])

    end

    function popupcallback(hObject,~,~)
    
    
    updateMap(crysnames{get(hObject,'value')});
    
    [~,Ival] = getCrystalInfo(crystal.heteronyms);
    end


    function updateMap(crysname)

        if nargin<1
            
            
%             crystal = getCrystalInfo('blglu');
%             crystal.M = [1,0,0,0,0,0;1,0,0,0,0,0;1,0,0,0,0,0;1,0,0,0,0,0;0,1,0,0,0,0;0,1,0,0,0,0;0,0,0,1,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,1,0,0,0;0,0,0,0,0,1;0,0,0,0,0,1;0,0,0,0,0,1;0,0,0,0,0,1;0,0,0,0,1,0;0,0,0,0,1,0];
            
            crystal.H = [-9.30702956532859e-17,-1.04077766150235e-16,1.00000000000000;9.30702956532859e-17,1.04077766150235e-16,-1.00000000000000;-0.760447467379201,0.614567348899053,-0.209825220180388;0.760447467379201,-0.614567348899053,0.209825220180388;-0.829652744433316,-0.486102158099949,0.274556033524742;0.829652744433316,0.486102158099949,-0.274556033524742;-0.604912610090194,0.488869455216894,0.628559774333202;0.604912610090194,-0.488869455216894,-0.628559774333202;0.148067858804268,0.916206732216432,0.372345448515765;-0.148067858804268,-0.916206732216432,-0.372345448515765];
            crystal.M = [1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;0,1,0,0;0,1,0,0;0,0,0,1;0,0,0,1;0,0,1,0;0,0,1,0];
            crystal.Mfull = crystal.M;
            
            [rvec,U,Cmat] = mapDrawer(crystal.H,crystal.M);
            
%             crystal.M = [1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0;0,1,0,0;0,1,0,0;0,0,0,1;0,0,0,1;0,0,1,0;0,0,1,0;0,0,0,3.16426789911892;0,0,0,3.16426789911892;0,0,0,3.16426789911892;0,0,0,3.16426789911892;0,0,0,2.38569814452229;0,0,0,2.38569814452229];
%             crystal.H = [-9.30702956532859e-17,-1.04077766150235e-16,1.00000000000000;9.30702956532859e-17,1.04077766150235e-16,-1.00000000000000;-0.760447467379201,0.614567348899053,-0.209825220180388;0.760447467379201,-0.614567348899053,0.209825220180388;-0.829652744433316,-0.486102158099949,0.274556033524742;0.829652744433316,0.486102158099949,-0.274556033524742;-0.604912610090194,0.488869455216894,0.628559774333202;0.604912610090194,-0.488869455216894,-0.628559774333202;0.148067858804268,0.916206732216432,0.372345448515765;-0.148067858804268,-0.916206732216432,-0.372345448515765;-0.411491282762496,-0.624112353132361,0.664197782951782;0.411491282762496,0.624112353132361,-0.664197782951782;0.960981254230386,0.180033558083600,0.210007016505529;-0.960981254230386,-0.180033558083600,-0.210007016505529;0.488869455216894,-0.395087389909806,0.777761280914777;-0.488869455216894,0.395087389909806,-0.777761280914777];
%             crystal.Mfull = [1,0,0,0,0,0;1,0,0,0,0,0;1,0,0,0,0,0;1,0,0,0,0,0;0,1,0,0,0,0;0,1,0,0,0,0;0,0,0,1,0,0;0,0,0,1,0,0;0,0,1,0,0,0;0,0,1,0,0,0;0,0,0,0,0,1;0,0,0,0,0,1;0,0,0,0,0,1;0,0,0,0,0,1;0,0,0,0,1,0;0,0,0,0,1,0];
            
            if ndims(U)>2
                hslider = uicontrol(hf,'style','slider','units','normalized','position',[0.7 0.7 0.2 0.03],'min',0,'max',1,'value',0.5,'sliderstep',[0.01 0.05],'callback',@sliderCallback);
                try
                    addlistener(hslider,'ContinuousValueChange',@ sliderCallback);
                catch
                    handle.addlistener(hslider,'ActionEvent',@ sliderCallback);
                end

            end
        else
            crystal = getCrystalInfo(crysname);
            [rvec,U,Cmat] = mapDrawer(crystal.H,crystal.M);
        end
        
        axes(ha1)
        h = imagesc(log(rvec{1}),log(rvec{2}),U(:,:,51));
        set(h,'cdata',squeeze(Cmat(:,:,51,:)))
        set(gca,'yDir','normal','fontsize',14)

        grid on
        
        hp = impoint(ha1,0,0);

        setColor(hp,'k')

        id = addNewPositionCallback(hp,@posChange);
        
        axes(ha2)
        hold on
        posChange([0 0])
        
        axes(ha1)
%         xlabel(['$$\ln \left(L_{',milstr{1},'}/L_{',milstr{3},'}\right)$$'],'interpreter','latex','fontname','calibri','fontsize',20,'interpreter','latex')
%         ylabel(['$$\ln \left(L_{',milstr{2},'}/L_{',milstr{3},'}\right)$$'],'interpreter','latex','fontname','calibri','fontsize',20,'interpreter','latex')
        

        
    end

    function sliderCallback(hObject,eventdata)
        curPos = getPosition(hp);
        axes(ha1)
        
        cla(ha1)
        
        sliderVal = get(hObject,'value');
        sliceI = round(sliderVal*(length(rvec{3})-1))+1;
        h = imagesc(log(rvec{1}),log(rvec{2}),U(:,:,sliceI));
        set(h,'cdata',squeeze(Cmat(:,:,sliceI,:)))
        set(gca,'yDir','normal','fontsize',14)
        hp = impoint(ha1,curPos);
        id = addNewPositionCallback(hp,@posChange);
        posChange(curPos);
    end

    
end


