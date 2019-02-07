clc
close all
clear all

maxIndex = 4;
for k=1:7
    figure(k)
    for i=1:maxIndex
        lunarAscent = load(strcat('/home/dominic/Software/tudatBundleTest/tudatBundle/tudatApplications/PropagationOptimizationAssignments/SimulationOutput/LunarAscent/dependentVariables',num2str(i-1),'.dat'));
        lunarAscentState = load(strcat('/home/dominic/Software/tudatBundleTest/tudatBundle/tudatApplications/PropagationOptimizationAssignments/SimulationOutput/LunarAscent/stateHistory',num2str(i-1),'.dat'));
        shapeParameters = load(strcat('/home/dominic/Software/tudatBundleTest/tudatBundle/tudatApplications/PropagationOptimizationAssignments/SimulationOutput/LunarAscent/thrustParameters',num2str(i-1),'.dat'));
        
        
        for j=1:3
            subplot(1,4,j)
            scatter(lunarAscent(:,1),lunarAscent(:,j+1),10*ones(size(lunarAscent(:,1))),shapeParameters(k)*ones(size(lunarAscent(:,1))))
            hold on
            if( i == 1 )
                xlabel('t [s]')
                if( j == 1 )
                    ylabel('Altitude [m]')
                elseif( j == 2)
                    ylabel('Speed w.r.t Moon CoM [m]')
                else
                    ylabel('Flight path angle [rad]')
                end
                grid on
            end
        end
        subplot(1,4,4)
        scatter(lunarAscentState(:,1),lunarAscentState(:,8),10*ones(size(lunarAscent(:,1))),shapeParameters(k)*ones(size(lunarAscent(:,1))))
        
        if( i==1)
            xlabel('t [s]')
            ylabel('Mass [kg]')
            grid on
        end
        
        if( i == maxIndex )
            colorbar
            if( k == 1 )
                suptitle( 'Thrust magnitude [N]' )
            elseif( k == 2 )
                suptitle( 'Node distance [s]' )
            elseif( k == 3 )
                suptitle( '$\theta_{1}$ [rad]' )
            elseif( k == 4 )
                suptitle( '$\theta_{2}$ [rad]' )
            elseif( k == 5 )
                suptitle( '$\theta_{3}$ [rad]' )
            elseif( k == 6 )
                suptitle( '$\theta_{4}$ [rad]' )
            elseif( k == 7 )
                suptitle( '$\theta_{5}$ [rad]' )
            end
        end
        
        hold on
    end
    
    set(figure(k), 'Units', 'normalized', 'Position', [0,0,1.0 1.0]);
    set(figure(k),'PaperUnits','centimeters','PaperPosition',[0 0 60 50]);
    set(figure(k),'PaperPositionMode','auto');
    
    saveas(figure(k),strcat('lunarAscentParameters',num2str(k)),'png');

end