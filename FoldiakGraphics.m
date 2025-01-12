classdef FoldiakGraphics
    methods (Static)
        function fig = showConnections(nX, nY, nOrient, nComplex, connections)
            SCALE_FACTOR = 5;
            fig = figure(Name="Connections to Complex Units");

            for iComplex = 1:nComplex
                subplot(2, 2, iComplex)
                hold on
            
                for orientation = 1:nOrient
                    U = zeros(nX * nY, 1);
                    V = zeros(nX * nY, 1);
            
                    for pos = 1:(nX * nY)
                        iConnection = (orientation - 1) * 64 + pos;
                        weight = connections(iConnection, iComplex);
            
                        % set vector based on orientation and synapse strength (connection
                        % weight)
                        switch orientation
                            case 1 %0˚
                                U(pos) = 1; V(pos) = 0;
                            case 2 %45˚
                                U(pos) = 1; V(pos) = 1;
                            case 3 %90˚
                                U(pos) = 0; V(pos) = 1;
                            case 4 %135˚
                                U(pos) = -1; V(pos) = 1;
                        end
                        % scale by connection strength and a constant scale
                        % factor (so arrows are visible in plot
                        U(pos) = U(pos) * weight * SCALE_FACTOR;
                        V(pos) = V(pos) * weight * SCALE_FACTOR;
            
                    end
            
                    X = repmat(1:8, 1, 8).';
                    Y = reshape(repmat(1:8, 8, 1), [], 1);
                    quiver(X, Y, U, V, 'off');
            
                end
            
                axis equal
                xlim([0.5 nX+0.5])
                ylim([0.5 nY+0.5])
                title(sprintf("Complex unit %1d", iComplex))
                hold off
            end
        end
    end
end