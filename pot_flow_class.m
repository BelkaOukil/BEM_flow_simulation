
classdef pot_flow_class
    
    properties
        U;
        H;
        N;
        Ng;
        r;
        th;
        sigma;
        %         thp;
        X;
        Y;
        nb_iso;
        interp;
    end
       
    methods
               
%%%%%%%%%%%%% Constructeur
        function obj = pot_flow_class(valH,valN,valNg,valR,valU,val_nb_iso)
            obj.H = valH;
            obj.U = valU;
            obj.N = valN;
            obj.Ng = valNg;
            obj.r = valR;
            obj.nb_iso = val_nb_iso;
            obj.th=0:pi/(valN-1):pi;
            obj.sigma=sin(obj.th)/pi;
            obj.X=linspace(-10,10,valNg);
            obj.Y=linspace(-10,10,valNg);
        end
%%%%%%%%%%%%%
               
%%%%%%%%%%%%% Noyeau
        function K = Kernel(obj,x,y)           
            K = log(sqrt((x+obj.r*cos(obj.th)).^2+(y-obj.r*sin(obj.th)).^2))...
                +...
                log(sqrt((x-2*obj.H+obj.r*cos(obj.th)).^2+(y-obj.r*sin(obj.th)).^2));           
        end
%%%%%%%%%%%%%
        
        
%%%%%%%%%%%%% Gradient du noyeau
        function K_d = Kernel_d(obj,x,y)            
            C=cos(obj.th);
            S=sin(obj.th);
            R = (sqrt(x^2+y^2));            
            K_d = ...
                (x.*C + y.*S) - obj.r*(C.*x/R + S.*y/R)...
                ./...
                ((x - obj.r*x/R).^2 + (y - obj.r*y/R).^2)...
                +...
                (x*C + y*S - obj.r*(C.*x/R + S.*y/R) - 2*C*obj.H)...
                ./...
                ((x - 2*obj.H - obj.r*x/R).^2 + (y - obj.r*y/R).^2);            
        end
%%%%%%%%%%%%%

%%%%%%%%%%%%% Calcul des sources
        function obj = calcul_source(obj)            
            RC = obj.r*cos(obj.th);
            RS = obj.r*sin(obj.th);
            RC = RC(2:obj.N-1);
            RS = RS(2:obj.N-1);
            Mat_A = pi*eye(obj.N-2);
            for i=1:obj.N-2
                for j=1:obj.N-2                   
                    Mat_A(i,j)= Mat_A(i,j)+...
                        pi/obj.N*...
                        (RC(i)*(RC(i)+RC(j)-2*obj.H)+RS(i)*(RS(i)-RS(j)))/((RC(i)+RC(j)-2*obj.H)^2+...
                        (RS(i)-RS(j))^2)...
                        -pi/obj.N*...
                        (RC(i)*(RC(i)+RC(j)-2*obj.H)+RS(i)*(RS(i)+RS(j)))/((RC(i)+RC(j)-2*obj.H)^2+...
                        (RS(i)+RS(j))^2);                    
                end
            end         
            obj.sigma(2:obj.N-1)=linsolve(Mat_A,(1/pi)*RS');
        end
%%%%%%%%%%%%%
       
%%%%%%%%%%%%% Intégrale de SIGMA
        function integ_sigma(obj)
            format long
            S = trapz(obj.th,obj.sigma.*obj.sigma);
            disp(S)
        end
%%%%%%%%%%%%%
               
%%%%%%%%%%%%% Calcul du POTENTIEL
        function P = calcul_potentiel(obj,x,y)
            RC = obj.r*cos(obj.th);
            RS = obj.r*sin(obj.th);
            K=zeros(1,obj.N);
            for j=1:obj.N
                K(1,j) = 0.5*log((x-2*obj.H+RC(j)).^2+(y-RS(j)).^2)-0.5*log((x-2*obj.H+RC(j)).^2+...
                         (y+RS(j)).^2)...                
                          +...
                         0.5*log((x-RC(j)).^2+(y-RS(j)).^2)-0.5*log((x-RC(j)).^2+(y+RS(j)).^2);
            end
                P = -trapz(obj.th,K.*sin(obj.th));
            obj.sigma
        end
%%%%%%%%%%%%%

%%%%%%%%%%%%% affichage POTENTIEL
        function affiche_phi(obj)
            phi = zeros(obj.Ng,obj.Ng);
            for i=1:obj.Ng
                for j=1:obj.Ng
                    
                    if  sqrt(obj.X(i).^2+obj.Y(j).^2)<obj.r*1.001                            
                        phi(j,i)=1/0;                       
                    elseif obj.X(i)>obj.H      
                        phi(j,i)=1/0;      
                    else      
                        phi(j,i) =  obj.U*obj.Y(j) + obj.calcul_potentiel(obj.X(i),obj.Y(j));
                    end
                end   
            end      
            figure(1)
            contour(obj.X,obj.Y,phi,obj.nb_iso);
            pbaspect([1 1 1])
            hold on
            line([obj.H obj.H],[-10 10])
            grid on
            viscircles([0 0],obj.r,'color','k');
        end
%%%%%%%%%%%%%
                
%%%%%%%%%%%%% Caclul de la fnct de courant
        function C = calcul_courant(obj,x,y)
            RC = obj.r*cos(obj.th);
            RS = obj.r*sin(obj.th);
            K_courant=zeros(1,obj.N);
            for j=1:obj.N
                K_courant(1,j) = new_atan2((y-RS(j)),(x-2*obj.H-RC(j)))-new_atan2((y+RS(j)),...
                                 (x-2*obj.H-RC(j)))...
                                 +...                 
                                 new_atan2((y-RS(j)),(x-RC(j)))-new_atan2((y+RS(j)),(x-RC(j)));         
            end
                C = -trapz(obj.th,K_courant.*obj.sigma);
        end
%%%%%%%%%%%%%

%%%%%%%%%%%%% Affichage de la fnct de courant
        function affiche_psy(obj)
            psy = zeros(obj.Ng,obj.Ng);
            for i=1:obj.Ng
                for j=1:obj.Ng                    
                    if  sqrt(obj.X(i).^2+obj.Y(j).^2)<obj.r*1.001                    
                        psy(j,i)=1/0;                      
                    elseif obj.X(i)>obj.H                      
                        psy(j,i)=1/0;                      
                    else                      
                        psy(j,i) =  -obj.U*obj.X(i) + obj.calcul_courant(obj.X(i),obj.Y(j)) ;                      
                    end                 
                end             
            end              
            figure(1)
            contour(obj.X,obj.Y,psy,obj.nb_iso);
            pbaspect([1 1 1])
            hold on
            line([obj.H obj.H],[-10 10])            
            grid on     
            viscircles([0 0],obj.r,'color','k')                  
        end
%%%%%%%%%%%%%

%%%%%%%%%%%%% Vitesse de glissement
        function V = vel_gliss(obj)
            RC = obj.r*cos(obj.th);
            RS = obj.r*sin(obj.th);  
            K_glissement = zeros(obj.Ng,obj.N);
            for i=1:obj.Ng
                for j=1:obj.N
                    K_glissement(i,j) = (obj.Y(i)-RS(j))*((1/((obj.H-RC(j))^2+(obj.Y(i)-RS(j))^2))+...
                                        (1/((-obj.H-RC(j))^2+(obj.Y(i)-RS(j))^2)))...
                                       -(obj.Y(i)+RS(j))*((1/((obj.H-RC(j))^2+(obj.Y(i)+RS(j))^2))+...
                                        (1/((-obj.H-RC(j))^2+(obj.Y(i)+RS(j))^2)));
                end
            end   
            V = obj.U + trapz(obj.th,obj.sigma.*K_glissement,2);
            plot(obj.Y,1./V)           
        end
%%%%%%%%%%%%%

%%%%%%%%%%%%% Affichage de la vitesse de glissement EXACTE sans paroi
        function gliss_exact(obj)
            for i=1:obj.Ng
                v_gliss_exact(i) = obj.U*(1+obj.r*((obj.H*obj.H-obj.Y(i)*obj.Y(i))/...
                                   (obj.H*obj.H+obj.Y(i)*obj.Y(i))^2));
            end
            plot(obj.Y,v_gliss_exact)
        end
%%%%%%%%%%%%%

    end
end

