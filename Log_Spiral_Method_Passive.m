clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Introduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Enok Cheon, Ph.D. Student in UIUC (Supervisor: Prof Timothy D. Stark)
Date: Completed on 23rd Feb, 2018
Description:
    The script was written to primarily solve CEE484 HW3 problem. The script allows users to 
    find a passive earth pressure of wall with nonlinear failure plain based log spiral method
    as described in the textbook by Terzaghi, Peck, Mesri. 
    The script allows to calculate passive pressure from the pressure from (1) earth, (2) water,
    (3) surcharge and (4) cohesion intercept for both (1) drained and (2) undrained condition.
    The script iterates the position of ptO and computes passive pressure force (Pp). Next, it 
    finds the critical, i.e. minimum, value of Pp.
    The script generates plots and export the critical computed values
Assumption and Limitations:
    (1) inclination of groundwater (GW) level and soil behind the wall is horizontal
    (2) wall inclination in vertical (perpendicular to the wall height) and the wall is not
    	embedded into the soil
    (3) the surcharge load is uniform surcharge applied; no point load or line load
    (4) Pp from earth pressure is assumed to apply on height above the wall base at value of 
        (1/3) of wall height
    (5) in drained condition, interface friction angle = (2/3)*(effective friction angle of soil)
        and wall adhesion is (2/3)*(effective cohesion intercept)
    (6) in undrained condition, wall adhesion is (2/3)*(undrained shear strength)   
    (7) in undrained condition, the soil in backfill are assumed to be all 
        saturated above and below the watertable. The water table level does not change. 
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wall dimension
H = 15;     % wall height - unit: ft

% soil properties
% drained
gamma_s = 130;	% unit weight - unit: pcf
phi_eff = 35;	% drained friction angle
c_eff = 0;		% drained cohesion - unit: psf

% undrained
gamma_s_sat = 130; % unit weight - unit: pcf
su = 2000;		% undrained shear strength - unit: psf

% ground water (GW)
hw = 10;	% wall height from base Pt A with GW level - unit: ft
gamma_w = 62.4;  % unit weight - unit: pcf

% unform surcharge
q = 1000;	% unit: psf

% anaylsis method
analysisMethod = 1;
% 0 = drained
% 1 = undrained

% steps 
startXforO = -15;
endXforO = 0;
stepsXforO = 0.5;
stepsXforLog = 1;

% format of numbering
format long
significantFigure = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pt A and B
ptA = [0,0];
ptB = [0,H];

dx = [startXforO:stepsXforO:endXforO];

%% drained analysis
if analysisMethod == 0

    Kp = (1+sind(phi_eff))/(1-sind(phi_eff)); % Rankine Passive Earth Pressure
    delta_w = (2/3)*phi_eff;    % drained - interface friction angle between soil and wall
    ad_w = (2/3)*c_eff;         % undrained - adhesion of soil and wall

    summaryData = []; % ptO, ptC, thetaMax, radius0, radiusMax, Pp_s, Pp_q, Pp_c, Pp
	for dxN = 1:length(dx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% Geometry
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pt 0
		ptO = [dx(dxN), H-dx(dxN)*tand(45-0.5*phi_eff)];
		
		% radius initial
		radius0 = sqrt(ptO(1)^2 + ptO(2)^2);

		% Find thetaMax
		vectorOA = ptO;
        vectorAX = [200, 0];
        lengthOA = sqrt(vectorOA(1)^2 + vectorOA(2)^2);
        lengthAX = sqrt(vectorAX(1)^2 + vectorAX(2)^2);
        thetaOA = acosd(dot(vectorOA,vectorAX)/(lengthOA*lengthAX));
		thetaMax = 180 - (thetaOA+45-(phi_eff/2));

		% radius final
		radiusMax = radius0*exp(tand(phi_eff)*thetaMax*pi/180);

		% ptC
		syms ptCx ptCy
        eqns = [radiusMax == sqrt((ptO(1)-ptCx).^2 + (ptO(2)-ptCy).^2),...
                ptCy == H-ptCx*tand(45-0.5*phi_eff)];
        S = solve(eqns, [ptCx ptCy]); 
        ptCAns = double(vpa([S.ptCx, S.ptCy]));
        
        if ptCAns(1,1) >= ptO(1)
            ptC = ptCAns(1,:);
        elseif ptCAns(2,1) >= ptO(1)
            ptC = ptCAns(2,:);
        end

        hr = H - ptC(2); % height of Rankine section

        % log spiral location between ptB and ptC
        dTheta = [0:stepsXforLog:floor(thetaMax), thetaMax];
        logSpiralPt = zeros(length(dTheta),5);
        for dThetaN = 1:length(dTheta)
            radiusTheta = radius0*exp(tand(phi_eff)*dTheta(dThetaN)*pi/180);
            deltax = radiusTheta*cosd(-(thetaMax-dTheta(dThetaN)+45-(phi_eff/2)));
            deltay = radiusTheta*sind(-(thetaMax-dTheta(dThetaN)+45-(phi_eff/2)));
            xcoord = ptO(1)+deltax;
            ycoord = ptO(2)+deltay;
        	logSpiralPt(dThetaN,:) = [dThetaN, dTheta(dThetaN), radiusTheta, xcoord, ycoord];
        end
        
        dw = H - hw; % GW level below the top of backfill

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pp soil
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rankine
        Pr1_s = 0.5*gamma_s*(dw^2)*Kp;
        Pr2_s = gamma_s*dw*(hr-dw)*Kp;
        Pr3_s = 0.5*(gamma_s-gamma_w)*((hr-dw)^2)*Kp;
        Pr_s = Pr1_s + Pr2_s + Pr3_s;
        
        Mr_s = Pr1_s*(hr-(2*dw/3)) + Pr2_s*0.5*(hr-dw) + Pr3_s*(hr-dw)/3;
        Lr_s = ptO(2)-((Mr_s/Pr_s)+ptC(2));

        %% Log Spiral
        % weight of soil above GW
        Wa_l = gamma_s*ptC(1)*dw;         % weight
        Ll_wa = 0.5*ptC(1) - ptO(1);   % distance from ptO

        % weight of soil below GW
        sumA = 0;
        sumAx = 0;
        for dThetaN = 2:length(dTheta)
            area_rec = (logSpiralPt(dThetaN,4) - logSpiralPt(dThetaN-1,4))*(H-dw-logSpiralPt(dThetaN-1,5));
            area_tri = 0.5*(logSpiralPt(dThetaN,4) - logSpiralPt(dThetaN-1,4))*abs(logSpiralPt(dThetaN,5) - logSpiralPt(dThetaN-1,5));
            sumA = sumA + area_tri + area_rec;

            len_rec = 0.5*(logSpiralPt(dThetaN,4) + logSpiralPt(dThetaN-1,4));
            len_tri = logSpiralPt(dThetaN-1,4) + (2/3)*(logSpiralPt(dThetaN,4) - logSpiralPt(dThetaN-1,4));
            sumAx = sumAx + area_rec*len_rec + area_tri*len_tri;
        end
        Wb_l = (gamma_s-gamma_w)*sumA;
        Ll_wb = (sumAx/sumA) - ptO(1);

        %% Pp on wall for soil
        % assume that passive force on wall applied on H/3 above the base of the wall
        L_sy = ptO(2) - (H/3);  
        Pp_s = (Wa_l*Ll_wa + Wb_l*Ll_wb + Pr_s*Lr_s)/(L_sy*cosd(delta_w)+ptO(1)*sind(delta_w));
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pp surcharge
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rankine
        Qr = q*Kp*hr;
        Lr_q = ptO(2) - (hr*0.5 + ptC(2));

        %% Log Spiral
        Ql = q*(ptC(1));
        Ll_q = 0.5*ptC(1) - ptO(1);

        %% Pp on wall for surcharge
        % assume that passive force on wall applied on H/2 above the base of the wall
        L_q = ptO(2) - (H/2);  
        Pp_q = (Qr*Lr_q + Ql*Ll_q)/(L_q*cosd(delta_w)+ptO(1)*sind(delta_w));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pp cohesion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rankine
        Pr_c = 2*c_eff*sqrt(Kp)*hr;
        Lr_c = ptO(2) - (hr*0.5 + ptC(2));

        %% Log Spiral
        Ml_c = (c_eff/(2*tand(phi_eff)))*(radiusMax^2 - radius0^2);

        %% adhesion on wall
        Ca = ad_w*H;

        %% Pp on wall for cohesion
        L_c = ptO(2) - (H/2);  
        Pp_c = (Pr_c*Lr_c + Ml_c - Ca*ptO(1))/(L_c*cosd(delta_w)+ptO(1)*sind(delta_w));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Summary of data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pp = Pp_s + Pp_q + Pp_c;
        summaryData(dxN,:) = [ptO, ptC, thetaMax, radius0, radiusMax, Pp_s, Pp_q, Pp_c, Pp];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
        hold on
        plot([ptA(1),ptB(1)],[ptA(2),ptB(2)],'k-','linewidth',5)
        plot(ptO(1),ptO(2),'kx')
        plot(ptC(1),ptC(2),'ko')
        plot(logSpiralPt(:,4),logSpiralPt(:,5),'g--')
        plot([-15,30], [H-(-15)*tand(45-0.5*phi_eff), H-(30)*tand(45-0.5*phi_eff)], 'k--');
        xlim([-15,30])
        xlabel('x coordiante [ft]')
        ylabel('y coordiante [ft]')
        pbaspect([1 1 1])
        
	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% critical Pp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the minimum value of computed Pp
    sorted_summaryData = sortrows(summaryData, 11);
    crit_summaryData = sorted_summaryData(1,:);
    
    % critical log spiral section
    dTheta = [0:stepsXforLog:floor(crit_summaryData(5)), crit_summaryData(5)];
    crit_logSpiralPt = zeros(length(dTheta),5);
    for dThetaN = 1:length(dTheta)
        radiusTheta = crit_summaryData(6)*exp(tand(phi_eff)*dTheta(dThetaN)*pi/180);
        deltax = radiusTheta*cosd(-(crit_summaryData(5)-dTheta(dThetaN)+45-(phi_eff/2)));
        deltay = radiusTheta*sind(-(crit_summaryData(5)-dTheta(dThetaN)+45-(phi_eff/2)));
        xcoord = crit_summaryData(1)+deltax;
        ycoord = crit_summaryData(2)+deltay;
        crit_logSpiralPt(dThetaN,:) = [dThetaN, dTheta(dThetaN), radiusTheta, xcoord, ycoord];
    end

    % account for water pressure on wall
    Pw = 0.5*gamma_w*hw^2;
    Lm_w = hw/3;

    % resultant moment with pivoting at ptA
    Pp_s = crit_summaryData(8);
    Pp_q = crit_summaryData(9);
    Pp_c = crit_summaryData(10);

    Lm_s = H/3;
    Lm_q = H/2;
    Lm_c = H/2;

    sumM_ptA = Pw*Lm_w + Pp_s*Lm_s + Pp_q*Lm_q + Pp_c*Lm_c;
    sumF = Pw + Pp_s + Pp_q + Pp_c;
    Lm = sumM_ptA/sumF;

    final_summaryData = [crit_summaryData(1:7), Pp_s, Lm_s, Pp_q, Lm_q, Pp_c, Lm_c, Pw, Lm_w, crit_summaryData(11), sumF, Lm, sumM_ptA];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot and display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(final_summaryData(1),final_summaryData(2),'rx')
    plot(final_summaryData(3),final_summaryData(4),'ro')
    plot(crit_logSpiralPt(:,4),crit_logSpiralPt(:,5),'r-')
    plot([-15,0],[0,0],'k-')
    plot([0,60],[H,H],'k-')
    legend('wall','trial ptO','trial ptC','trial log spiral','45-(phi''/2) line')

    figure(2)
    hold on
    plot([ptA(1),ptB(1)],[ptA(2),ptB(2)],'k-','linewidth',5)
    plot(final_summaryData(1),final_summaryData(2),'rx')
    plot(final_summaryData(3),final_summaryData(4),'ro')
    plot(crit_logSpiralPt(:,4),crit_logSpiralPt(:,5),'r-')
    plot([-15,30], [H-(-15)*tand(45-0.5*phi_eff), H-(30)*tand(45-0.5*phi_eff)], 'k--');
    plot([-15,0],[0,0],'k-')
    plot([0,60],[H,H],'k-')
    legend('wall','crit. ptO','crit. ptC','crit. log spiral','45-(phi''/2) line')
    xlim([-15,30])
    xlabel('x coordiante [ft]')
    ylabel('y coordiante [ft]')
    pbaspect([1 1 1])

    figure(3)
    hold on
    plot(summaryData(:,1),summaryData(:,11),'k-')
    plot(final_summaryData(1),final_summaryData(16),'rx')
    legend('Pp','critical Pp')
    xlabel('x coordiante [ft]')
    ylabel('Pp [lbs]')
%     xlim([-15,-10])

    % disp(final_summaryData)
    disp('critical horizontal distance of ptO from the wall in drained condition: (positive horizontal direction is toward right)')
    disp(final_summaryData(1))
    disp('resultant passive force (earth and water) on the wall in drained condition:')
    disp(final_summaryData(16))
    disp('resultant moment about base of the wall due to resultant passive force (earth and water) on the wall in drained condition:')
    disp(final_summaryData(19))

%% undrained analysis
elseif analysisMethod == 1

    Kp = 1;     % Rankine Passive Earth Pressure
    ad_w = (2/3)*su;    % undrained - adhesion of soil and wall

    summaryData = []; % ptO, ptC, thetaMax, radius0, Pp_s, Pp_q, Pp_c, Pp
    for dxN = 1:length(dx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Geometry
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pt 0
        ptO = [dx(dxN), H-dx(dxN)];
        
        % radius initial
        radius0 = sqrt(ptO(1)^2 + ptO(2)^2);

        % Find thetaMax
        vectorOA = ptO;
        vectorAX = [200, 0];
        lengthOA = sqrt(vectorOA(1)^2 + vectorOA(2)^2);
        lengthAX = sqrt(vectorAX(1)^2 + vectorAX(2)^2);
        thetaOA = acosd(dot(vectorOA,vectorAX)/(lengthOA*lengthAX));
        thetaMax = 180 - (thetaOA+45);

        % ptC
        syms ptCx ptCy
        eqns = [radius0 == sqrt((ptO(1)-ptCx).^2 + (ptO(2)-ptCy).^2),...
                ptCy == H-ptCx];
        S = solve(eqns, [ptCx ptCy]); 
        ptCAns = double(vpa([S.ptCx, S.ptCy]));
        
        if ptCAns(1,1) >= ptO(1)
            ptC = ptCAns(1,:);
        elseif ptCAns(2,1) >= ptO(1)
            ptC = ptCAns(2,:);
        end

        hr = H - ptC(2); % height of Rankine section

        % circular location between ptB and ptC
        dTheta = [0:stepsXforLog:floor(thetaMax), thetaMax];
        cirPt = zeros(length(dTheta),4);
        for dThetaN = 1:length(dTheta)
            deltax = radius0*cosd(-(thetaMax-dTheta(dThetaN)+45));
            deltay = radius0*sind(-(thetaMax-dTheta(dThetaN)+45));
            xcoord = ptO(1)+deltax;
            ycoord = ptO(2)+deltay;
            cirPt(dThetaN,:) = [dThetaN, dTheta(dThetaN), xcoord, ycoord];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pp soil
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rankine
        Pr_s = 0.5*(gamma_s_sat-gamma_w)*(hr^2);
        Lr_s = ptO(2)-((hr/3)+ptC(2));

        %% circular section
        % weight of soil
        sumA = 0;
        sumAx = 0;
        for dThetaN = 2:length(dTheta)
            area_rec = (cirPt(dThetaN,3) - cirPt(dThetaN-1,3))*(H-cirPt(dThetaN-1,4));
            area_tri = 0.5*(cirPt(dThetaN,3) - cirPt(dThetaN-1,3))*abs(cirPt(dThetaN,4) - cirPt(dThetaN-1,4));
            sumA = sumA + area_tri + area_rec;

            len_rec = 0.5*(cirPt(dThetaN,3) + cirPt(dThetaN-1,3));
            len_tri = cirPt(dThetaN-1,3) + (2/3)*(cirPt(dThetaN,3) - cirPt(dThetaN-1,3));
            sumAx = sumAx + area_rec*len_rec + area_tri*len_tri;
        end
        W_s = (gamma_s_sat-gamma_w)*sumA;
        Lc_s = (sumAx/sumA) - ptO(1);

        %% Pp on wall for soil
        % assume that passive force on wall applied on H/3 above the base of the wall
        L_sy = ptO(2) - (H/3);  
        Pp_s = (W_s*Lc_s + Pr_s*Lr_s)/L_sy;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pp surcharge
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rankine
        Qr = q*hr;
        Lr_q = ptO(2) - (hr*0.5 + ptC(2));

        %% circular section
        Qc = q*(ptC(1));
        Lc_q = 0.5*ptC(1) - ptO(1);

        %% Pp on wall for surcharge
        % assume that passive force on wall applied on H/2 above the base of the wall
        L_q = ptO(2) - (H/2);  
        Pp_q = (Qc*Lc_q + Qr*Lr_q)/L_q;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pp cohesion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rankine
        Pr_c = 2*su*hr;
        Lr_c = ptO(2) - (hr*0.5 + ptC(2));

        %% Log Spiral
        Mc_c = su*(radius0^2)*thetaMax*(pi/180);

        %% adhesion on wall
        Ca = ad_w*H;

        %% Pp on wall for cohesion
        L_c = ptO(2) - (H/2);  
        Pp_c = (Pr_c*Lr_c + Mc_c - Ca*ptO(1))/L_c;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Summary of data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pp = Pp_s + Pp_q + Pp_c;
        summaryData(dxN,:) = [ptO, ptC, thetaMax, radius0, Pp_s, Pp_q, Pp_c, Pp];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
        hold on
        plot([ptA(1),ptB(1)],[ptA(2),ptB(2)],'k-','linewidth',5)
        plot(ptO(1),ptO(2),'kx')
        plot(ptC(1),ptC(2),'ko')
        plot(cirPt(:,3),cirPt(:,4),'g--')
        plot([-15,30], [H-(-15), H-(30)], 'k--');
        xlim([-15,15])
        ylim([-5,30])
        xlabel('x coordiante [ft]')
        ylabel('y coordiante [ft]')
        pbaspect([1 1 1])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% critical Pp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the minimum value of computed Pp
    sorted_summaryData = sortrows(summaryData, 10);
    crit_summaryData = sorted_summaryData(1,:);
    
    % critical circular section
    dTheta = [0:stepsXforLog:floor(crit_summaryData(5)), crit_summaryData(5)];
    crit_cirPt = zeros(length(dTheta),4);
    for dThetaN = 1:length(dTheta)
        deltax = crit_summaryData(6)*cosd(-(crit_summaryData(5)-dTheta(dThetaN)+45));
        deltay = crit_summaryData(6)*sind(-(crit_summaryData(5)-dTheta(dThetaN)+45));
        xcoord = crit_summaryData(1)+deltax;
        ycoord = crit_summaryData(2)+deltay;
        crit_cirPt(dThetaN,:) = [dThetaN, dTheta(dThetaN), xcoord, ycoord];
    end

    % account for water pressure on wall
    Pw = 0.5*gamma_w*hw^2;
    Lm_w = hw/3;

    % resultant moment with pivoting at ptA
    Pp_s = crit_summaryData(7);
    Pp_q = crit_summaryData(8);
    Pp_c = crit_summaryData(9);

    Lm_s = H/3;
    Lm_q = H/2;
    Lm_c = H/2;

    sumM_ptA = Pw*Lm_w + Pp_s*Lm_s + Pp_q*Lm_q + Pp_c*Lm_c;
    sumF = Pw + Pp_s + Pp_q + Pp_c;
    Lm = sumM_ptA/sumF;

    final_summaryData = [crit_summaryData(1:6), Pp_s, Lm_s, Pp_q, Lm_q, Pp_c, Lm_c, Pw, Lm_w, crit_summaryData(10), sumF, Lm, sumM_ptA];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot and display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(final_summaryData(1),final_summaryData(2),'rx')
    plot(final_summaryData(3),final_summaryData(4),'ro')
    plot(crit_cirPt(:,3),crit_cirPt(:,4),'r-')
    plot([-15,0],[0,0],'k-')
    plot([0,60],[H,H],'k-')
    legend('wall','trial ptO','trial ptC','trial log spiral','45-(phi''/2) line')

    figure(2)
    hold on
    plot([ptA(1),ptB(1)],[ptA(2),ptB(2)],'k-','linewidth',5)
    plot(final_summaryData(1),final_summaryData(2),'rx')
    plot(final_summaryData(3),final_summaryData(4),'ro')
    plot(crit_cirPt(:,3),crit_cirPt(:,4),'r-')
    plot([-15,30], [H-(-15), H-(30)], 'k--')
    plot([-15,0],[0,0],'k-')
    plot([0,60],[H,H],'k-')
    legend('wall','crit. ptO','crit. ptC','crit. log spiral','45-(phi''/2) line')
    xlim([-15,15])
    ylim([-5,30])
    xlabel('x coordiante [ft]')
    ylabel('y coordiante [ft]')
    pbaspect([1 1 1])

    figure(3)
    hold on
    plot(summaryData(:,1),summaryData(:,10),'k-')
    plot(final_summaryData(1),final_summaryData(15),'rx')
    legend('Pp','critical Pp')
    xlabel('x coordiante [ft]')
    ylabel('Pp [lbs]')
%     xlim([-7,-3])

%     disp(final_summaryData)
    
    disp('critical horizontal distance of ptO from the wall in undrained condition: (positive horizontal direction is toward right)')
    disp(final_summaryData(1))
    disp('resultant passive force (earth and water) on the wall in undrained condition:')
    disp(final_summaryData(16))
    disp('resultant moment about base of the wall due to resultant passive force (earth and water) on the wall in undrained condition:')
    disp(final_summaryData(18))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analysisMethod == 0
    dlmwrite(strcat('CEE484_HW3_summary Data of Log Spiral Passive_drained.csv'),summaryData, 'delimiter', ',', 'precision', significantFigure);
    dlmwrite(strcat('CEE484_HW3_summary Data of critical Log Spiral Passive_drained.csv'), final_summaryData, 'delimiter', ',', 'precision', significantFigure);
    dlmwrite(strcat('CEE484_HW3_cooridante of Log Spiral points_drained.csv'), crit_logSpiralPt, 'delimiter', ',', 'precision', significantFigure);
elseif analysisMethod == 1
    dlmwrite(strcat('CEE484_HW3_summary Data of Log Spiral Passive_undrained.csv'),summaryData, 'delimiter', ',', 'precision', significantFigure);
    dlmwrite(strcat('CEE484_HW3_summary Data of critical Log Spiral Passive_undrained.csv'), final_summaryData, 'delimiter', ',', 'precision', significantFigure);
    dlmwrite(strcat('CEE484_HW3_cooridante of Circular points_undrained.csv'), crit_cirPt, 'delimiter', ',', 'precision', significantFigure);
end

disp('done')
