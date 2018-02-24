# log-spiral-method---passive-force-for-retaining-structure

Summary: MATLAB script to run log spiral method to calculate the critical passive pressure for retaining structure

Author: Enok Cheon, Ph.D. Student in UIUC (Supervisor: Prof Timothy D. Stark)

Date: Completed on 23rd Feb, 2018

Description:

The script was written to primarily solve CEE484 HW3 problem. The script allows users to find a passive earth pressure of wall with nonlinear failure plain based log spiral method as described in the textbook by Terzaghi, Peck, Mesri. The script allows to calculate passive pressure from the pressure from (1) earth, (2) water, (3) surcharge and (4) cohesion intercept for both (1) drained and (2) undrained condition. The script iterates the position of ptO and computes passive pressure force (Pp). Next, it finds the critical, i.e. minimum, value of Pp. The script generates plots and export the critical computed values
    
Assumption and Limitations:

(1) inclination of groundwater (GW) level and soil behind the wall is horizontal

(2) wall inclination in vertical (perpendicular to the wall height) and wall is not embedded into the ground

(3) the surcharge load is uniform surcharge applied; no point load or line load

(4) Pp from earth pressure is assumed to apply on height above the wall base at value of (1/3) of wall height

(5) in drained condition, interface friction angle = (2/3)*(effective friction angle of soil) and wall adhesion is (2/3)*(effective cohesion intercept)

(6) in undrained condition, wall adhesion is (2/3)*(undrained shear strength)

(7) in undrained condition, the soil in backfill are assumed to be all saturated above and below the watertable. The water table level does not change.
