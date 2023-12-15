# Lunar-Gravity-Assist-Patched-Conic-Simulation
This repo contains the MATLAB code used as part of my masters using a novel patched conics approximation method the V infinity Sphere to look at earth-earth lunar gravity assists.

The aim of this project was to analyse the use for Lunar gravity assists to transfer satellites from initial earth orbits available when launching from proposed sites in the UK to a Geostationary Orbit. This was done using the patched conics approximation and a set of equations called the "V infinity sphere" to plan and assess the specific parameters needed for a lunar gravity assists to move from a intercept orbit to a target orbit with a perigee at geostationary altitude and zero degrees inclination. Then using the vector difference between the orbits used the delta V needed to execute a mission plan is calculated to be compared to transfer methods not using lunar gravity assists and asses the feasibility of such a lunar assist mission.

This program uses MATLAB to take the initial conditions from three launch sites in the UK which based on information available from sources about each site gives an initial orbital parameter a potential payload with initial orbital parameters as follows,

1. Sutherland Spaceport, 300x300Km at 83 degrees
2. Spaceport Cornwall, 230x230KM at 70 degrees
3. SaxaVord Spaceport, 300x300 at 62 degrees

From these initial orbits a lunar intercept orbit is planned, then using the lunar intercept parameters the post gravity assist orbit is calculated using the V infinity method the post gravity assist orbit is calculated, finally a circularised geostationary orbit is drawn. The vector velocity difference for the shared point between each orbit is found as the delta V needed to complete the burn. This script will generate the V-infinity sphere graph used to show the two different gravity assist approaches and then a 3D graph showing the orbital paths of each orbit used.

To operate the simulation MATLAB version R2022b will be needed.
1. Download MSc_Project_Propogator_MK2 and supporting functions (COE2RV, RV2COE and TRUPRO)
2. Ensure file pathways allow the main propagator to be able to access the support functions
3. Run main propagator file comments within the file outline operation to generate data and graphs for each launch site

Important data generated is the Transfer time for the Anti-planet assist and planet assist shown as Transfer_TimeS and Transfer_TimeL respectively the delta V needed to complete the mission for both assist types shown as SDV_1 and LDV_1 for the TLI burn and SDV_2 and LDV_2 for the geostationary circularisation burn with the mission totals as L and S DV_T.
