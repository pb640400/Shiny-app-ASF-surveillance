This repository contains code for an R Shiny app designed to evaluate barn-level pre-movement active surveillance protocols for movement to a harvest facility of finishing pigs located within a Control Area established in response to an outbreak of African Swine Fever (ASF).

Active surveillance protocols are evaluated using a within-barn disease transmission and active surveillance simulation model. Users may choose from a variety of ASF virus strain and within-barn virus transmission scenarios. Additionally, users have multiple active surveillance protocol options including different numbers of samples and testing days prior to pig movement. Simulation output is presented to the app user in a variety of graphs and tables.

The file App.R is the primary file. The file transmission-asf-debug-gammaall.c must be compiled prior to use. The code depends on multiple R packages including shiny, shinydashboard, shinyWidgets, shinydashboardPlus, shinyBS, shinyjs,shinythemes, ggplot2, dplyt, DT, reshape2, and lubridate.

The app may be accessed by the link below:
https://sumn.shinyapps.io/ASF_hetl_g1/

\n
For details on the disease transmission simulation model, see: 
Ssematimba A, Malladi S, Bonney PJ, Charles KM, Boyer TC, Goldsmith T, Cardona CJ, Corzo CA, Culhane MR. African swine fever detection and transmission estimates using homogeneous versus heterogeneous model formulation in stochastic simulations within pig premises. Open Veterinary Journal. 2022;12(6):787-96.

For details on the active surveillance simulation model, see: 
Bonney PJ, Malladi S, Ssematimba A, O’Hara KC, Remmenga MD, Farr M, Leonard M, Alexander CY, Blair B, Martin SW, Culhane MR. Simulation of Premovement Active Surveillance Protocols for Moving Finishing Pigs to a Harvest Facility from a Control Area during an Outbreak of African Swine Fever in the United States. Transboundary and Emerging Diseases. 2024;2024(1):6657600.
