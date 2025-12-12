# Mechanistic-Model-of-Lyme-Disease-in-Maryland

Mechanistic modeling of Lyme disease in Maryland under current and projected climate scenarios.

## Content
- `codes/` — R scripts relevant to this study  
- `data/`  /` — Maryland temperature, county cases, population, land area  

## Install required packages once
install.packages(c("deSolve","ggplot2","dplyr","patchwork","sf","tigris",
                   "readr","viridis","lhs","sensitivity"))

> ⏱ Note: Full simulation may take several minutes to complete, depending on system specs.

## 📚 References:
- Ogden NH et al. *Risk maps for range expansion of the Lyme disease vector, Ixodes scapularis, in Canada now and with climate change*. Int J Health Geogr. 2008;7:24.  
- Worton AJ et al. *GIS-ODE: linking dynamic population models with GIS to predict pathogen vector abundance across a country under climate change scenarios*. J R Soc Interface. 2024;21(217):20240004.  
- Couper LI, MacDonald AJ, Mordecai EA. *Impact of prior and projected climate change on US Lyme disease incidence*. Glob Change Biol. 2021;27(4):738–754.
- Radolf JD, Caimano MJ, Stevenson B, Hu LT. *Of ticks, mice and men: understanding the dual-host lifestyle of Lyme disease spirochaetes*. Nat Rev Microbiol. 2012;10(2):87–99.  
- Eisen RJ, Eisen L, Beard CB. *County-scale distribution of Ixodes scapularis and Ixodes pacificus (Acari: Ixodidae) in the continental United States*. J Med Entomol. 2016;53(2):349–386.
- Mordecai EA, Caldwell JM, Grossman MK, et al. Thermal biology of mosquito‐borne disease. Ecol Lett. 2019;22(10):1690–708.
- Bayoh MN, Thomas CJ, Lindsay SW. Mapping distributions of chromosomal forms of Anopheles gambiae in West Africa using climate data. Med Vet Entomol. 2001;15(3):267–74.
