# Lyme-Disease-Model-in-Maryland

Mechanistic modeling of Lyme disease in Maryland under current and projected climate scenarios.

## 📁 Project Contents
- `codes/` contains R scripts for this study
- **`data/`** – Input datasets  
  (case reports, temperature records, county boundaries, population data) 

## Install required packages once
install.packages(c("deSolve","ggplot2","dplyr","patchwork","sf","tigris",
                   "readr","viridis","lhs","sensitivity"))

> ⏱ Note: Full simulation may take several minutes to complete, depending on system specs.

## 📚 References:

### Tick Biology & Ecology
- Radolf JD, Caimano MJ, Stevenson B, Hu LT. *Of ticks, mice and men: understanding the dual-host lifestyle of Lyme disease spirochaetes*. Nat Rev Microbiol. 2012;10(2):87–99.
- Eisen RJ, Eisen L, Beard CB. *County-scale distribution of Ixodes scapularis and Ixodes pacificus (Acari: Ixodidae) in the continental United States*. J Med Entomol. 2016;53(2):349–386.
- Steere, A. C., Strle, F., Wormser, G. P., et al. *Lyme borreliosis*. Nat Rev Dis Primers. 2016;2:16090.

### Climate-Vector Modeling & Lyme Disease
- Ogden NH et al. *Risk maps for range expansion of the Lyme disease vector, Ixodes scapularis, in Canada now and with climate change*. Int J Health Geogr. 2008;7:24.
- Worton AJ et al. *GIS-ODE: linking dynamic population models with GIS to predict pathogen vector abundance across a country under climate change scenarios*. J R Soc Interface. 2024;21(217):20240004.
- Couper LI, MacDonald AJ, Mordecai EA. *Impact of prior and projected climate change on US Lyme disease incidence*. Glob Change Biol. 2021;27(4):738–754.
- Bayoh MN, Thomas CJ, Lindsay SW. *Mapping distributions of chromosomal forms of Anopheles gambiae in West Africa using climate data*. Med Vet Entomol. 2001;15(3):267–74.

### Mathematical Modeling of Vector-Borne Diseases
- Eikenberry SE, Gumel AB. *Mathematics of malaria and climate change*. In: Mathematics of Planet Earth: Protecting Our Planet, Learning from the Past, Safeguarding for the Future. Cham: Springer International Publishing; 2019. pp. 77-108.
- Okuneye K, Eikenberry SE, Gumel AB. *Weather-driven malaria transmission model with gonotrophic and sporogonic cycles*. J Biol Dyn. 2019;13(sup1):288-324.
- Okuneye K, Gumel AB. *Analysis of a temperature-and rainfall-dependent model for malaria transmission dynamics*. Math Biosci. 2017 May;287:72-92.
- Husar K, Pittman DC, Rajala J, Mostafa F, Allen LJ. *Lyme disease models of tick-mouse dynamics with seasonal variation in births, deaths, and tick feeding*. Bull Math Biol. 2024 Mar;86(3):25.

