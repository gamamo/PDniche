### This is the repository of the article:

*Bird species contributions to plant phylogenetic diversity in mutualistic networks along a human impact gradient

### Contact:
- Gabriel M. Moulatlet (mandaprogabriel@gmail.com)
﻿
### Abstract: 
1. Human activities are increasingly altering biodiversity, reshaping the network of interactions that sustain ecological systems. Bird-plant seed-dispersal networks organization and function are fundamentally affected as impacts on either plant or bird species accumulate.  Yet it remains largely unexplored how human activities affect bird species’ contributions to plant phylogenetic diversity and whether these contributions driven by birds’ functional traits and by geographical and evolutionary properties
2. We measured bird species’ contributions to plant phylogenetic diversity as the bird species’ phylogenetic interaction niche (PDniche) and as its contribution to plant phylogenetic diversity through species interactions (PDcontrib). Moreover, we explored how human impacts, functional traits and geographical and evolutionary properties influence these metrics at different network levels, from individual local networks to the global meta-network representing interactions across all local networks.
3. Using 508 seed dispersal networks, across five continents and 11 biogeographical regions and a human impact gradient, we show that traits related to bird species dispersal, such as body mass, climatic niche position and range sizes, were significantly associated with birds’ PDniche and PDcontrib at both local network and global meta-network levels.
4. Human impacts showed contrary effect directions depending on the network level, being positively related to PDniche and PDcontrib at the local network and negatively related to birds’ PDniche at the meta-network level.
5. Taken together, our results showed that the traits determining bird species PDniche are relative to the levels of network organisation and to species trait specificity. Also, human impacts are an important driver of PDniche and PDcontrib, suggesting that predicting species contribution to the phylogenetic diversity of a community requires combined biological and human impact assessments.

### Responsible for writing code: 
- Gabriel M. Moulatlet (mandaprogabriel@gmail.com)
﻿
### Folders and files:
- There is one data folder and one code file.
﻿
	* Data folder:
		* 	`DataForAnalysis_v21jul26_encodingFix_utf_NoDupli.csv`, with the following columns:
    		* `Species`: Bird species
			* `plant_id`: Plant species
   			* `database`: individual network identification
    		* `ref`: Original reference from where the network was obtained
    		* `interaction`: 1 if the bird and plant species interact with each other and 0 if not.
    		* `lat` and `long`: latitude and longitude coordinates
    		* `Trophic.niche`: Bird species diet, information obtained from AVONET
    		* `Realm`: Biogegraphic realm affiliation of each species.
      	* `Avonet_realmFixed.csv`, with the following columns:
      		*`Species`: Bird species
      	    *`Trophic.niche`: Bird species diet, information obtained from AVONET
      	  	*`Realm`: Biogegraphic realm affiliation of each species.
    	
  		* 	Auxiliary data sources (which we are not allowed to make avilable) buy are used in our analysis AVONET, Biogeographic Realms and the Human Impact Index. They can be found in the orginal publications:
		
			- Tobias, J. A., C. Sheard, A. L. Pigot, A. J. M. Devenish, J. Yang, F. Sayol, M. H. C. Neate-Clegg, N. Alioravainen, T. L. Weeks, R. A. Barber, and Others. 2022. AVONET: morphological, ecological and geographical data for all birds. Ecology letters 25:581–597. Wiley Online Library.
    - 	 Holt, B. G., J.-P. Lessard, M. K. Borregaard, S. A. Fritz, M. B. Araújo, D. Dimitrov, P.-H. Fabre, C. H. Graham, G. R. Graves, K. A. Jønsson, and Others. 2013. An update of Wallace’s zoogeographic regions of the world. Science 339:74–78. American Association for the Advancement of Science.
			- Theobald, D.M., Oakleaf, J.R., Moncrieff, G., Voigt, M., Kiesecker, J. & Kennedy, C.M. (2025). Global extent and change in human modification of terrestrial ecosystems from
1990 to 2022. Scientific Data 12(489). [doi:10.1038/s41597-025-04892-2](https://doi.org/10.1038/s41597-025-04892-2)

	* Code file:
	    * `1-run_phylo_niche.R`: contains the codes to calculate the metrics PDniche and PDcontrib
      * `2-HumanImpacts_Overview.R`: contains the codes to calculate the overall human impacts on PDniche and PDcontrib
      * `3-AnalysisMCMC.R`: contains the codes to run MCMCglmm models
      * `4-Figures_MCMC.R`: contains the codes to obtains the figures and tables resulting of MCMCglmm models

* Software version:
   - R version 4.5.1 (2025-06-13 ucrt)
   - Platform: x86_64-w64-mingw32/x64
   - Running under: Windows 11 x64 (build 26200)
   - Rstudio 2026.08.2

### Acknowledgments:
- We are thankful to Tom Bradfer-Lawrence, Clementine Durand-Bessart, and Manuel Nogales for sharing their compiled datasets of seed dispersal networks. GM was supported by a SECIHTI (Mexico) postdoctoral grant (CVU 1135090). FV and WD thank the Instituto de Ecologia A.C. for continuous institutional support. SBM thanks the Dorothea Schlözer Postdoctoral Programme of the Georg-August-Universität Göttingen for their support. IM-C acknowledges funding from the Spanish Ministry for Science, Innovation and Universities (grant number PID2023-152329OB-I00).

