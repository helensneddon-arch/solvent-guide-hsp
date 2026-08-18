# solvent-guide-hsp
Python code for solvent finding, heatmap generation and circle diagram creation, with Parameter set up Excel files for both the pharmaceutical and paints and coatings solvent guides

The attached files are intended to support the publication "Expansion of a Paints and Coatings Solvent guide. Combining Solvent Guides and Hansen Solubility Parameters for Greener Solvent and Solvent Mixture Selection"
Erica Kemp, Harry Maslen, Helen F. Sneddon*

Published in RSC Sustainability, 2026 <insert doi>

SolventFinder.py can be used to generate csv files listing binary mixtures as potential mixtures to solvate a given substrate, or replace a given solvent, in order of increasing RED which can be filtered, for example by Red, Amber, Green score of the worst scoring component, or by biobased content etc.

RAG_Heatmap_Generator.py can be used to generate "clouds" of Red, Amber or Green binary mixtures of solvent in HSP space, similar to figures 10 and 11 in the paper.

Circular_diagram_generator.py can be used to generate diagrams showing binary solvent mixture alternatives for screening as potential mixtures to solvate a given substrate, or replace a given solvent, similar to figures 5-8 in the paper, using an input from the output of SolventFinder.py.

The Paper and supporting code are published under CC-BY

The authors acknowledge the use of Claude Opus 4.6 (Anthropic) for assistance with Python code development, and the code builds on and modifies the code from Albrecht and Strube shared in https://git.uibk.ac.at/c7511030/solvent-finder and F. Albrecht and O. I. Strube, Journal of Coatings Technology and Research, 2025, 22, 1263–1268.

For any queries which are not answered by the publication, or the supporting information to the publication, please contact helen.sneddon@york.ac.uk
