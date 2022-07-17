# Jet V2 Measureument
Tristan Protzman
Lehigh University
December 28, 2021

# Code Design
* As usual, the code is divided into two parts
* `qa` processes PicoDST files to populate the histograms for QA and the trees for the analysis
* `analysis` processes the tree generated by `qa` to extract the $v_{2}$ signal

## QA Details
* Histograms:
    * $v_z$ vs $v_{z,\;zdc}$
    * $v_r$
    * Track $p_T$
    * Track $\eta$ and $\phi$
    * Centrality16
    * RefMult3 vs TOFMult
    * BBC Rate, East and West
    * nMips vs RefMult3
* Event Plane Histograms - All second order
    * East/West uncorrected
    * East/West psi corrected
    * East/West psi and shifted
    * Difference in EP angle vs centrality
* Event Tree:
    * $v_x$, $v_y$, $v_z$, $v_{z,\;vpd}$
    * RefMult3
    * TOFMatch
    * TOFMult
    * BBC rate, East/West
    * Centrality16
    * Fully corrected event plane, East and West
* Jet Tree:
    * Number of jets
    * Number of constituents
    * Jet $p_T$
    * Median subtracted $p_T$
    * $\eta$ and $\phi$
    * Jet energy
    * Charged z
    * Jet area



# To Do
* `qa.cxx:144`: implement pileup cuts correctly
* `qa.cxx:181`: difference in angle not calculated correctly
* `analysis`: Seg faulting, no clue why :(
