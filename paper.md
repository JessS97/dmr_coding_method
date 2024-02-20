title: 'DNA Data Storage: Error Correction with Dynamic Mapping Rule'
tags:
  - DNA data storage
  - error correction
  - base mapping
  - reed solomon
  - random mask
  - 
authors:
  - name: Adrian M. Price-Whelan
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Author with no affiliation
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - given-names: Ludwig
    dropping-particle: van
    surname: Beethoven
    affiliation: 3
affiliations:
 - name: Fraunhofer Institute for Toxicology and Experimental Medicine (ITEM), Hannover, Germany
   index: 1
date: 23 January 2023
bibliography: paper.bib
---

# Summary

- DNA storage
- Fehlerkorrektur
- DMR Verfahren mehr eingehen

# Statement of need

- DNA speicher 
- verschiedene Fehlerkorrekturverfahren vorhanden
- möglichst gute Korrekturmethode finden 
- andere nennen --> meist nur subs fehler, oft mit RS gekoppelt
- 
DMR is a method for encoding bits in DNA and then decoding them. The method is based on a mapping table, which is used for the translation and also gives 
the dynamic mapping rule its name. With the help of this mapping table, a self-correction of the DNA sequence is to be achieved. In addition, the coding 
was combined with RS so that the correction with the DMR scheme can be validated with the help of RS and thus leads to a better correction than the stand-alone 
use of DMR or RS. 
In addition, the use of the DMR scheme offers the possibility to localize incorrect areas of the DNA by comparison with the mapping possibilities and thus 
to carry out targeted corrective procedures at the affected sites. This should make it possible to correct substitution errors as well as insertion and 
deletion errors. For the correction, the levels must be adapted and expanded accordingly so that an optimal correction can be achieved. The current correction 
is only suitable for the correction of substitution errors, but it is possible to extend it for the correction of insertion and deletion errors. 


- Bild: Edit RS und DMR Vergleich \autoref{fig:error_correction}


# Mathematics



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Comparison of DMR and RS correction with different settings and error rates of the bitarray segmentation coding method.\label{fig:error_correction}]
(DMR(n=20) and RS(n=50) correction with bitarray segmentation.png){ width=100% }


# Acknowledgements

This work was funded by the Fraunhofer PREPARE Programm and developed within the BIOSYNTH project under the project number 40-03168-2420 as part of the data 
coding research.

# References