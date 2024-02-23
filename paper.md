---
title: 'DNA Data Storage: Error Correction with Dynamic Mapping Rule'
tags:
  - DNA data storage
  - error correction
  - base mapping
  - reed solomon
  - random mask

authors:
  - name: Jessica Schindler
    equal-contrib: true
    affiliation: 1 
  - name: Babak Saremi
    equal-contrib: true
    affiliation: 1
  - name: Omar El Menuawy
    affiliation: 1
  - name: Tamara Hadzic
    affiliation: 1
  - name: Lena Wiese
    corresponding: true 
    affiliation: 1

affiliations:
 - name: Fraunhofer Institute for Toxicology and Experimental Medicine (ITEM), Hannover, Germany
   index: 1

date: 23 January 2023
bibliography: paper.bib
---

# Summary

The amount of generated data has multiplied rapidly in recent years. Conventional storage methods will eventually reach their limits and may no longer 
be sufficient to store the amount of data. Therefore, a new long-term storage method that can primarily store archival data is needed in order to be able 
to store the increasing amount of data in the coming years. DNA offers the advantages of long-term storage, high storage density and 
low energetic cost, causing the scientific community to attempt to make DNA usable as a storage method. 
This program contains the coding and error correction procedure of our self-developed "dynamic mapping rule" DMR method for storing data in DNA. 
The error correction of the method is based on the comparison of the faulty sequence with the mapping scheme and a level-based correction, which specifically 
addresses the known locations where the mapping does not match. By combining this with the Reed Solomon (RS) code, the valid segments of the DMR correction can 
be confirmed. Simulations have shown that the correction with the DMR scheme works better than the correction with RS alone. 


# Statement of need 

It is estimated that the amount of data will be up to 180 zettabytes in 2025, an increase 
of up to 23% from 2020 to 2025 [@IDC:2021]. Based on these forecasts, a new storage option is being sought, whereby the storage of data in DNA is being 
considered. DNA offers the advantages of extremely high information density and longevity, enabling space-saving and resource-conserving long-term storage. 
To store data in DNA, it must first be translated into DNA using a transcription scheme so that the encoded sequence can then be synthesized and sequenced.
During synthesis and sequencing, errors can occur, therefore they have to be corrected before decoding to ensure that an error-free output file can be 
restored. 

The popular error-correcting coding system is Reed Solomon (RS). It is used in a wide range of applications in digital communications and data storage and find
use in mass storage devices (hard disk drives, DVD, barcode tags), wireless and mobile communications units, satellite links, digital TV, digital video 
broadcasting (DVB), and modem technologies like xDSL [@Shrivastava:2013]. To increase the efficiency of the correction in data storage applications,
the RS code is used in combination with the own correction method. For example, Erlich and Zielinski [@Erlich:2017], Press et. al [@Press:2020] and 
Blawat et. al [@Blawat:2016] use RS within their developed correction methods.

For efficient decoding of data, it is necessary to determine an optimal correction method. This requires different methods to be developed, tested and 
compared with each other. For this reason, we developed the "dynamic mapping rule" (DMR) method, which provides coding and decoding including error 
correction of data. The method is based on a mapping table, which is used for the translation and also gives the dynamic mapping rule its name. 
With the help of this mapping table, a self-correction of the DNA sequence is achieved. For encoding, the input sequence was divided into segments, 
translated and reassembled using defined DNA sequences known as spacers. The spacers are placed between the individual segments so that the segments 
can be separated and corrected during decoding using the spacer sequences. The selection of the segment size is variable, but it is kept small 
so that the correction with DMR works better and the runtime is reduced. In addition, the coding 
was combined with RS so that the correction with the DMR scheme can be validated with the help of RS and thus leads to a better correction than the 
stand-alone use of DMR or RS. Because RS can only correct a certain number of errors, which depends on the encoded error correction symbols, it is not able 
to correct errors if the maximum error limit is exceeded. However, the use of DMR enables the correction of error overflows and ensures that the number of 
errors remains within the correction limits of RS, thus restoring RS functionality. Furthermore, the use of the DMR scheme offers the possibility to localize 
incorrect areas of the DNA by comparison with the mapping possibilities and thus carry out targeted corrective procedures at the affected sites. 
This should make it possible to correct substitution as well as insertion and deletion errors. 

The correction in the DMR scheme is executed with different designed levels, which proceed one after the other and try out different methods for the correction.
In the first level, for example, the incorrect sections are determined and an attempt is made to replace incorrect bases at precisely these points and thus 
repair the errors. To do this, the mapping table is consulted and all possible combinations for the previous and next 2-mers are searched for. Of these 
identified 2-mers, only those that show at least 1 matching base to the faulty 2-mere are tested. The possible sequences are checked using the RS code 
to correct them. This step determines whether the possibilities found are correct or not. The subsequent levels are designed in such a way that more and more 
possibilities are found if no correction has been made in the previous level. For the correction, the levels must be adapted and expanded accordingly so that 
an optimal correction can be achieved. \autoref{fig:level_design} shows the schematic illustration of the described DMR correction.     

The current code was developed internally as part of the BIOSYNTH project and is now being published as open software. The current designed levels are only 
suitable for the correction of substitution errors, but it is possible to extend them for the correction of insertion and deletion errors. 
The developed DMR method was compared with the RS method by encoding the small Fraunhofer logo with the segmented bit array encoding method using different 
amounts of error-correcting symbols (example in \autoref{fig:level_design}). Substitution errors were then added to the encoded DNA and decoded. 
A total of 20 trials for DMR and 50 trials for RS 
were performed per setting. \autoref{fig:error_correction} shows the results of this comparison and shows that the correction with the DMR scheme 
produces higher edit distances than the correction with RS, which proves the more effective correction of substitution errors. 


# Figures

![Schematic illustration of the dmr error correction. Firstly, the spacers are removed from the DNA sequence and the segments are initially scanned for correct 
and incorrect sequences. The incorrect segments then enter the DMR correction process. The levels that are passed through attempt to create segments that match 
the mapping scheme, these segments are then verified with RS. If the check with RS does not work, an attempt is made to correct the incorrect segment with the 
next level. If the check works, all corrected possible segments are analysed. If they are all the same or one segment has the majority, it is accepted as 
corrected. If none of the segments stand out by majority, the segment with the greatest similarity to the incorrect segment is assumed to be corrected. 
Created with BioRender.com \label{fig:level_design}](dmr_level_design.png)

![Comparison of DMR and RS correction with different settings and error rates of the bit array segmentation coding method. The simulation of the DMR correction
is shown on the left and was performed with 20 trials. The simulation with RS is shown on the right and was performed with 50 trials. In both diagrams, fully 
corrected images with an edit distance of 1 are shown in blue and erroneous images with a smaller edit distance are shown in red. The coding was carried out
with the help of a bit array segmentation. Here, the image is first read in and then the pixels are divided into segments. The individual segments are then 
compressed and translated. \label{fig:error_correction}](DMR_RS_correction_bitarray_segmentation.png)

# Acknowledgements

This work was funded by the Fraunhofer PREPARE Programm and developed within the BIOSYNTH project under the project number 40-03168-2420 as part of the data 
coding research.

# References