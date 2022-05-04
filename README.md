# FUSED
FUSion Error prone repair Detection

<img src="Microhomology_Error_Prone_Repair_Spectrum.svg" alt="error prone repair cartoon" height="600" align="center"/>  

<br></br>
<br></br>  

## Description

FUSED is a collection of functions to determine the presence of hallmarks of common DNA repair mechanisms including NHEJ, MMEJ, SSA, BIR, and MMBIR. 

### NHEJ

For non-homologous end-joining, FUSED looks for homology between 1-4nt (nt==nucleotide) in length between the left fusion gene and the right fusion gene in a window up to 4nt away from the breakpoint. It requires the homology to be anchored *at* the breakpoint.

### SSA

For single-strand annealing, FUSED looks for homology between 30 and 70nt in length in a window around the breakpoint up to 70nt away on both sides of the breakpoint.

### MMEJ

For microhomology-mediated end joining, FUSED looks for homology between 1-20nt in length between the left fusion gene and the right fusion gene in a window up to 50nt away from the breakpoint. 

### BIR

For break-induced replication, FUSED looks for homology between 70-100nt in length in a window 100nt upstream and downstream of the breakpoint.

### MMBIR

For microhomology-mediated break-induced replication, FUSED looks for homology between 1-20nt in length in a window 20nt upstream and downstream of the breakpoint. It requires the homology to be anchored *at* the breakpoint, expects a gapped alignment of at least a single nt but with at least 80% overall homology.

## Input

FUSED expects a fusion prediction file from either: STAR, ARRIBA, or ORIEN derived pipelines. It utilizes a BSgenome object (or string) to pull matching sequence information from the coordinates of the predicted fusion. It can evaluate 5 different possible error repair mechanisms: NHEJ, MMEJ, SSA, BIR, or MMBIR.

## Usage

```R
devtools::install_github("databio/FUSED")
library(FUSED)
inputDT <- fread("ORIEN-FusionBreakpoints.txt")
calcHomology(inputDT, mechanism="NHEJ", bsg="BSgenome.Hsapiens.UCSC.hg38")
```
