# Import pathways and generate a reference network

Utility to create pathway lists as igraph objects and interaction
networks from Reactome, KEGG, and other pathway databases.

## Usage

``` r
loadPathways(db, organism = "hsapiens", id_type = "ENTREZID", lcc = TRUE, ...)
```

## Arguments

- db:

  String indicating the database name. Please, check the
  [`pathways`](https://rdrr.io/pkg/graphite/man/pathways.html) function
  from graphite to list the available datasets.

- organism:

  A string indicating the source organism. Please, check the
  [`pathways`](https://rdrr.io/pkg/graphite/man/pathways.html) function
  from graphite to list the available datasets (default = "hsapiens")

- id_type:

  Gene ID type. The default is set to "ENTREZID" (standard SEM fitting
  nomenclature). A common choice could be "SYMBOL", for official gene
  symbols.

- lcc:

  A logical value. If TRUE (default), the reference network will only
  include the largest connected component. It will include all
  disconnected components otherwise.

- ...:

  Currently ignored.

## Value

A list of 2 objects:

1.  a list of pathways ad igraph objects;

2.  the union of graphs in the pathway list.

## Details

This function uses `graphite` to download and preprocess network data
from pathway databases. The output is then created using igraph and
SEMgraph utilities.

## References

Sales G, Calura E, Cavalieri D, Romualdi C (2012). graphite - a
Bioconductor package to convert pathway topology to gene network. BMC
Bioinformatics.
\<https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-20\>

## Author

Fernando Palluzzi <fernando.palluzzi@gmail.com>

## Examples

``` r
# \dontrun{

# Create KEGG reference pathway list and reference network for Homo sapiens
kegg.hs <- loadPathways("kegg", "hsapiens", "ENTREZID")
#> 
#> Delete pathway 35 of 332: D-Amino acid metabolism
#> Delete pathway 53 of 332: Linoleic acid metabolism
#> Delete pathway 66 of 332: Vitamin B6 metabolism
#> Delete pathway 69 of 332: Biotin metabolism
#> Delete pathway 75 of 332: Nitrogen metabolism
#> Delete pathway 82 of 332: Metabolic pathways
#> Delete pathway 153 of 332: Cornified envelope formation
#> Delete pathway 230 of 332: Type I diabetes mellitus
#> Delete pathway 242 of 332: Vitamin digestion and absorption
#> Delete pathway 328 of 332: Valine, leucine and isoleucine biosynthesis
#> Delete pathway 329 of 332: Neomycin, kanamycin and gentamicin biosynthesis
#> Delete pathway 330 of 332: Sulfur cycle
#> Delete pathway 331 of 332: Protein digestion and absorption
#> Delete pathway 332 of 332: Nicotine addiction
#> Frequency distribution of graph components
#> 
#>    n.nodes n.graphs
#> 1        2       29
#> 2        3        7
#> 3        4        7
#> 4        5        1
#> 5        6        6
#> 6        7        1
#> 7        8        2
#> 8       10        4
#> 9       15        1
#> 10      16        1
#> 11      20        1
#> 12    5962        1
#> 
#> Percent of vertices in the giant component: 92.3 %
#> 
#>   is.simple      is.dag is.directed is.weighted 
#>        TRUE       FALSE        TRUE       FALSE 
#> 
#> which.mutual.FALSE  which.mutual.TRUE 
#>              51339               1754 

# Inspect results
names(kegg.hs$pathways)
#>   [1] "Glycolysis / Gluconeogenesis"                                           
#>   [2] "Citrate cycle (TCA cycle)"                                              
#>   [3] "Pentose phosphate pathway"                                              
#>   [4] "Pentose and glucuronate interconversions"                               
#>   [5] "Fructose and mannose metabolism"                                        
#>   [6] "Galactose metabolism"                                                   
#>   [7] "Ascorbate and aldarate metabolism"                                      
#>   [8] "Fatty acid biosynthesis"                                                
#>   [9] "Fatty acid elongation"                                                  
#>  [10] "Fatty acid degradation"                                                 
#>  [11] "Steroid biosynthesis"                                                   
#>  [12] "Primary bile acid biosynthesis"                                         
#>  [13] "Ubiquinone and other terpenoid-quinone biosynthesis"                    
#>  [14] "Steroid hormone biosynthesis"                                           
#>  [15] "Oxidative phosphorylation"                                              
#>  [16] "Arginine biosynthesis"                                                  
#>  [17] "Purine metabolism"                                                      
#>  [18] "Caffeine metabolism"                                                    
#>  [19] "Pyrimidine metabolism"                                                  
#>  [20] "Alanine, aspartate and glutamate metabolism"                            
#>  [21] "Glycine, serine and threonine metabolism"                               
#>  [22] "Cysteine and methionine metabolism"                                     
#>  [23] "Valine, leucine and isoleucine degradation"                             
#>  [24] "Lysine degradation"                                                     
#>  [25] "Arginine and proline metabolism"                                        
#>  [26] "Histidine metabolism"                                                   
#>  [27] "Tyrosine metabolism"                                                    
#>  [28] "Phenylalanine metabolism"                                               
#>  [29] "Tryptophan metabolism"                                                  
#>  [30] "Phenylalanine, tyrosine and tryptophan biosynthesis"                    
#>  [31] "beta-Alanine metabolism"                                                
#>  [32] "Taurine and hypotaurine metabolism"                                     
#>  [33] "Phosphonate and phosphinate metabolism"                                 
#>  [34] "Selenocompound metabolism"                                              
#>  [35] "Glutathione metabolism"                                                 
#>  [36] "Starch and sucrose metabolism"                                          
#>  [37] "N-Glycan biosynthesis"                                                  
#>  [38] "Mucin type O-glycan biosynthesis"                                       
#>  [39] "Various types of N-glycan biosynthesis"                                 
#>  [40] "Mannose type O-glycan biosynthesis"                                     
#>  [41] "Amino sugar and nucleotide sugar metabolism"                            
#>  [42] "Glycosaminoglycan degradation"                                          
#>  [43] "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate"
#>  [44] "Glycosaminoglycan biosynthesis - heparan sulfate / heparin"             
#>  [45] "Biosynthesis of various nucleotide sugars"                              
#>  [46] "Glycerolipid metabolism"                                                
#>  [47] "Inositol phosphate metabolism"                                          
#>  [48] "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis"                 
#>  [49] "Glycerophospholipid metabolism"                                         
#>  [50] "Ether lipid metabolism"                                                 
#>  [51] "Arachidonic acid metabolism"                                            
#>  [52] "alpha-Linolenic acid metabolism"                                        
#>  [53] "Sphingolipid metabolism"                                                
#>  [54] "Glycosphingolipid biosynthesis - lacto and neolacto series"             
#>  [55] "Glycosphingolipid biosynthesis - globo and isoglobo series"             
#>  [56] "Glycosphingolipid biosynthesis - ganglio series"                        
#>  [57] "Pyruvate metabolism"                                                    
#>  [58] "Glyoxylate and dicarboxylate metabolism"                                
#>  [59] "Propanoate metabolism"                                                  
#>  [60] "Butanoate metabolism"                                                   
#>  [61] "One carbon pool by folate"                                              
#>  [62] "Thiamine metabolism"                                                    
#>  [63] "Riboflavin metabolism"                                                  
#>  [64] "Nicotinate and nicotinamide metabolism"                                 
#>  [65] "Pantothenate and CoA biosynthesis"                                      
#>  [66] "Lipoic acid metabolism"                                                 
#>  [67] "Folate biosynthesis"                                                    
#>  [68] "Retinol metabolism"                                                     
#>  [69] "Porphyrin metabolism"                                                   
#>  [70] "Terpenoid backbone biosynthesis"                                        
#>  [71] "Sulfur metabolism"                                                      
#>  [72] "Aminoacyl-tRNA biosynthesis"                                            
#>  [73] "Metabolism of xenobiotics by cytochrome P450"                           
#>  [74] "Drug metabolism - cytochrome P450"                                      
#>  [75] "Drug metabolism - other enzymes"                                        
#>  [76] "Biosynthesis of unsaturated fatty acids"                                
#>  [77] "Carbon metabolism"                                                      
#>  [78] "2-Oxocarboxylic acid metabolism"                                        
#>  [79] "Fatty acid metabolism"                                                  
#>  [80] "Biosynthesis of amino acids"                                            
#>  [81] "Nucleotide metabolism"                                                  
#>  [82] "Biosynthesis of cofactors"                                              
#>  [83] "Biosynthesis of nucleotide sugars"                                      
#>  [84] "EGFR tyrosine kinase inhibitor resistance"                              
#>  [85] "Endocrine resistance"                                                   
#>  [86] "Antifolate resistance"                                                  
#>  [87] "Platinum drug resistance"                                               
#>  [88] "mRNA surveillance pathway"                                              
#>  [89] "RNA degradation"                                                        
#>  [90] "Viral life cycle - HIV-1"                                               
#>  [91] "PPAR signaling pathway"                                                 
#>  [92] "Homologous recombination"                                               
#>  [93] "Fanconi anemia pathway"                                                 
#>  [94] "MAPK signaling pathway"                                                 
#>  [95] "ErbB signaling pathway"                                                 
#>  [96] "Ras signaling pathway"                                                  
#>  [97] "Rap1 signaling pathway"                                                 
#>  [98] "Calcium signaling pathway"                                              
#>  [99] "cGMP-PKG signaling pathway"                                             
#> [100] "cAMP signaling pathway"                                                 
#> [101] "Cytokine-cytokine receptor interaction"                                 
#> [102] "Viral protein interaction with cytokine and cytokine receptor"          
#> [103] "Chemokine signaling pathway"                                            
#> [104] "NF-kappa B signaling pathway"                                           
#> [105] "HIF-1 signaling pathway"                                                
#> [106] "FoxO signaling pathway"                                                 
#> [107] "Phosphatidylinositol signaling system"                                  
#> [108] "Sphingolipid signaling pathway"                                         
#> [109] "Phospholipase D signaling pathway"                                      
#> [110] "Neuroactive ligand-receptor interaction"                                
#> [111] "Hormone signaling"                                                      
#> [112] "Neuroactive ligand signaling"                                           
#> [113] "Cell cycle"                                                             
#> [114] "Oocyte meiosis"                                                         
#> [115] "p53 signaling pathway"                                                  
#> [116] "Sulfur relay system"                                                    
#> [117] "SNARE interactions in vesicular transport"                              
#> [118] "Autophagy - other"                                                      
#> [119] "Mitophagy - animal"                                                     
#> [120] "Autophagy - animal"                                                     
#> [121] "Protein processing in endoplasmic reticulum"                            
#> [122] "Endocytosis"                                                            
#> [123] "Phagosome"                                                              
#> [124] "Peroxisome"                                                             
#> [125] "Efferocytosis"                                                          
#> [126] "mTOR signaling pathway"                                                 
#> [127] "PI3K-Akt signaling pathway"                                             
#> [128] "AMPK signaling pathway"                                                 
#> [129] "Apoptosis"                                                              
#> [130] "Longevity regulating pathway"                                           
#> [131] "Longevity regulating pathway - multiple species"                        
#> [132] "Apoptosis - multiple species"                                           
#> [133] "Ferroptosis"                                                            
#> [134] "Necroptosis"                                                            
#> [135] "Cellular senescence"                                                    
#> [136] "Cardiac muscle contraction"                                             
#> [137] "Adrenergic signaling in cardiomyocytes"                                 
#> [138] "Vascular smooth muscle contraction"                                     
#> [139] "Wnt signaling pathway"                                                  
#> [140] "Notch signaling pathway"                                                
#> [141] "Hedgehog signaling pathway"                                             
#> [142] "TGF-beta signaling pathway"                                             
#> [143] "Axon guidance"                                                          
#> [144] "VEGF signaling pathway"                                                 
#> [145] "Apelin signaling pathway"                                               
#> [146] "Osteoclast differentiation"                                             
#> [147] "Hippo signaling pathway"                                                
#> [148] "Hippo signaling pathway - multiple species"                             
#> [149] "Focal adhesion"                                                         
#> [150] "ECM-receptor interaction"                                               
#> [151] "Cell adhesion molecules"                                                
#> [152] "Adherens junction"                                                      
#> [153] "Tight junction"                                                         
#> [154] "Gap junction"                                                           
#> [155] "Signaling pathways regulating pluripotency of stem cells"               
#> [156] "Complement and coagulation cascades"                                    
#> [157] "Platelet activation"                                                    
#> [158] "Antigen processing and presentation"                                    
#> [159] "Neutrophil extracellular trap formation"                                
#> [160] "Renin-angiotensin system"                                               
#> [161] "Toll-like receptor signaling pathway"                                   
#> [162] "NOD-like receptor signaling pathway"                                    
#> [163] "RIG-I-like receptor signaling pathway"                                  
#> [164] "Cytosolic DNA-sensing pathway"                                          
#> [165] "C-type lectin receptor signaling pathway"                               
#> [166] "JAK-STAT signaling pathway"                                             
#> [167] "Natural killer cell mediated cytotoxicity"                              
#> [168] "IL-17 signaling pathway"                                                
#> [169] "Th1 and Th2 cell differentiation"                                       
#> [170] "Th17 cell differentiation"                                              
#> [171] "T cell receptor signaling pathway"                                      
#> [172] "B cell receptor signaling pathway"                                      
#> [173] "Fc epsilon RI signaling pathway"                                        
#> [174] "Fc gamma R-mediated phagocytosis"                                       
#> [175] "TNF signaling pathway"                                                  
#> [176] "Leukocyte transendothelial migration"                                   
#> [177] "Intestinal immune network for IgA production"                           
#> [178] "Circadian rhythm"                                                       
#> [179] "Circadian entrainment"                                                  
#> [180] "Thermogenesis"                                                          
#> [181] "Long-term potentiation"                                                 
#> [182] "Synaptic vesicle cycle"                                                 
#> [183] "Neurotrophin signaling pathway"                                         
#> [184] "Retrograde endocannabinoid signaling"                                   
#> [185] "Glutamatergic synapse"                                                  
#> [186] "Cholinergic synapse"                                                    
#> [187] "Serotonergic synapse"                                                   
#> [188] "GABAergic synapse"                                                      
#> [189] "Dopaminergic synapse"                                                   
#> [190] "Long-term depression"                                                   
#> [191] "Olfactory transduction"                                                 
#> [192] "Taste transduction"                                                     
#> [193] "Phototransduction"                                                      
#> [194] "Inflammatory mediator regulation of TRP channels"                       
#> [195] "Regulation of actin cytoskeleton"                                       
#> [196] "Insulin signaling pathway"                                              
#> [197] "Insulin secretion"                                                      
#> [198] "GnRH signaling pathway"                                                 
#> [199] "Ovarian steroidogenesis"                                                
#> [200] "Progesterone-mediated oocyte maturation"                                
#> [201] "Estrogen signaling pathway"                                             
#> [202] "Melanogenesis"                                                          
#> [203] "Prolactin signaling pathway"                                            
#> [204] "Thyroid hormone synthesis"                                              
#> [205] "Thyroid hormone signaling pathway"                                      
#> [206] "Adipocytokine signaling pathway"                                        
#> [207] "Oxytocin signaling pathway"                                             
#> [208] "Glucagon signaling pathway"                                             
#> [209] "Regulation of lipolysis in adipocytes"                                  
#> [210] "Renin secretion"                                                        
#> [211] "Aldosterone synthesis and secretion"                                    
#> [212] "Relaxin signaling pathway"                                              
#> [213] "Cortisol synthesis and secretion"                                       
#> [214] "Parathyroid hormone synthesis, secretion and action"                    
#> [215] "GnRH secretion"                                                         
#> [216] "Type II diabetes mellitus"                                              
#> [217] "Insulin resistance"                                                     
#> [218] "Non-alcoholic fatty liver disease"                                      
#> [219] "AGE-RAGE signaling pathway in diabetic complications"                   
#> [220] "Cushing syndrome"                                                       
#> [221] "Growth hormone synthesis, secretion and action"                         
#> [222] "Alcoholic liver disease"                                                
#> [223] "Maturity onset diabetes of the young"                                   
#> [224] "Aldosterone-regulated sodium reabsorption"                              
#> [225] "Endocrine and other factor-regulated calcium reabsorption"              
#> [226] "Vasopressin-regulated water reabsorption"                               
#> [227] "Proximal tubule bicarbonate reclamation"                                
#> [228] "Salivary secretion"                                                     
#> [229] "Gastric acid secretion"                                                 
#> [230] "Pancreatic secretion"                                                   
#> [231] "Carbohydrate digestion and absorption"                                  
#> [232] "Fat digestion and absorption"                                           
#> [233] "Bile secretion"                                                         
#> [234] "Mineral absorption"                                                     
#> [235] "Cholesterol metabolism"                                                 
#> [236] "Cobalamin transport and metabolism"                                     
#> [237] "Alzheimer disease"                                                      
#> [238] "Parkinson disease"                                                      
#> [239] "Amyotrophic lateral sclerosis"                                          
#> [240] "Huntington disease"                                                     
#> [241] "Spinocerebellar ataxia"                                                 
#> [242] "Prion disease"                                                          
#> [243] "Pathways of neurodegeneration - multiple diseases"                      
#> [244] "Cocaine addiction"                                                      
#> [245] "Amphetamine addiction"                                                  
#> [246] "Morphine addiction"                                                     
#> [247] "Alcoholism"                                                             
#> [248] "Bacterial invasion of epithelial cells"                                 
#> [249] "Vibrio cholerae infection"                                              
#> [250] "Epithelial cell signaling in Helicobacter pylori infection"             
#> [251] "Pathogenic Escherichia coli infection"                                  
#> [252] "Shigellosis"                                                            
#> [253] "Salmonella infection"                                                   
#> [254] "Pertussis"                                                              
#> [255] "Legionellosis"                                                          
#> [256] "Yersinia infection"                                                     
#> [257] "Leishmaniasis"                                                          
#> [258] "Chagas disease"                                                         
#> [259] "African trypanosomiasis"                                                
#> [260] "Malaria"                                                                
#> [261] "Toxoplasmosis"                                                          
#> [262] "Amoebiasis"                                                             
#> [263] "Staphylococcus aureus infection"                                        
#> [264] "Tuberculosis"                                                           
#> [265] "Hepatitis C"                                                            
#> [266] "Hepatitis B"                                                            
#> [267] "Measles"                                                                
#> [268] "Human cytomegalovirus infection"                                        
#> [269] "Influenza A"                                                            
#> [270] "Human papillomavirus infection"                                         
#> [271] "Human T-cell leukemia virus 1 infection"                                
#> [272] "Kaposi sarcoma-associated herpesvirus infection"                        
#> [273] "Herpes simplex virus 1 infection"                                       
#> [274] "Epstein-Barr virus infection"                                           
#> [275] "Human immunodeficiency virus 1 infection"                               
#> [276] "Coronavirus disease - COVID-19"                                         
#> [277] "Pathways in cancer"                                                     
#> [278] "Transcriptional misregulation in cancer"                                
#> [279] "Viral carcinogenesis"                                                   
#> [280] "Chemical carcinogenesis - DNA adducts"                                  
#> [281] "Proteoglycans in cancer"                                                
#> [282] "MicroRNAs in cancer"                                                    
#> [283] "Chemical carcinogenesis - receptor activation"                          
#> [284] "Chemical carcinogenesis - reactive oxygen species"                      
#> [285] "Colorectal cancer"                                                      
#> [286] "Renal cell carcinoma"                                                   
#> [287] "Pancreatic cancer"                                                      
#> [288] "Endometrial cancer"                                                     
#> [289] "Glioma"                                                                 
#> [290] "Prostate cancer"                                                        
#> [291] "Thyroid cancer"                                                         
#> [292] "Basal cell carcinoma"                                                   
#> [293] "Melanoma"                                                               
#> [294] "Bladder cancer"                                                         
#> [295] "Chronic myeloid leukemia"                                               
#> [296] "Acute myeloid leukemia"                                                 
#> [297] "Small cell lung cancer"                                                 
#> [298] "Non-small cell lung cancer"                                             
#> [299] "Breast cancer"                                                          
#> [300] "Hepatocellular carcinoma"                                               
#> [301] "Gastric cancer"                                                         
#> [302] "Central carbon metabolism in cancer"                                    
#> [303] "Choline metabolism in cancer"                                           
#> [304] "PD-L1 expression and PD-1 checkpoint pathway in cancer"                 
#> [305] "Asthma"                                                                 
#> [306] "Autoimmune thyroid disease"                                             
#> [307] "Inflammatory bowel disease"                                             
#> [308] "Systemic lupus erythematosus"                                           
#> [309] "Rheumatoid arthritis"                                                   
#> [310] "Allograft rejection"                                                    
#> [311] "Graft-versus-host disease"                                              
#> [312] "Hypertrophic cardiomyopathy"                                            
#> [313] "Arrhythmogenic right ventricular cardiomyopathy"                        
#> [314] "Dilated cardiomyopathy"                                                 
#> [315] "Diabetic cardiomyopathy"                                                
#> [316] "Viral myocarditis"                                                      
#> [317] "Lipid and atherosclerosis"                                              
#> [318] "Fluid shear stress and atherosclerosis"                                 
kegg.hs$network
#> IGRAPH dba13c8 DN-- 5962 53093 -- 
#> + attr: name (v/c), label (v/c)
#> + edges from dba13c8 (vertex names):
#>  [1] 10327->9104      124  ->120356740 124  ->1312      124  ->220074   
#>  [5] 124  ->316       124  ->216       124  ->220       124  ->8854     
#>  [9] 124  ->339761    125  ->120356740 125  ->1312      125  ->220074   
#> [13] 125  ->316       125  ->216       125  ->220       125  ->8854     
#> [17] 125  ->339761    126  ->120356740 126  ->1312      126  ->220074   
#> [21] 126  ->316       126  ->216       126  ->220       126  ->8854     
#> [25] 126  ->339761    127  ->120356740 127  ->1312      127  ->220074   
#> [29] 127  ->316       127  ->216       127  ->220       127  ->8854     
#> + ... omitted several edges

# }
```
