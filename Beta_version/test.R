### 05/03/2023
# Set working directory
setwd("/home/spu")

# Load library
library(tidyverse)

# Read in data
Paths_and_targets<-readxl::read_xlsx("Paths_and_targets.xlsx")
MOA_BOX<-readxl::read_xlsx("MOA_BOX.xlsx")
SJ_MOAset_full_20220906<-readxl::read_xlsx("SJ_MOAset_full_20220906.xlsx")

# Change Paths_and_targets to have one column containing all targets
Paths_and_targets<-Paths_and_targets %>%
  unite("targets",...6:...394,sep=",",na.rm=TRUE)

# Consolidate the gene_symbols, gene_ids, chembl_ids, pubchem_ids, dois, pubmed_ids, target_scores, selectivities, active_strengths columns of SJ_MOAset_full_20220906
SJ_MOAset_full_20220906<-SJ_MOAset_full_20220906 %>% 
  unite("gene_symbols",starts_with("gene_symbols"),sep=",",na.rm=TRUE) %>% 
  unite("gene_ids",starts_with("gene_ids"),sep=",",na.rm=TRUE) %>%
  unite("chembl_ids",starts_with("chembl_ids"),sep=",",na.rm=TRUE) %>%
  unite("pubchem_ids",starts_with("pubchem_ids"),sep=",",na.rm=TRUE) %>%
  unite("dois",starts_with("dois"),sep=",",na.rm=TRUE) %>%
  unite("pubmed_ids",starts_with("pubmed_ids"),sep=",",na.rm=TRUE) %>%
  unite("target_scores",starts_with("target_scores"),sep=",",na.rm=TRUE) %>% 
  unite("selectivities",starts_with("selectivities"),sep=",",na.rm=TRUE) %>%
  unite("active_strengths",starts_with("active_strengths"),sep=",",na.rm=TRUE)

# Save the string-type consolidated columns
save("Paths_and_targets","SJ_MOAset_full_20220906",file="Cleaned_SJ_Paths_and_targets_data_string.RData")

# Change the concatenated columns to list columns
## Paths_and_targets
Paths_and_targets<-Paths_and_targets %>%
  mutate(across(contains("targets"), ~ str_split(.x,pattern=",")))

## SJ_MOAset
concat_columns<-c("gene_symbols","gene_ids","chembl_ids","pubchem_ids","dois","pubmed_ids","target_scores","selectivities","active_strengths")
SJ_MOAset_full_20220906<-SJ_MOAset_full_20220906 %>%
  mutate(across(all_of(concat_columns), ~ str_split(.x,pattern=",")))
  
# Save RData (for list columns)
save("Paths_and_targets","SJ_MOAset_full_20220906",file="Cleaned_SJ_Paths_and_targets_data_list.RData")

### 05/04/2023
# Make the pathways a column (list of lists) in the SJ_MOAset_full_20220906 dataframe
Paths_and_targets_list<-Paths_and_targets %>% 
  pull('targets')
names(Paths_and_targets_list)<-Paths_and_targets %>% pull('Pathway')
Paths_and_targets_list_names<-names(Paths_and_targets_list)
Paths_and_targets_flat<-unlist(Paths_and_targets_list) %>% unique()

SJ_MOAset_full_20220906_targets<-SJ_MOAset_full_20220906 %>%
  pull('gene_symbols') 
SJ_MOAset_full_20220906_targets_flat<-SJ_MOAset_full_20220906_targets %>%
  unlist() %>%
  unique()
SJ_MOAset_full_20220906_target_names<-SJ_MOAset_full_20220906 %>% pull('SJN')

SJ_MOAset_full_20220906_targets_flat<-unlist(SJ_MOAset_full_20220906_targets) %>% unique()
SJ_MOAset_paths<-map(set_names(SJ_MOAset_full_20220906_targets_flat), ~list())

# save flat data
save(Paths_and_targets_flat,SJ_MOAset_full_20220906_targets_flat,file="flat_data.RData")

node_types<-data.frame(node=SJ_MOAset_full_20220906_targets_flat,type="gene")
node_types<-rbind(node_types,data.frame(node=Paths_and_targets_list_names,type="pathway"))
# save node types
saveRDS(node_types,file="node_types.rds")

node_types_full<-rbind(node_types,data.frame(node=names(SJ_CTP_list),type="compound"))
node_types_full<-rbind(node_types_full,data.frame(node=setdiff(Paths_and_targets_flat,SJ_MOAset_full_20220906_targets_flat),type="gene_no_compound"))
saveRDS(node_types_full,file="node_types_full.rds")

# for each pathway, if it has a gene in it, add that pathway name to the gene in SJ_MOAset_paths
counter<-0
for(i in Paths_and_targets_list){
  counter<-counter+1
  for(j in i){
    SJ_MOAset_paths[[j]]<-append(SJ_MOAset_paths[[j]],Paths_and_targets_list_names[counter])
  }
}

SJ_MOA_genes<-SJ_MOAset_full_20220906$gene_symbols
pathways_col<-vector(mode="list",length=length(SJ_MOA_genes))

for(i in 1:length(SJ_MOA_genes)){
  pathways_col[[i]]<-SJ_MOAset_paths[intersect(SJ_MOA_genes[[i]],names(SJ_MOAset_paths))]
}

# Add pathways to the SJ_MOAset_full_20220906
SJ_MOAset_full_20220906_paths<-SJ_MOAset_full_20220906 %>% 
  mutate(pathways = pathways_col)

# Save SJ_MOAset with pathways added
save(SJ_MOAset_full_20220906_paths,file="Cleaned_SJ_MOAset_paths.RData")

# Make genes and compounds list
genes_and_compounds_list<-vector(mode="list",length=length(SJ_MOAset_full_20220906_targets_flat))
names(genes_and_compounds_list)<-SJ_MOAset_full_20220906_targets_flat

for(i in SJ_MOAset_full_20220906_targets_flat){
  for(j in 1:nrow(SJ_MOAset_full_20220906)){
    if(i %in% (SJ_MOAset_full_20220906$gene_symbols[[j]] %>% unlist())){
      genes_and_compounds_list[[i]]<-append(genes_and_compounds_list[[i]],SJ_MOAset_full_20220906$SJN[j])
    }
  }
}

# save genes_and_compounds
saveRDS(genes_and_compounds_list,file="genes_and_compounds_list.rds")

compound_gene_strengths<-data.frame(SJN=character(),genes=character(),active_strengths=numeric())
for(i in 1:nrow(SJ_MOAset_full_20220906)){
  genes<-SJ_MOAset_full_20220906$gene_symbols[i][[1]]
  active_strengths<-SJ_MOAset_full_20220906$active_strengths[i][[1]]
  SJN<-SJ_MOAset_full_20220906$SJN[i]
  for(j in 1:length(genes)){
    compound_gene_strengths<-rbind(compound_gene_strengths,data.frame(SJN=SJN,genes=genes[j],active_strengths=active_strengths[j]))
  }
}

# save compound_gene_strengths
saveRDS(compound_gene_strengths,file="compound_gene_strengths.rds")

# save each gene as a row, with compound:strength as list
compound_gene_strengths<-compound_gene_strengths %>% 
  unite("compound_strengths",all_of(c("SJN","active_strengths")),sep=" : ",na.rm=TRUE)
gene_compound_strengths<-data.frame(gene=SJ_MOAset_full_20220906_targets_flat,compound_strengths=rep(NA,length=length(SJ_MOAset_full_20220906_targets_flat)))
for(i in gene_compound_strengths$gene){
  matches<-compound_gene_strengths %>% filter(genes==i)
  strengths<-matches %>% pull('compound_strengths') %>% paste(collapse=" , ")
  gene_compound_strengths[gene_compound_strengths$gene==i,c("compound_strengths")]<-strengths
}

# save gene-compound strengths
saveRDS(gene_compound_strengths,file="gene_compound_strengths.rds")

## 
# Make S4 object for Compound
setClass("Compound",representation(SJNB="numeric",
                                   SJN="numeric",
                                   BATCHNUMBER="numeric",
                                   moabox_id="character",
                                   inchi_key="character",
                                   smiles="character",
                                   cas_number="character",
                                   moa="character",
                                   gene_symbols="list",
                                   gene_ids="list",
                                   chembl_ids="list",
                                   pubchem_ids="list",
                                   member_status="character",
                                   dois="list",
                                   pubmed_ids="character",
                                   target_scores="list",
                                   selectivities="list",
                                   active_strengths="list",
                                   UNIQUEID="character",
                                   VENDORNAME="character",
                                   BATCHCOMMENT="character",
                                   MOLSMILES="character",
                                   CAS="character",
                                   LOT="numeric",
                                   SALT_NAME="character",
                                   SALT_EQUIVALENTS="numeric",
                                   SOLVATE_NAME="character",
                                   SOLVATE_EQUIVALENTS="numeric",
                                   BATCH_INTERNAL_ID="numeric",
                                   MolecularWeight="numeric",
                                   SynonymPrimary="character",
                                   Synonym="character",
                                   MCE_Pathway="character",
                                   MCE_Research_Area="character",
                                   MOABOX_ID="character",
                                   Primary_Target="character",
                                   MOA="character"
                                   ))

# Make S4 class for PathsTargets
setClass("PathsTargets",representation(Pathway="character",
                                       Description="character",
                                       targets="list"
                                       ))

PathsTargets_list<-list()
for(i in 1:nrow(Paths_and_targets)){
  TP<-new("PathsTargets",
          Pathway=Paths_and_targets$Pathway[i],
          Description=Paths_and_targets$Description[i],
          targets=Paths_and_targets$targets[i])
  PathsTargets_list<-append(PathsTargets_list,TP)
}
names(PathsTargets_list)<-Paths_and_targets$Pathway

# save PathsTargets list
save(PathsTargets_list,file="PathsTargets_list.RData")
saveRDS(PathsTargets_list,file="PathsTargets_list.rds")

# Make S4 class for CompoundTargetsPaths
setClass("CompoundTargetsPaths", slots = c(
  compound_description = "list",
  genes = "list",
  pathways="list",
  chemical_id = "list",
  MOA_info = "list"
))

CompoundTargetsPaths_StandardFormat<-new("CompoundTargetsPaths",
              compound_description = list(
                SJNB = "character",
                SJN = "character",
                BATCHNUMBER = "numeric",
                moabox_id = "character",
                inchi_key = "character",
                smiles = "character",
                cas_number = "character",
                moa = "character"
              ),
              genes = list(
                gene_symbols = list(),
                gene_ids = list(),
                target_scores = list(),
                selectivities = list(),
                active_strengths = list()
              ),
              pathways = list(
                pathways = list()
              ),
              chemical_id = list(
                chembl_id = list()
              ),
              MOA_info = list(
                UNIQUEID = "character",
                VENDORNAME = "character",
                BATCHCOMMENT = "character",
                MOLSMILES = "character",
                CAS = "character",
                LOT = "numeric",
                BATCHPROJECT = "character",
                SALT_NAME = "character",
                SALT_EQUIVALENTS = "numeric",
                SOLVATE_NAME = "character",
                SOLVATE_EQUIVALENTS = "numeric",
                BATCH_INTERNAL_ID = "numeric",
                MolecularWeight = "numeric",
                SynonymPrimary = "character",
                Synonym = "character",
                MCE_Pathway = "character",
                MCE_Research_Area = "character",
                MOABOX_ID = "character",
                Primary_Target = "character",
                MOA = "character"
              ))

test_CompoundTargetsPaths<-new("CompoundTargetsPaths",
                               compound_description = list(
                                 SJNB = SJ_MOAset_full_20220906_paths$SJNB[1],
                                 SJN = SJ_MOAset_full_20220906_paths$SJN[1],
                                 BATCHNUMBER = SJ_MOAset_full_20220906_paths$BATCHNUMBER[1],
                                 moabox_id = SJ_MOAset_full_20220906_paths$moabox_id[1],
                                 inchi_key = SJ_MOAset_full_20220906_paths$inchi_key[1],
                                 smiles = SJ_MOAset_full_20220906_paths$smiles[1],
                                 cas_number = SJ_MOAset_full_20220906_paths$cas_number[1],
                                 moa = SJ_MOAset_full_20220906_paths$moa[1]
                               ),
                               genes = list(
                                 gene_symbols = SJ_MOAset_full_20220906_paths$gene_symbols[1],
                                 gene_ids = SJ_MOAset_full_20220906_paths$gene_ids[1],
                                 target_scores = SJ_MOAset_full_20220906_paths$target_scores[1],
                                 selectivities = SJ_MOAset_full_20220906_paths$selectivities[1],
                                 active_strengths = SJ_MOAset_full_20220906_paths$active_strengths[1]
                               ),
                               pathways = list(
                                 pathways = SJ_MOAset_full_20220906_paths$pathways[1][[1]]
                               ),
                               chemical_id = list(
                                 chembl_ids = SJ_MOAset_full_20220906_paths$chembl_ids[1]
                               ),
                               MOA_info = list(
                                 UNIQUEID = SJ_MOAset_full_20220906_paths$UNIQUEID[1],
                                 VENDORNAME = SJ_MOAset_full_20220906_paths$VENDORNAME[1],
                                 BATCHCOMMENT = SJ_MOAset_full_20220906_paths$BATCHCOMMENT[1],
                                 MOLSMILES = SJ_MOAset_full_20220906_paths$MOLSMILES[1],
                                 CAS = SJ_MOAset_full_20220906_paths$CAS[1],
                                 LOT = SJ_MOAset_full_20220906_paths$LOT[1],
                                 BATCHPROJECT = SJ_MOAset_full_20220906_paths$BATCHPROJECT[1],
                                 SALT_NAME = SJ_MOAset_full_20220906_paths$SALT_NAME[1],
                                 SALT_EQUIVALENTS = SJ_MOAset_full_20220906_paths$SALT_EQUIVALENTS[1],
                                 SOLVATE_NAME = SJ_MOAset_full_20220906_paths$SOLVATE_NAME[1],
                                 SOLVATE_EQUIVALENTS = SJ_MOAset_full_20220906_paths$SOLVATE_EQUIVALENTS[1],
                                 BATCH_INTERNAL_ID = SJ_MOAset_full_20220906_paths$BATCH_INTERNAL_ID[1],
                                 MolecularWeight = SJ_MOAset_full_20220906_paths$MolecularWeight[1],
                                 SynonymPrimary = SJ_MOAset_full_20220906_paths$SynonymPrimary[1],
                                 Synonym = SJ_MOAset_full_20220906_paths$Synonym[1],
                                 MCE_Pathway = SJ_MOAset_full_20220906_paths$MCE_Pathway[1],
                                 MCE_Research_Area = SJ_MOAset_full_20220906_paths$MCE_Research_Area[1],
                                 MOABOX_ID = SJ_MOAset_full_20220906_paths$MOABOX_ID[1],
                                 Primary_Target = SJ_MOAset_full_20220906_paths$Primary_Target[1],
                                 MOA = SJ_MOAset_full_20220906_paths$MOA[1]
                               ))

# loop for making S4 CompoundTargetsPaths objects - full SJ data
SJ_CTP_list<-list()
for(i in 1:nrow(SJ_MOAset_full_20220906_paths)){
  SJ_CTP<-new("CompoundTargetsPaths",
           compound_description = list(
             SJNB = SJ_MOAset_full_20220906_paths$SJNB[i],
             SJN = SJ_MOAset_full_20220906_paths$SJN[i],
             BATCHNUMBER = SJ_MOAset_full_20220906_paths$BATCHNUMBER[i],
             moabox_id = SJ_MOAset_full_20220906_paths$moabox_id[i],
             inchi_key = SJ_MOAset_full_20220906_paths$inchi_key[i],
             smiles = SJ_MOAset_full_20220906_paths$smiles[i],
             cas_number = SJ_MOAset_full_20220906_paths$cas_number[i],
             moa = SJ_MOAset_full_20220906_paths$moa[i]
           ),
           genes = list(
             gene_symbols = SJ_MOAset_full_20220906_paths$gene_symbols[i],
             gene_ids = SJ_MOAset_full_20220906_paths$gene_ids[i],
             target_scores = SJ_MOAset_full_20220906_paths$target_scores[i],
             selectivities = SJ_MOAset_full_20220906_paths$selectivities[i],
             active_strengths = SJ_MOAset_full_20220906_paths$active_strengths[i]
           ),
           pathways = list(
             pathways = SJ_MOAset_full_20220906_paths$pathways[i]
           ),
           chemical_id = list(
             chembl_ids = SJ_MOAset_full_20220906_paths$chembl_ids[i]
           ),
           MOA_info = list(
             UNIQUEID = SJ_MOAset_full_20220906_paths$UNIQUEID[i],
             VENDORNAME = SJ_MOAset_full_20220906_paths$VENDORNAME[i],
             BATCHCOMMENT = SJ_MOAset_full_20220906_paths$BATCHCOMMENT[i],
             MOLSMILES = SJ_MOAset_full_20220906_paths$MOLSMILES[i],
             CAS = SJ_MOAset_full_20220906_paths$CAS[i],
             LOT = SJ_MOAset_full_20220906_paths$LOT[i],
             BATCHPROJECT = SJ_MOAset_full_20220906_paths$BATCHPROJECT[i],
             SALT_NAME = SJ_MOAset_full_20220906_paths$SALT_NAME[i],
             SALT_EQUIVALENTS = SJ_MOAset_full_20220906_paths$SALT_EQUIVALENTS[i],
             SOLVATE_NAME = SJ_MOAset_full_20220906_paths$SOLVATE_NAME[i],
             SOLVATE_EQUIVALENTS = SJ_MOAset_full_20220906_paths$SOLVATE_EQUIVALENTS[i],
             BATCH_INTERNAL_ID = SJ_MOAset_full_20220906_paths$BATCH_INTERNAL_ID[i],
             MolecularWeight = SJ_MOAset_full_20220906_paths$MolecularWeight[i],
             SynonymPrimary = SJ_MOAset_full_20220906_paths$SynonymPrimary[i],
             Synonym = SJ_MOAset_full_20220906_paths$Synonym[i],
             MCE_Pathway = SJ_MOAset_full_20220906_paths$MCE_Pathway[i],
             MCE_Research_Area = SJ_MOAset_full_20220906_paths$MCE_Research_Area[i],
             MOABOX_ID = SJ_MOAset_full_20220906_paths$MOABOX_ID[i],
             Primary_Target = SJ_MOAset_full_20220906_paths$Primary_Target[i],
             MOA = SJ_MOAset_full_20220906_paths$MOA[i]
           ))
  SJ_CTP_list<-append(SJ_CTP_list,SJ_CTP)
}
names(SJ_CTP_list)<-SJ_MOAset_full_20220906_paths$SJN

# Save full SJ CTP list
save(SJ_CTP_list,file="SJ_CTP_list.RData")
saveRDS(SJ_CTP_list,file="SJ_CTP_list.rds")

## chek that genes and active strengths are matched
gene_strength_count<-data.frame()
for(i in SJ_CTP_list){
  gene_count<-i@genes$gene_symbols %>% unlist() %>% length()
  strength_count<-i@genes$active_strengths %>% unlist() %>% length()
  gene_strength_count<-rbind(gene_strength_count,data.frame(genes=gene_count,strengths=strength_count))
}
gene_strength_count$is_equal<-gene_strength_count$genes==gene_strength_count$strengths

# loop for making smaller S4 CompoundTargetsPaths object for testing
CTP_list<-list()
for(i in 1:20){
  CTP<-new("CompoundTargetsPaths",
                compound_description = list(
                  SJNB = SJ_MOAset_full_20220906_paths$SJNB[i],
                  SJN = SJ_MOAset_full_20220906_paths$SJN[i],
                  BATCHNUMBER = SJ_MOAset_full_20220906_paths$BATCHNUMBER[i],
                  moabox_id = SJ_MOAset_full_20220906_paths$moabox_id[i],
                  inchi_key = SJ_MOAset_full_20220906_paths$inchi_key[i],
                  smiles = SJ_MOAset_full_20220906_paths$moabox_id[i],
                  cas_number = SJ_MOAset_full_20220906_paths$cas_number[i],
                  moa = SJ_MOAset_full_20220906_paths$moa[i]
                ),
                genes = list(
                  gene_symbols = SJ_MOAset_full_20220906_paths$gene_symbols[i],
                  gene_ids = SJ_MOAset_full_20220906_paths$gene_ids[i],
                  target_scores = SJ_MOAset_full_20220906_paths$target_scores[i],
                  selectivities = SJ_MOAset_full_20220906_paths$selectivities[i],
                  active_strengths = SJ_MOAset_full_20220906_paths$active_strengths[i]
                ),
                pathways = list(
                  pathways = SJ_MOAset_full_20220906_paths$pathways[i]
                ),
                chemical_id = list(
                  chembl_ids = SJ_MOAset_full_20220906_paths$chembl_ids[i]
                ),
                MOA_info = list(
                  UNIQUEID = SJ_MOAset_full_20220906_paths$UNIQUEID[i],
                  VENDORNAME = SJ_MOAset_full_20220906_paths$VENDORNAME[i],
                  BATCHCOMMENT = SJ_MOAset_full_20220906_paths$BATCHCOMMENT[i],
                  MOLSMILES = SJ_MOAset_full_20220906_paths$MOLSMILES[i],
                  CAS = SJ_MOAset_full_20220906_paths$CAS[i],
                  LOT = SJ_MOAset_full_20220906_paths$LOT[i],
                  BATCHPROJECT = SJ_MOAset_full_20220906_paths$BATCHPROJECT[i],
                  SALT_NAME = SJ_MOAset_full_20220906_paths$SALT_NAME[i],
                  SALT_EQUIVALENTS = SJ_MOAset_full_20220906_paths$SALT_EQUIVALENTS[i],
                  SOLVATE_NAME = SJ_MOAset_full_20220906_paths$SOLVATE_NAME[i],
                  SOLVATE_EQUIVALENTS = SJ_MOAset_full_20220906_paths$moabox_id[i],
                  BATCH_INTERNAL_ID = SJ_MOAset_full_20220906_paths$BATCH_INTERNAL_ID[i],
                  MolecularWeight = SJ_MOAset_full_20220906_paths$MolecularWeight[i],
                  SynonymPrimary = SJ_MOAset_full_20220906_paths$SynonymPrimary[i],
                  Synonym = SJ_MOAset_full_20220906_paths$Synonym[i],
                  MCE_Pathway = SJ_MOAset_full_20220906_paths$MCE_Pathway[i],
                  MCE_Research_Area = SJ_MOAset_full_20220906_paths$MCE_Research_Area[i],
                  MOABOX_ID = SJ_MOAset_full_20220906_paths$MOABOX_ID[i],
                  Primary_Target = SJ_MOAset_full_20220906_paths$Primary_Target[i],
                  MOA = SJ_MOAset_full_20220906_paths$MOA[i]
                ))
  CTP_list<-append(CTP_list,CTP)
}

names(CTP_list)<-SJ_MOAset_full_20220906_paths$SJN[1:20]

# save CTP_list
save(CTP_list,file="CTP_list.RData")

## Make a network object from a CompoundTargetsPaths
setClass("CompoundTargetsPathsNetwork", slots = c(
  nodes = "list",
  edges = "list"
))

# Make a sample CompoundTargetsPathsNetwork
CTPN_list<-list()
for(i in 1:20){
  CTPN<-new("CompoundTargetsPathsNetwork",
            nodes=list(genes=CTP_list[[i]]@genes$gene_symbols,
            pathways=CTP_list[[i]]@pathways),
            edges=list(edge_weights=CTP_list[[i]]@genes$selectivities))
  CTPN_list<-append(CTPN,CTPN_list)
}

# Make edgelist from CompoundTargetsPaths 
get_CompoundTargetsPaths_edgelist<-function(CTP_object){
  CTP_SJN<-CTP_object@compound_description$SJN
  gene_paths<-CTP_object@pathways$pathways[[1]]
  edge_list<-data.frame()
  counter<-0
  for(i in gene_paths){
    counter<-counter+1
    if(length(i)>0){
      gene_name<-names(gene_paths)[counter]
      paths<-unlist(i)
      curr_edges<-data.frame(edge1=rep(gene_name,length(paths)),edge2=paths)
      edge_list<-rbind(edge_list,curr_edges)
    }
  }
  return(edge_list)
}

CTP_list_edgelist<-data.frame()

for(i in CTP_list){
  CTP_list_edgelist<-rbind(CTP_list_edgelist,get_CompoundTargetsPaths_edgelist(i))
}

# Save CTP_list_edgelist
save(CTP_list_edgelist,file="CTP_list_edgelist.RData")

## Take input data and turn into network
SJNs_in<-names(CTP_list)[1:20]
SJ_CTP_list[SJNs_in]
# generate edgelists
for(i in 1:length(SJ_CTP_list[SJNs_in])){
  
}

# igraph
library(igraph)
CTP_list_edgelist_noSJN<-CTP_list_edgelist[,c("edge1","edge2")]
ig_object<-graph_from_edgelist(as.matrix(CTP_list_edgelist),directed=FALSE)
ig_object_V<-V(ig_object) %>% names()
ig_object_V_types<-dplyr::left_join(data.frame(node=ig_object_V),node_types) %>% pull('type')
ig_object_V_colors<-case_when(ig_object_V_types == 'gene' ~ '#8DD3C7', ig_object_V_types == 'pathway' ~ '#FFFFB3')

ig_object<-ig_object %>% 
  set_vertex_attr("node_type",value=ig_object_V_types) %>%
  set_vertex_attr("color",value=ig_object_V_colors)
plot(ig_object,edge.arrow.size=4, vertex.size=20,
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.4, vertex.label.dist=4, edge.curved=10)

# visNetwork
require(visNetwork, quietly = TRUE)
# minimal example
nodes <- data.frame(id = ig_object_V,
                    label = ig_object_V,
                    group = ig_object_V_types,
                    title = ig_object_V,
                    color = ig_object_V_colors)
edges <- data.frame(from = CTP_list_edgelist$edge1, to = CTP_list_edgelist$edge2)
CTP_network_test<-visNetwork(nodes, edges, width = "100%") %>%
  visGroups(groupname = "gene", color = "#8DD3C7") %>%
  visGroups(groupname = "pathway", color = "#FFFFB3") %>%
  visLegend(width = 0.1, position = "right", main = "group")

# nodes value will affect node size (1:3?)
# edges label, length, width, title, font.color, font.size
# visEdges
visEdges(width=strengths)

CTP_network_test %>% visSave(file = "CTP_network_test.html")


setwd("/home")