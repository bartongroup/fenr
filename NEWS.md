## Version 0.1.5

 - added vignette
 - removeLazyData
 - DESCRIPTION updaates
 - replacing url_exists with RCurl::url.exists
 - added NEWS.md
 
## Version 0.1.6

 - bug fix in fetch_kegg
 - removed `min_count` and `fdr_limit` arguments from `functional_enrichment`; filtering can be done afterwards
 - added a small Shiny app as an example of `fenr`
 - updates to vignette

## Version 0.1.7

 - small fixes
 - link to a separate GitHub Shiny app added
 - added support for WikiPathways
 - improved robustness
 - more tests
 
 ## Version 0.1.8
 
 - `fetch_reactome` provides two ways of retrieving data, via one downloadable file or via APIs
 
 ## Version 0.1.9
 
  - Replacing ontologyIndex::get_ontology with a simpler parser
  - Replacing KEGGREST with own simple API parsers
  - Applying BiocStyle to the vignette
  - improving test_functional_enrichment
  
## Version 0.1.10

 - Style changes for BiocCheck
 - Adding more tests
 - Fixing a bug in `parse_kegg_genes`
 
 ## Version 0.1.11
 
 - Significant speed-up of enrichment by using Rfast::Hash in place of R lists
 - KEGG improvements, recognizing flat file genes with no gene synonym
 - Additional tests for Reactome
 - Minor improvements and fixes

 
