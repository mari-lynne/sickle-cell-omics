install.packages("rentrez")
library(rentrez)

genes <- c("FGR2A", "IFNG", "ABCC1")
term_of_interest <- "typhoid fever"

get_pubmed_links <- function(gene, term) {
  query <- paste(gene, term, sep=" AND ")
  ids <- entrez_search(db="pubmed", term=query)$ids
  
  if (length(ids) == 0) {
    return(data.frame(Gene=gene, Link=NA, stringsAsFactors=FALSE))
  } else {
    links <- lapply(ids, function(id) paste0("https://pubmed.ncbi.nlm.nih.gov/", id))
    return(data.frame(Gene=rep(gene, length(links)), Link=unlist(links), stringsAsFactors=FALSE))
  }
}

all_links <- lapply(genes, get_pubmed_links, term=term_of_interest)
all_links <- do.call(rbind, all_links)
print(all_links)



