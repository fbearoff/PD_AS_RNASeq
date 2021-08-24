library(rentrez)

#100 most recent entries of term plus top N genes by padj
retrieve_PMIDs <- function(term_name, N){
    search_result <<- list()
    for (gene in z[order(padj)][1:N, gene_symbol]){
             r_search <- entrez_search(db="pubmed", term=paste0(term_name,"[TIAB] AND ", gene, "[TIAB]"),retmax=100)
             search_result <<- append(search_result, setNames(list(as.character(r_search$ids)), gene))
             # search_count <<- rbind(search_count, gene = r_search$count)
    }
    return(search_result)
}

r_fetch <- entrez_fetch(db="pubmed", id = r_search$ids, rettype="abstract")
