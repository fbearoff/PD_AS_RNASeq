library(rentrez)

#100 most recent entries of term plus top N genes by padj
retrieve_PMIDs <- function(term_name, N){
    search_result <<- list()
    search_count <<- data.table()
    for (gene in z[order(padj)][abs(log2FoldChange)>2][1:N, gene_symbol]){
             r_search <- entrez_search(db="pubmed", term=paste0(term_name,"[TIAB] AND ", gene, "[TIAB]"),retmax=100)
             search_result <<- append(search_result, setNames(list(as.character(r_search$ids)), gene))
             search_count <<- rbind(search_count, data.table(Gene=gene, Count=r_search$count))
    }
    search_count <<- search_count[order(Gene)]
        pubmed_search_plot <<- ggplot(search_count,
                                  aes(x=Gene, y=Count)) +
                           geom_col() +
                           labs(title = paste0("Pubmed Entries for \"",
                                               term_name,
                                               "\" and Top ",
                                               N,
                                               " Differentlially Expressed Genes")
                               ) +
                           theme(axis.text.x=element_text(angle=45,
                                                          vjust=1,
                                                          hjust=1,
                                                          face=ifelse(search_count$Count>0,
                                                                      "bold",
                                                                      "plain")
                                                          )
                                )
    pubmed_search_plot
    print(search_result)
    print(search_count)
}

search_result[lapply(search_result,length)>0]
rlist::list.save(search_result[lapply(search_result,length)>0], 'Parkinson.yaml')



# r_fetch <- entrez_fetch(db="pubmed", id = r_search$ids, rettype="abstract")
