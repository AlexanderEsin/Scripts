library(CHNOSZ)
taxdir<-("/users/aesin/desktop/Clean_Proteomes/taxdmp")
nodes=getnodes(taxdir)
names=getnames(taxdir)
binomials<-scan(file="/users/aesin/desktop/Clean_proteomes/Taxa_names.tsv",what="characters",sep="\n")

getid <- function (name, taxdir, names = NULL)
{
  if (is.null(names))
    names <- getnames(taxdir)
  return(as.numeric(names$id[match(name,names$name)]))
}

species_up <- function (id, taxdir, nodes = NULL) 
{
  if (is.null(nodes)) 
    nodes <- getnodes(taxdir)
  thisid <- id
  current_rank <- getrank(thisid,nodes=nodes)
  if (current_rank != "species") {
    id <- parent(thisid, rank="species",nodes=nodes)
    thisid <- id
    while (thisid != 1) {
      thisid <- parent(thisid, taxdir = taxdir, nodes = nodes)
      id <- c(id, thisid)
    }
  }
  else {
    while (thisid != 1) {
    thisid <- parent(thisid, taxdir = taxdir, nodes = nodes)
    id <- c(id, thisid)
    }
  }
  return(id)
}
print("Hello")

z=1
alltaxa<-character()
for (x in binomials) {
  taxa<-sciname(species_up(getid(x,names=names),nodes=nodes),names=names)
  taxa<-paste(c(x,taxa),collapse="\t")
  alltaxa <- c(alltaxa,taxa)
  print(z)
  z=z+1
}

write(alltaxa,file="/users/aesin/desktop/Clean_proteomes/Taxonomy_all.tsv",sep="\n")

