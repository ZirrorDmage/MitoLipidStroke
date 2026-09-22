
rm(list=ls());gc()
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})
conflicted::conflict_prefer_matching("^filter$|^select$|^arrange$|^mutate$|", "dplyr",quiet = T)

####---- 项目根目录定位（脚本位于 <项目根>/Rscript/，脚本内全部使用相对路径）----####
.script_dir <- function() {
  a <- commandArgs(FALSE)
  i <- grep("^--file=", a)
  if (length(i) > 0) return(dirname(normalizePath(sub("^--file=", "", a[i[1]]), mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  for (k in seq_len(sys.nframe())) {
    f <- tryCatch(sys.frame(k)$ofile, error = function(e) NULL)
    if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f, mustWork = FALSE)))
  }
  if (basename(getwd()) == "Rscript") getwd() else file.path(getwd(), "Rscript")
}
SCRIPT_DIR <- .script_dir()
ROOT_CAND <- unique(c(SCRIPT_DIR, dirname(SCRIPT_DIR), getwd(), dirname(getwd())))
ROOT_HIT <- ROOT_CAND[file.exists(file.path(ROOT_CAND, "00.rawdata"))]
root <- if (nzchar(Sys.getenv("PROJECT_ROOT"))) {
  normalizePath(Sys.getenv("PROJECT_ROOT"), mustWork = FALSE)
} else if (length(ROOT_HIT) > 0) {
  normalizePath(ROOT_HIT[1], mustWork = FALSE)
} else if (basename(SCRIPT_DIR) == "Rscript") {
  normalizePath(dirname(SCRIPT_DIR), mustWork = FALSE)
} else {
  normalizePath(SCRIPT_DIR, mustWork = FALSE)
}
options(project.root = root)   # 供脚本中 rm(list = ls()) 之后继续使用
cat("项目根目录：", root, "\n")

setwd(root)

dir <- "11.network"
if(!dir.exists(dir)) { dir.create(dir)}
setwd(dir)



######miRNA---------
library(multiMiR)
hubgene <- read.csv('../06.expressionROC/keygenes.csv')[,1]
gene2mir <- get_multimir(org     = 'hsa',
                         target  = hubgene,
                         table   = 'predicted',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff= 500000)
table(gene2mir@data$database)
table(gene2mir@data$database, gene2mir@data$target_symbol)
saveRDS(gene2mir,"gene2mir.rds")

gene_sum <- gene2mir@summary  
table(gene_sum$target_symbol)
database_result<- gene_sum[6:(ncol(gene_sum)-2)]

get_count<- function(df,num){
  all_result <- c()
  combn_result<- combn(colnames(database_result),num)
  for (i in 1:ncol(combn_result)){
    tmp<- database_result[combn_result[,i]]
    
    tmp_string<- paste0("subset(tmp,",paste0(combn_result[,i]," !=0",collapse = " & "),")")
    nrow_com<- nrow(eval(parse(text = tmp_string)))
    tmp_name<- paste0(combn_result[,i],collapse = " & ")
    names(nrow_com) <- tmp_name
    all_result<- c(all_result,nrow_com)
  }
  return(sort(all_result,decreasing = T))
}
###
get_count(gene2mir,2)


mirdb = gene2mir@data[gene2mir@data$database=="mirdb",]
mirdb$miRNA <- paste0(mirdb$mature_mirna_id,"_",mirdb$target_symbol)
mirdb2<- rstatix::filter(mirdb,!duplicated(mirdb$miRNA))
table(mirdb2$target_symbol)
# CYP1B1   DEGS1   HSDL2 OSBPL1A    PTEN 
# 22      19      20      35      97 


diana_microt = gene2mir@data[gene2mir@data$database=="diana_microt",]
diana_microt$miRNA <- paste0(diana_microt$mature_mirna_id,"_",diana_microt$target_symbol)
diana_microt2<- rstatix::filter(diana_microt,!duplicated(diana_microt$miRNA))
table(diana_microt2$target_symbol)
# CYP1B1   DEGS1   HSDL2 OSBPL1A    PTEN 
# 40       4      11      14     176 

write.csv(mirdb2,"00.mirdb_mirna.csv",row.names = F) 
write.csv(diana_microt2,"00.diana_microt_mirna.csv",row.names = F) 


comon <- data.frame(symbol=base::intersect(mirdb2$miRNA,diana_microt2$miRNA))
mirdb2 <- mirdb2%>%dplyr::select(mature_mirna_id,target_symbol,miRNA)
diana_microt2 <- diana_microt2%>%dplyr::select(mature_mirna_id,target_symbol,miRNA)

names(diana_microt2) <- names(mirdb2)
mirna<- rbind(mirdb2,diana_microt2)
mirna <- mirna[which(mirna$miRNA%in%comon$symbol),]
mirna2<- rstatix::filter(mirna,!duplicated(mirna$miRNA))
mirna2<-mirna2[,c("mature_mirna_id","target_symbol","miRNA")]
colnames(mirna2)<-c("mirna","mrna","mirna_mrna")

write.csv(mirna2,"01.mirna.csv",row.names = F)

table(mirna2$mrna)
# CYP1B1  DEGS1  HSDL2   PTEN 
# 7      3      6     49 
length(unique(mirna2$mirna))
# 63

# #lncRNA-------------------------------------------------------------------
rm(list=ls())

library(dplyr)
dir.create("lncRNA")
mirna2 <- read.csv("01.mirna.csv")
for( RNA in mirna2$mirna){
  
  file=paste(getwd(),"/lncRNA/",RNA,".txt",sep="")
  
  link=paste("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA=",RNA,
             "&clipExpNum=1&degraExpNum=0&pancancerNum=1&programNum=1&program=None&target=all&cellType=all",sep="")
  download.file(link,file)
  Sys.sleep(1)
}

merge <- c()
for (RNA in mirna2$mirna) {
  aa <- read.table(paste(getwd(), "/lncRNA/",RNA,".txt",sep = "") ,header=T,
                   quote="",sep="\t",dec=".",
                   comment.char="#",na.strings =c("NA"),fill=T )
  merge <- rbind(merge,aa)
}

merge <- na.omit(merge)
table(merge$geneType)
merge <- merge%>%filter(geneType=="lncRNA")

write.csv(merge,"02.lncRNA.csv")

lncRNA <- merge%>%dplyr::select(miRNAname,geneName)
lncRNA$miRNA_lncRNA <- paste0(lncRNA$miRNAname,"_",lncRNA$geneName)
lncRNA2<- rstatix::filter(lncRNA,!duplicated(lncRNA$miRNA_lncRNA))
length(unique(lncRNA2$geneName))# 202
length(unique(lncRNA2$miRNAname))# 25
names(lncRNA2)[1:2] <- c("mirna","lncRNA")


####netplot-------

library(tidyverse)
library(magrittr)
library(tidygraph)
library(ggraph)
library(igraph)
library(MetBrewer)
library(ggforce)

names(lncRNA2)
lncRNA1 <- lncRNA2 %>%  group_by(mirna) %>%     
  slice_head(n = 5) %>%    # 选每个组的前  行
  ungroup()
length(unique(lncRNA1$lncRNA))

mi<-mirna2
lnc<-lncRNA1
mi<-mi[,2:1]
colnames(mi)<-c('val1','val2')
lnc<-lnc[,1:2]
colnames(lnc)<-c('val1','val2')

ce<-rbind(mi,lnc)
write.csv(ce,'03.cerna_plot.csv',row.names = F,quote = F)


lncname <- unique(lnc$val2) ## 71
miname <- unique(mi$val2) ## 63
keygene <-  read.csv('../06.expressionROC//keygenes.csv')[,1]

df1 <- ce %>%
  dplyr::select(1) %>%
  dplyr::rename("Var" = "val1") %>%
  bind_rows(
    ce %>%
      dplyr::select(2) %>%
      dplyr::rename("Var" = "val2")
  ) %>%
  separate(col = Var, into = c("gene", "Type"), sep = " ") %>%
  mutate(
    Type = case_when(
      gene %in% keygene ~ "mRNA",
      # grepl("hsa-miR-", gene) ~ "miRNA",
      gene %in% lncname ~ "lncRNA",
      TRUE ~ "miRNA"
    )
  )

nodes <- df1 %>% distinct()

edges <- ce %>%
  separate(col=val1, into = c("val1"), sep = " ") %>%
  separate(col=val2, into = c("val2"), sep = " ") %>%
  set_colnames(c("from","to"))

graph <- graph_from_data_frame(edges, nodes, directed = FALSE)
V(graph)$degree <- igraph::degree(graph, mode = "all")
tidy_graph <- tidygraph::as_tbl_graph(graph)


color_mapping <- c("miRNA" = "#0daddd",
                   "lncRNA" =  "#7F7FFF",
                   "mRNA" = "#f35e5a")

set.seed(12345)
p <- ggraph::ggraph(graph = tidy_graph,layout = "kk") +# layout = 'kk', maxiter = 50
  geom_edge_link(width = 0.5, color = "#bac4d0", alpha = 0.8) +
  geom_node_point(aes(size = sqrt(degree), color = Type, shape = Type),
                  alpha = 0.7, show.legend = TRUE) +
  # geom_node_text(aes(label = name,
  #                    filter = !name %in% keygene),
  #                repel = TRUE,
  #                size = 2.5,
  #                fontface = "bold",
  #                color = "gray10") +
  geom_node_text(aes(label = name ),#,filter = name %in% keygene
                 repel = TRUE,
                 size = 2.5,  #
                 fontface = "bold",
                 color = "gray10") +
  scale_color_manual(values = color_mapping) +
  scale_edge_width(range = c(0.5, 1.5)) +
  scale_size(range = c(5, 10)) +
  ggraph::theme_graph() +
  theme(legend.position = "right")

p
library(Cairo)
CairoPDF("04.lnc_mi_mRNA_net.pdf", width=13, height=12)
print(p)
dev.off()
png("04.lnc_mi_mRNA_net.png", width = 13, height = 12,res = 600,units = 'in')
print(p)
dev.off()

num_nodes <- length(V(graph))

num_edges <- length(E(graph))

cat("节点数量:", num_nodes, "\n") # 138
cat("边数量:", num_edges, "\n") # 188






#####TF-------
rm(list = ls())
tfraw <- read.csv("05.TF_mRNA_RegNetwork.csv")
tfraw$tf <- paste0(tfraw$ID,"_",tfraw$Target)
tfraw2<- rstatix::filter(tfraw,!duplicated(tfraw$tf))
table(tfraw2$Target)
# CYP1B1   DEGS1   HSDL2 OSBPL1A    PTEN 
# 7       6       6       5      20 

tf <- tfraw2%>%select(ID,Target)
colnames(tf)<-c('val1','val2')



library(tidyverse)
library(magrittr)
library(tidygraph)
library(ggraph)
library(igraph)
library(MetBrewer)
library(ggforce)

keygene <-  read.csv('../06.expressionROC/keygenes.csv')[,1]

df1 <- tf%>% 
  select(1) %>% 
  rename("Var" = "val1") %>% 
  mutate(across(Var, as.character)) %>%  #
  bind_rows(
    tf %>% 
      select(2) %>% 
      rename("Var" = "val2") %>% 
      mutate(across(Var, as.character))  #
  ) %>% 
  separate(col = Var, into = c("gene", "Type"), sep = " ") %>% 
  mutate(
    Type = case_when(
      gene %in% keygene ~ "mRNA",
      TRUE ~ "TF"
    )
  )

nodes <- df1 %>% distinct()

edges <- tf %>%
  separate(col=val1, into = c("val1"), sep = " ") %>%
  separate(col=val2, into = c("val2"), sep = " ") %>%
  set_colnames(c("from","to"))

graph <- graph_from_data_frame(edges, nodes, directed = FALSE)
V(graph)$degree <- igraph::degree(graph, mode = "all")
tidy_graph <- tidygraph::as_tbl_graph(graph)


color_mapping <- c("TF" = "#38bdf8",
                   "mRNA" = "#f35e5a")

set.seed(6666)
p <- ggraph::ggraph(graph = tidy_graph, layout = "fr") +
  geom_edge_link(width = 0.5, color = "#bac4d0", alpha = 0.5) +
  geom_node_point(aes(size = sqrt(degree), color = Type, shape = Type),
                  alpha = 0.7, show.legend = TRUE) +
  
  # 
  geom_node_text(aes(label = name,
                     filter = !name %in% keygene),
                 repel = TRUE,
                 size = 2.5,
                 fontface = "bold",
                 color = "gray10") +
  
  # 
  geom_node_text(aes(label = name,
                     filter = name %in% keygene),
                 repel = TRUE,
                 size = 3,  # 
                 fontface = "bold",
                 color = "gray10") +
  
  scale_color_manual(values = color_mapping) +
  scale_edge_width(range = c(0.5, 1.5)) +
  scale_size(range = c(5, 10)) +
  ggraph::theme_graph() +
  theme(legend.position = "right")
p

library(Cairo)
CairoPDF("06.tf_net.pdf", width=7, height=6)
print(p)
dev.off()
png("06.tf_net.png", width = 7, height = 6,res = 600,units = 'in')
print(p)
dev.off()

num_nodes <- length(V(graph))

num_edges <- length(E(graph))

cat("节点数量:", num_nodes, "\n") # 36
cat("边数量:", num_edges, "\n") # 44




