#' @importFrom ape is.binary.tree ace multi2di root pic ace compute.brlen is.ultrametric bind.tree dist.nodes read.tree write.tree

## function does fast estimation of ML ancestral states using  
## written by Liam J. Revell 2012, 2013, 2015
fastAnc<-function(tree,x,vars=FALSE,CI=FALSE,...){
  if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
  if(length(class(tree)>1)) class(tree)<-"phylo"
  if(hasArg(anc.states)) anc.states<-list(...)$anc.states
  else anc.states<-NULL
  if(!is.null(anc.states)){
    nodes<-as.numeric(names(anc.states))
    tt<-tree
    for(i in 1:length(nodes)){
      M<-matchNodes(tt,tree,method="distances",quiet=TRUE)
      ii<-M[which(M[,2]==nodes[i]),1]
      tt<-bind.tip(tt,nodes[i],edge.length=0,where=ii)
    }
    x<-c(x,anc.states)
  } else tt<-tree
  if(!is.binary.tree(tt)) btree<-multi2di(tt)
  else btree<-tt
  M<-btree$Nnode
  N<-length(btree$tip.label)
  anc<-v<-vector()
  for(i in 1:M+N){
    a<-multi2di(root(btree,node=i))
    anc[i-N]<-ace(x,a,method="pic")$ace[1]
    names(anc)[i-N]<-i
    if(vars||CI){
      picx<-pic(x,a,rescaled.tree=TRUE)
      b<-picx$rescaled.tree
      d<-which(b$edge[,1]==(length(b$tip.label)+1))
      v[i-N]<-(1/b$edge.length[d[1]]+1/b$edge.length[d[2]])^(-1)*mean(picx$contr^2)
      names(v)[i-N]<-names(anc)[i-N]
    }
  }
  if(!is.binary.tree(tree)||!is.null(anc.states)){
    ancNames<-matchNodes(tree,btree,method="distances",quiet=TRUE)
    anc<-anc[as.character(ancNames[,2])]
    names(anc)<-ancNames[,1]
    if(vars||CI){ 
      v<-v[as.character(ancNames[,2])]
      names(v)<-ancNames[,1]
    }
  }
  obj<-list(ace=anc)
  if(vars) obj$var<-v
  if(CI){ 
    obj$CI95<-cbind(anc-1.96*sqrt(v),anc+1.96*sqrt(v))
    rownames(obj$CI95)<-names(anc)
  }
  if(length(obj)==1) obj<-obj$ace
  class(obj)<-"fastAnc"
  obj
}


# function to match nodes between trees
# written by Liam J. Revell 2012, 2013, 2015
matchNodes<-function(tr1,tr2,method=c("descendants","distances"),...){
  if(!inherits(tr1,"phylo")||!inherits(tr1,"phylo")) stop("tr1 & tr2 should both be objects of class \"phylo\".")
  if(hasArg(quiet)) quiet<-list(...)$quiet
  else quiet<-FALSE
  method<-method[1]
  method<-matchType(method,c("descendants","distances"))
 if(method=="distances"){
    if(hasArg(tol)) tol<-list(...)$tol
    else tol<-1e-6
    if(hasArg(corr)) corr<-list(...)$corr
    else corr<-FALSE
    if(corr) tr1$edge.length<-tr1$edge.length/max(nodeHeights(tr1))
    if(corr) tr2$edge.length<-tr2$edge.length/max(nodeHeights(tr2))
    D1<-dist.nodes(tr1)[1:length(tr1$tip),1:tr1$Nnode+length(tr1$tip)]
    D2<-dist.nodes(tr2)[1:length(tr2$tip),1:tr2$Nnode+length(tr2$tip)]
    rownames(D1)<-tr1$tip.label
    rownames(D2)<-tr2$tip.label
    common.tips<-intersect(tr1$tip.label,tr2$tip.label)
    D1<-D1[common.tips,]
    D2<-D2[common.tips,]
    Nodes<-matrix(NA,tr1$Nnode,2,dimnames=list(NULL,c("tr1","tr2")))
    for(i in 1:tr1$Nnode){
      if(corr) z<-apply(D2,2,function(X,y) cor(X,y),y=D1[,i])
      else z<-apply(D2,2,function(X,y) 1-sum(abs(X-y)),y=D1[,i])
      Nodes[i,1]<-as.numeric(colnames(D1)[i])
      if(any(z>=(1-tol))){
        a<-as.numeric(names(which(z>=(1-tol))))
        if(length(a)==1) Nodes[i,2]<-a
        else {
          Nodes[i,2]<-a[1]
          if(!quiet) warning("polytomy detected; some node matches may be arbitrary")
        }
      }
    }
  }
  return(Nodes)
}

# function adds a new tip to the tree
# written by Liam J. Revell 2012, 2013, 2014, 2015
bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL,position=0,...){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  use.edge.length<-if(is.null(tree$edge.length)) FALSE else TRUE
  if(use.edge.length==FALSE) tree<-compute.brlen(tree)
   if(is.null(where)) where<-length(tree$tip.label)+1
  if(where<=length(tree$tip.label)&&position==0){
    pp<-1e-12
    if(tree$edge.length[which(tree$edge[,2]==where)]<=1e-12){
      tree$edge.length[which(tree$edge[,2]==where)]<-2e-12
      ff<-TRUE
    } else ff<-FALSE
  } else pp<-position
  if(is.null(edge.length)&&is.ultrametric(tree)){
    H<-nodeHeights(tree)
    if(where==(length(tree$tip.label)+1)) edge.length<-max(H)
    else edge.length<-max(H)-H[tree$edge[,2]==where,2]+position
  }
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=tip.label,
            edge.length=edge.length,
            Nnode=1)
  class(tip)<-"phylo"
  obj<-bind.tree(tree,tip,where=where,position=pp)
  if(where<=length(tree$tip.label)&&position==0){
    nn<-obj$edge[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label)),1]
    obj$edge.length[which(obj$edge[,2]==nn)]<-obj$edge.length[which(obj$edge[,2]==nn)]+1e-12
    obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label))]<-0
    obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tree$tip.label[where]))]<-0
  }
  root.time<-if(!is.null(obj$root.time)) obj$root.time else NULL
  obj<-untangle(obj,"read.tree")
  if(!is.null(root.time)) obj$root.time<-root.time
  if(!use.edge.length) obj$edge.length<-NULL
  obj
}

## function finds the height of a given node
## written by Liam Revell 2014, 2015, 2016
nodeheight<-function(tree,node,...){
  if(hasArg(root.edge)) root.edge<-list(...)$root.edge
  else root.edge<-FALSE
  if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
  else ROOT<-0 
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(node==(length(tree$tip.label)+1)) h<-0
  else {
    a<-setdiff(c(getAncestors(tree,node),node),length(tree$tip.label)+1)
    h<-sum(tree$edge.length[sapply(a,function(x,e) which(e==x),e=tree$edge[,2])])
  }
  h+ROOT
}

# returns the heights of each node
# written by Liam J. Revell 2011, 2012, 2013, 2015, 2016
nodeHeights<-function(tree,...){
  if(hasArg(root.edge)) root.edge<-list(...)$root.edge
  else root.edge<-FALSE
  if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
  else ROOT<-0 
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  if(attr(tree,"order")!="cladewise"||is.null(attr(tree,"order"))) t<-reorder(tree)
  else t<-tree
  root<-length(t$tip.label)+1
  X<-matrix(NA,nrow(t$edge),2)
  for(i in 1:nrow(t$edge)){
    if(t$edge[i,1]==root){
      X[i,1]<-0.0
      X[i,2]<-t$edge.length[i]
    } else {
      X[i,1]<-X[match(t$edge[i,1],t$edge[,2]),2]
      X[i,2]<-X[i,1]+t$edge.length[i]
    }
  }
  if(attr(tree,"order")!="cladewise"||is.null(attr(tree,"order")))
    o<-apply(matrix(tree$edge[,2]),1,function(x,y) which(x==y),y=t$edge[,2])
  else o<-1:nrow(t$edge)
  return(X[o,]+ROOT)
}

# function 'untangles' (or attempts to untangle) a tree with crossing branches
# written by Liam J. Revell 2013, 2015
untangle<-function(tree,method=c("reorder","read.tree")){
  if(inherits(tree,"multiPhylo")){
    tree<-lapply(tree,untangle,method=method)
    class(tree)<-"multiPhylo"
  } else {
    if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
    obj<-attributes(tree)
    method<-method[1]
    if(method=="reorder") tree<-reorder(reorder(tree,"pruningwise"))
    else if(method=="read.tree"){
      tree<-read.tree(text=write.tree(tree))
    }
    ii<-!names(obj)%in%names(attributes(tree))
    attributes(tree)<-c(attributes(tree),obj[ii])
  }
  tree
}

# match type
# written by Liam J. Revell 2012
matchType<-function(type,types){
  for(i in 1:length(types))
    if(all(strsplit(type,split="")[[1]]==strsplit(types[i],split="")[[1]][1:length(strsplit(type,split="")[[1]])]))
      type=types[i]
    return(type)
}

## function gets ancestor node numbers, to be used internally by 
## written by Liam J. Revell 2014
getAncestors<-function(tree,node,type=c("all","parent")){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  type<-type[1]
  if(type=="all"){
    aa<-vector()
    rt<-length(tree$tip.label)+1
    currnode<-node
    while(currnode!=rt){
      currnode<-getAncestors(tree,currnode,"parent")
      aa<-c(aa,currnode)
    }
    return(aa)
  } else if(type=="parent"){
    aa<-tree$edge[which(tree$edge[,2]==node),1]
    return(aa)
  } else stop("do not recognize type")
}