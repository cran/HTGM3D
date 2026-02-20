#' HTGM3Ddriver
#'
#' @import minimalistGODB
#' @import GoMiner
#' @import HTGM
#' @import HTGM2D
#' @import randomGODB
#' @import grDevices
#' @import stats
#' @import R2HTML
#' @import rgl
#' @import vprint
#' @import stringr
#'
#' @description driver to invoke HTGM3D()
#'
#' @param dir character string full path name to the directory acting as result repository
#' @param geneList character vector of user-supplied genes of interest
#' @param GOGOA3 return value of subsetGOGOA()
#' @param thresh1  numerical acceptance threshold for individual ontologies
#' @param thresh3 numerical acceptance threshold for joint ontology
#' @param mn integer min category size threshold passed to trimGOGOA3()
#' @param mx integer max category size threshold passed to trimGOGOA3()
#' @param pcgMN integer param passed to pruneCatGenes
#' @param pcgMX integer param passed to pruneCatGenes
#' @param verbose integer vector representing vprint classes
#'
#' @details
#' suggested standardized class codes for vprint()
#' -1 = developer debugging only
#' 0 = constitutively turned on
#' 1 = help for new user
#' 2 = follow progress of long computation
#' 3 = primary results
#' 4 = meta information (e.g. dims of a matrix before and after trimming)
#' 5 = warnings
#' 6 = errors
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#'
#' #load("data/cluster52.RData")
#' geneList<-cluster52
#' dir<-tempdir()
#' l3<-HTGM3Ddriver(dir,geneList,GOGOA3,thresh1=3,thresh3=3,
#'  mn=2,mx=20000,pcgMN=2,pcgMX=200,verbose=1:5)
#' }
#' 
#' @return returns matrix containing information that provides the input needed for running plot3d()
#' 
#' @export
HTGM3Ddriver<-
  function(dir,geneList,GOGOA3,thresh1,thresh3,mn,mx,pcgMN,pcgMX,verbose) {
    
    stamp<-gsub(":","_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
    title<-"HTGM3D"
    subd<-sprintf("%s/%s_%s",dir,title,stamp)
    dir.create(subd)
    
    meta<-list()
    meta$subd<-subd
    meta$geneList<-geneList
    #meta$species<-"human"
    meta$species<-GOGOA3$species
    meta$thresh1<-thresh1
    meta$thresh3<-thresh3
    meta$min<-mn
    meta$max<-mx
    meta$pcgMN<-pcgMN
    meta$pcgMX<-pcgMX
    
    save(meta,file=sprintf("%s/meta.RData",subd))
    
    ontology<-"biological_process"
    
    vprint(-1,verbose,c("BEFORE gene list DB",mn,mx,length(geneList),dim(GOGOA3$ontologies[["biological_process"]])))
    
    vprint(-1,verbose,"geneListDistHitters before preprocessDB")
    geneListDistHitters(geneList,GOGOA3,ontology,FALSE)
    pp<-preprocessDB(geneList,GOGOA3,ontology,mn,mx,thresh=.5,verbose)
    vprint(-1,verbose,c("names pp",names(pp)))
    geneList<-pp$sampleList
    GOGOA3<-pp$GOGOA3
    
    vprint(-1,verbose,c("AFTER gene list DB",mn,mx,length(geneList),dim(GOGOA3$ontologies[["biological_process"]])))
    
    vprint(-1,verbose,"geneListDistHitters after preprocessDB")
    geneListDistHitters(geneList,GOGOA3,ontology,FALSE)
    
    l3<-list()
    
    mat<-HTGM3D(subd,geneList,GOGOA3,thresh1,thresh3,mn,mx,pcgMN,pcgMX,verbose)
    
    # list l contains the table() results of frequency of appearance of each category in joint histogram
    l<-list()
    
    l[[names(GOGOA3$ontologies)[1]]]<-sort(table(mat[,"cat1"]),decreasing=TRUE)
    l[[names(GOGOA3$ontologies)[2]]]<-sort(table(mat[,"cat2"]),decreasing=TRUE)
    l[[names(GOGOA3$ontologies)[3]]]<-sort(table(mat[,"cat3"]),decreasing=TRUE)
    
    vprint(4,verbose,"Frequency of Appearance of Each Category in Joint Histogram:")
    for(i in 1:3) {
      vprint(4,verbose,names(GOGOA3$ontologies)[i])
      vprint(4,verbose,sort(l[[i]],decreasing=TRUE))
    }
    
    # list l2 contains the axis coordinate number to be used for the plotting order for each category
    #x_l<-l
    #save(x_l,file="data/x_l.RData")
    l2<-catNum3(l)
    
    vprint(4,verbose,"Axis Coordinate Number for the Plotting Order of each Category per Ontology:")
    for(i in 1:3) {
      vprint(4,verbose,names(GOGOA3$ontologies)[i])
      vprint(4,verbose,l2[[i]])
    }
    
    #x_mat<-mat
    #save(x_mat,file="data/x_mat.RData")
    # re-number categories to show up in 3D plot in meaningful order
    
    mat3d<-plot3Dmat(mat,l2)
    
    # more informative colnames for mat3d
    colnames(mat3d)[2]<-names(GOGOA3$ontologies[1])
    colnames(mat3d)[4]<-names(GOGOA3$ontologies[2])
    colnames(mat3d)[6]<-names(GOGOA3$ontologies[3])
    
    f<-sprintf("%s/%s",subd,"HTGM3D.html")
    if(file.exists(f))
      unlink(f)
    mat3d2<-cbind(1:nrow(mat3d),mat3d)
    colnames(mat3d2)<-c('#',colnames(mat3d))
    rownames(mat3d2)<-NULL
    HTML(mat3d2,file=f)
    
    sys<-sprintf("open -a Safari.app %s",f)
    system(sys)
    
    l3$mat3d<-mat3d
    l3$jointFreq<-l
    l3$plot3Dindex<-l2
    
    #x_mat3d<-mat3d
    #save(x_mat3d,file="data/x_mat3d.RData")
    save(l3,file=sprintf("%s/l3.RData",subd))
    
    vprint(1,verbose,"Summary spreadsheet should already be open for you in Safari browser")
    vprint(1,verbose,"You can use super cool interative 3D visualization by copy and paste into R Console:")
    vprint(1,verbose,"interactWithGraph3D(l3$mat3d)")
    
    return(l3)
  }

#' catNum3
#' 
#' @description assign the axis coordinate number to be used for the plotting position for each category
#' 
#' @param l list each component corresponds to an ontology branch, and contains the decreasing sorted
#'  tabulation of output of the number of times that a category appears in a triplet
#'
#' @examples
#' #load("data/x_l.RData")
#' catNum3(x_l)
#' 
#' @details
#' a component of l is like:
#'  GO_0005515__protein_binding GO_0042802__identical_protein_binding   GO_0005178__integrin_binding  
#'  38                          4                                        3                                     
#' 
#' #' a component of l1 is like:
#'  GO_0005515__protein_binding GO_0042802__identical_protein_binding  GO_0005178__integrin_binding
#'  1                           2                                      3 
#' 
#' @return returns a list each component corresponds to an ontology branch, and contains a vector of category plotting positions
#' 
#' @export
catNum3<-
  function(l) {
    l1<-list()
    
    for(i in names(l)) {
      v<-1:length(l[[i]])
      names(v)<-names(l[[i]])
      l1[[i]]<-v
    }
    
    return(l1)
  }

#' HTGM3D
#'
#' @description compute matrix to use as input to plot3d()
#'
#' @param dir character string full path name to the directory acting as result repository
#' @param geneList character vector of user-supplied genes of interest
#' @param GOGOA3 return value of subsetGOGOA()
#' @param thresh1  numerical acceptance threshold for individual ontologies
#' @param thresh3 numerical acceptance threshold for joint ontology
#' @param mn integer min category size threshold passed to trimGOGOA3() 
#' @param mx integer max category size threshold passed to trimGOGOA3() 
#' @param pcgMN integer param passed to pruneCatGenes
#' @param pcgMX integer param passed to pruneCatGenes
#' @param verbose integer vector representing vprint classes
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' 
#' geneList<-cluster52
#' dir<-tempdir()
#' mat3d<-HTGM3D(dir,geneList,GOGOA3,thresh1=3,thresh3=3,mn=2,
#'  mx=10000,pcgMN=2,pcgMX=200,verbose=1:5)
#' }
#' 
#' @return returns matrix containing information that provides the input needed for running plot3d()
#' 
#' @export
HTGM3D<-
  function(dir,geneList,GOGOA3,thresh1,thresh3,mn,mx,pcgMN,pcgMX,verbose) {  	
    m<-list()
    
    vprint(4,verbose,"Number of Terms in Ontologies Before and After Pruning:")
    for(ontology in names(GOGOA3$ontologies)) {
      vprint(4,verbose,ontology)
      m[[ontology]]<-catGenes(geneList,GOGOA3,ontology)
      vprint(4,verbose,"HTGM3D after catGenes - cats matching gene list")
      vprint(4,verbose,dim(m[[ontology]]))
      
      #x_cg<-catGenes(geneList,GOGOA3,ontology)
      #save(x_cg,file="data/x_cg.RData")
      m[[ontology]]<-pruneCatGenes(catGenes(geneList,GOGOA3,ontology),pcgMN,pcgMX)
      vprint(4,verbose,"HTGM3D after catGenes AND pruneCatGenes - cats matching gene list")
      vprint(4,verbose,c(mn,mx))
      vprint(4,verbose,dim(m[[ontology]]))
      
      #x_m11<-m
      #save(x_m11,file="data/x_m11.RData")
    }
    
    mat<-Jaccard3(dir,m,thresh1,thresh3,verbose)
    if(nrow(mat)==0 | is.null(dim(mat))) {
      vprint(5,verbose,"HTGM3D RETURNING NULL")
      return(NULL)
    }
    
    mat<-insertCatSize(mat,GOGOA3)
    
    #x_jmat<-mat
    #save(x_jmat,file="data/x_jmat.RData")
    
    if(nrow(mat)==0 | is.null(dim(mat))) {
      vprint(5,verbose,"HTGM3D RETURNING NULL")
      return(NULL)
    }
    
    return(mat[order(mat[,"lint"],decreasing=TRUE),])
  }

#' pruneCatGenes
#' 
#' @description eliminate those categories to which no genes map
#' 
#' @param m the return value of catGenes()
#' @param mn integer min category size threshold passed to trimGOGOA3()
#' @param mx integer max category size threshold passed to trimGOGOA3()
#' 
#' @examples
#' #load("data/x_cg.RData")
#' m<-pruneCatGenes(x_cg,2,200)
#' 
#' @return returns pruned version of matrix m
#' 
#' @export
pruneCatGenes<-
  function(m,mn=2,mx=200) {
    zero_rows <- (rowSums(m) < mn) | (rowSums(m) > mx)
    
    return(m[!zero_rows,])
  }

#'Jaccard3
#'
#' @description compute the number of genes in the intersection of categories from 3 ontology branches
#' 
#' @param dir character string full path name to the directory acting as result repository
#' @param m return value of pruneCatGenes()
#' @param thresh1 parameter passed to Jaccard3()
#' @param thresh3 parameter passed to Jaccard3()
#' @param verbose integer vector representing vprint classes
#'
#' @examples
#' #load("data/x_m11.RData")
#' Jaccard3(tempdir(),x_m11,thresh1=3,thresh3=3,verbose=1:5)
#' 
#' @return returns matrix tabulating genes in the intersection of categories from 3 ontology branches
#'  also has side effect of saving files containing those genes 
#' 
#' @export
Jaccard3<-
  function(dir,m,thresh1=3,thresh3=3,verbose=2) {
    subdir<-sprintf("%s/%s",dir,"hyperGenes")
    
    if(!dir.exists(subdir))
      dir.create(subdir)
    
    cn<-c("cat1","cat2","cat3","lint","n1","n2","n3","hyperGenes")
    mat<-matrix(ncol=length(cn))
    rownames(mat)<-NA
    colnames(mat)<-cn
    count<-0
    nr<-nrow(m[[1]])
    for(cat1 in rownames(m[[1]])) {
      #count<-count+1
      #print(c("IN JACCARD3D",count,nr))
      if(sum(m[[1]][cat1,])<=thresh1)
        next
      n<-0
      for(cat2 in rownames(m[[2]])) {
        if(sum(m[[2]][cat2,])<=thresh1)
          next
        for(cat3 in rownames(m[[3]])) {
          if(sum(m[[3]][cat3,])<=thresh1)
            next
          
          w<-which(m[[1]][cat1,]&m[[2]][cat2,]&m[[3]][cat3,]==TRUE)					
          if(length(w)<=thresh3)
            next
          
          count<-count+1
          vprint(2,verbose,c("IN JACCARD3D",count,nr))
          
          label<-sprintf("%s__%s__%s.txt",substr(cat1,1,10),substr(cat2,1,10),substr(cat3,1,10))
          # names(w) are the desired common gene names
          writeLines(sort(names(w)),sprintf("%s/%s",subdir,label))				
          
          v<-c(cat1,cat2,cat3,length(w),sum(m[[1]][cat1,]),sum(m[[2]][cat2,]),sum(m[[3]][cat3,]),"HyperGenes")
          rn<-rownames(mat)
          mat<-rbind(mat,v)
          rownames(mat)<-c(rn,label)
        } # for(cat3 in rownames(m[[3]]))
      } # for(cat2 in rownames(m[[2]]))
    } # for(cat1 in rownames(m[[1]]))
    mat[,"hyperGenes"]<-sprintf('<a href = "%s/%s/%s">%s</a>',dir,"HyperGenes",rownames(mat),rownames(mat))
    
    return(mat[-1,]) # first row of mat is just a bunch of NA's
  }

#' plot3Dmat
#' 
#' @description compute x,y,z coordinates for each triplet
#' 
#' @param mat return value of HTGM3D()
#' @param l return value of catNum3()
#'
#' @examples
#' #load("data/x_mat.RData")
#' #load("data/x_l.RData")
#' p3<-plot3Dmat(x_mat,x_l)
#' 
#' @return augmented version of matrix containing x,y,z coordinates for each triplet
#' 
#' @export
plot3Dmat<-
  function(mat,l) {
    
    # prepare matrix for submitting to plot3D
    mat3d<-matrix(nrow=nrow(mat),ncol=5)
    colnames(mat3d)<-c("x","y","z","xyz","col")
    for(r in 1:nrow(mat)) {
      x<-l[[1]][mat[r,"cat1"]]
      y<-l[[2]][mat[r,"cat2"]]
      z<-l[[3]][mat[r,"cat3"]]
      xyz<-sprintf("%s_%s_%s",x,y,z)
      mat3d[r,]<-c(x,y,z,xyz,mat[r,"lint"])
    }
    
    #return(cbind(mat[,1:3],mat3d,mat[,4:(ncol(mat))]))
    return(cbind(mat[,1:6],mat3d,mat[,7:(ncol(mat))]))
  }

#' blackBodyRadiationColors
#' 
#' @description set up color scale for black body spectrum
#' 
#' @param x numeric should be between 0 (black) and 1 (white)
#' @param max_value numeric maximum value to be used for scaling
#'
#' @examples
#' colors.blackBody <- rev(blackBodyRadiationColors(seq(0.3,1,length.out=20)))
#' 
#' @details
#' I obtained this by copy and paste from internet (reference unknown)
#' 
#' @return returns no value, but has side effect of generating color map
#' 
#' @export
blackBodyRadiationColors<-
  function(x, max_value=1) {
    # x should be between 0 (black) and 1 (white)
    # if large x come out too bright, constrain the bright end of the palette
    #     by setting max_value lower than 1
    foo <- colorRamp(c(rgb(0,0,0),rgb(1,0,0),rgb(1,1,0),rgb(1,1,1)))(x*max_value)/255
    apply(foo,1,function(bar)rgb(bar[1],bar[2],bar[3]))
  }

#' interactWithGraph3D
#' 
#' @description rotate the 3d graph and/or select a point within the 3D graph, and annotate that point
#' 
#' @param mat3d component of HTGM3Ddriver() output list
#' @param maxfract numeric upper threshold for category size to display
#' @param newWindow Boolean if TRUE open new window to avoid over writing current window 
#' @param verbose integer vector representing vprint classes
#' 
#' @examples
#' if(interactive()){
#' #load("data/x_mat3d.RData")
#' interactWithGraph3D(x_mat3d)
#' }
#' 
#' @return returns -1 if OS is OSX 26.2 which has glitch with rgl package
#' 
#' @export 
interactWithGraph3D<-
  function(mat3d,maxfract=1.00,newWindow=TRUE,verbose=TRUE) {
    
    OS<-whichOS()
    
    range<-list()
    range$x<-range(as.numeric(mat3d[,"x"]))
    range$y<-range(as.numeric(mat3d[,"y"]))
    range$z<-range(as.numeric(mat3d[,"z"]))
    
    r<-max(range$x[2]-range$x[1],range$y[2]-range$y[1],range$z[2]-range$z[1])
    
    # set newWindow == TRUE to avoid over writing current window 
    if(newWindow)
      open3d()
    colors.blackBody <- rev(blackBodyRadiationColors(seq(0.3,1,length.out=20)))
    w<-which((as.numeric(mat3d[,2])<maxfract) & (as.numeric(mat3d[,4])<maxfract) & (as.numeric(mat3d[,6])<maxfract))
    mat3d<-mat3d[w,]
    mx<-as.integer(max(mat3d[,"col"]))
    
    lcb<-length(colors.blackBody)
    
    oldPar<-par3d("windowRect")
    # cannot use on.exit(), must be performed explicitly upon user entering 'x'
    par3d("windowRect"=c(50,50,1000,1000)) # adjust position and size of plot3d window
    
    if(OS=="OK"){
    ids<-plot3d(mat3d[,c("x","y","z")],xlab="molecular_function",ylab="cellular_component",
                zlab="biological_process",col=colors.blackBody[as.integer(lcb*as.integer(mat3d[,"col"])/mx)],size=10)
    
    ng<-0
    w<-NULL
    while(TRUE) {
      if(length(w)>0)
        option<-readline(prompt="Enter 'i' for interactive or a digit for line number or 'g' for gene list or 'x' to exit : ")
      else
        option<-readline(prompt="Enter 'i' for interactive or a digit for line number or 'x' to exit : ")
      if(option=="x") {
        par3d("windowRect"=oldPar)
        stop("user selected 'x' to exit")
      }
      
      if(str_detect(option, "^[:digit:]+$")) {
        w<-as.integer(option)
        sp<-as.integer(c(mat3d[w,"x"],mat3d[w,"y"],mat3d[w,"z"]))
        vprint(4,verbose,c("sp",sp))
        graphIt(mat3d,sp,w,r,verbose=TRUE)
      }
      
      if(option=="i") {
        sp<-selectpoints3d(ids["data"])
        vprint(4,verbose,c("sp",sp))
        sp3<-sprintf("%s_%s_%s",sp[1],sp[2],sp[3])
        w<-which(mat3d[,"xyz"]==sp3)
        graphIt(mat3d,sp,w,r,verbose=TRUE)
      }
      
      if(option=="g") {
        showGenes(mat3d,w,range,npad=25*ng)
        ng<-ng+1
      }
      
      next
    }
    }
    else if(OS=="BAD") { # just bare bones os 26.2, bypass interactive options
      print("glitch with rgl package in mac os Tahoe 26.2, just bare bones, bypass extra interactive options")
      print("3D will come up in web browser rather than in quartz window")
      print("I will restore full functionality when developers fix the gitch")
      plot3d(mat3d[,c("x","y","z")],xlab="molecular_function",ylab="cellular_component",
             zlab="biological_process",col=colors.blackBody[as.integer(lcb*as.integer(mat3d[,"col"])/mx)],size=10)
    }
    else { #total failure for os 26.3
      print("glitch with rgl package in mac os Tahoe 26.3, aborting until system developers fix")
    }
  }

#' graphIt
#' 
#' @description annotate a selected point in the 3d graphic
#' 
#' @param mat3d component of HTGM3Ddriver() output list
#' @param sp integer vector containing c(x,y,z) coordinated of point 
#' @param w integer line number within mat3d
#' @param r numeric max value of x,y,z ranges
#' @param verbose integer vector representing vprint classes
#' 
#' @return returns no value, but has side effect of annotating the 3D graph
#' 
#' @export 
graphIt<-
  function(mat3d,sp,w,r,verbose) {
    arrow3d(sp-r/10, sp, type = "rotation",  col = "green")
    
    text3d(sp[1],sp[2],sp[3],texts=mat3d[w,5],col='red',cex=1,pos=3)
    text3d(sp[1],sp[2],sp[3],texts=mat3d[w,3],col='red',cex=1,pos=3,offset=2)
    text3d(sp[1],sp[2],sp[3],texts=mat3d[w,1],col='red',cex=1,pos=3,offset=3.5)
    
    vprint(4,verbose,c(mat3d[w,c("lint","n1","n2","n3")]))
  }

#' showGenes
#' 
#' @description open gene list in a textEdit window
#' 
#' @param mat3d component of HTGM3Ddriver() output list
#' @param w integer line number within mat3d
#' @param range list of ranges
#' @param npad integer number of blanks for padding
#' 
#' @return returns no value, but has side effect of opening textEdit
#' 
#' @export 
showGenes<-
  function(mat3d,w,range,npad) {
    ss<-strsplit(mat3d[w,"hyperGenes"],'"')
    #sys<-sprintf("open %s",ss[[1]][2])
    #system(sys)
    
    x<-readLines(ss[[1]][2])
    lx<-length(x)
    ll<-range$z[1]
    for(i in 1:lx) {
      text3d(range$x[2],ll,range$z[1],sprintf("%s%s",str_pad(" ",npad),x[i]),pos=4,cex=.75)
      ll<-ll+1
    }
  }

#' insertCatSize
#' 
#' @description compute fraction of total genes in the entire ontology that map to each category,
#' and insert column into matrix mat 
#' 
#' @param mat return value of Jaccard3()
#' @param GOGOA3  return value of subsetGOGOA()
#' 
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' 
#' mat2<-insertCatSize(x_mat,GOGOA3)
#' }
#' 
#' @return returns augmented version of matrix mat
#' 
#' @export
insertCatSize<-
  function(mat,GOGOA3) {
    mat2<-mat
    cn<-names(GOGOA3$ontologies)
    for(i in 1:3) {
      tcats<-GOGOA3$stats$tcats[[cn[i]]]
      ngene<-vector("numeric",nrow(mat))
      names(ngene)<-mat[,i]
      for(j in 1:nrow(mat)) {
        ngene[j]<-tcats[names(ngene)[j],2]
      }
      mat2<-cbind(mat2[,1:(2*i-1)],round(ngene,digits=3),mat2[,(2*i):ncol(mat2)])
    }
    colnames(mat2)[1]<-"cat1"
    return(mat2)
  }

#' subMatDiffs
#' 
#' @description retrieve submatrices of common rownames and colnames, and report differences
#' 
#' @param m1 matrix
#' @param m2 matrix
#' @param verbose integer vector representing vprint classes
#'
#' @examples
#' rn<-c("a","b","c")
#' cn<-c("a","b","c","d")
#' m1<-matrix(1:12,nrow=length(rn),ncol=length(cn))
#' rownames(m1)<-rn
#' colnames(m1)<-cn
#' m2<-m1
#' subMatDiffs(m1,m2)
#' 
#' m3<-m1
#' m3[1,1]<-0
#' m3[1,2]<-0
#' subMatDiffs(m1,m3)
#' 
#' m4<-m3
#' colnames(m4)<-c("aa","b","c","d")
#' subMatDiffs(m1,m4)
#' 
#' @details
#' compare submatrices of m1 and m2 that have common rownames and colnames
#' 
#' @return returns no values
#' 
#' @export
subMatDiffs<-
  function(m1,m2,verbose=3) {
    rn<-intersect(rownames(m1),rownames(m2))
    cn<-intersect(colnames(m1),colnames(m2))
    vprint(3,verbose,"dims of original and submatrix:")
    vprint(3,verbose,dim(m1))
    vprint(3,verbose,dim(m2))
    vprint(3,verbose,c(length(rn),length(cn)))
    
    vprint(3,verbose,"setdiffs of original matrices rownames and colnames:")
    vprint(3,verbose,setdiff(rownames(m1),rownames(m2)))
    vprint(3,verbose,setdiff(rownames(m2),rownames(m1)))
    vprint(3,verbose,setdiff(colnames(m1),colnames(m2)))
    vprint(3,verbose,setdiff(colnames(m2),colnames(m1)))
    
    different_elements <- m1[rn,cn] != m2[rn,cn]
    w <- which(different_elements,arr.ind=TRUE)
    lw<-nrow(w)
    if(lw>0) {
      vprint(3,verbose,"differences in submatrices:")
      vprint(3,verbose,w)
      
      for(i in 1:nrow(w)) {
        r<-rn[w[i,"row"]]
        c<-cn[w[i,"col"]]
        vprint(3,verbose,c(i,r,c,m1[r,c],m2[r,c]))
      }
    }
    
    tot<-length(rn)*length(cn)
    vprint(3,verbose,c("summary: tot same diff: ",tot,tot-lw,lw))
  }

#' whichOS
#' 
#' @import salesforcer
#' 
#' @description determine if host machine is running OSX Tahoe 26.2
#' 
#' @details
#' as of Feb 6 2026 Tahoe 26.2 does not interact properly with rgl Xwindows package
#' 
#' @returns "BAD" Tahoe 26.2, "OK" otherwise
#' 
#' @export
whichOS<-
  function() {
    if(salesforcer::get_os()=="osx") {
      output <- system("sw_vers", intern = TRUE)
      ss<-strsplit(output[2],"\t\t")
      if(ss[[1]][2]=="26.2")
        return("BAD")
      if(ss[[1]][2]=="26.3")
        return("REALBAD")
    }
    return("OK")
  }


