# GPEDI
Software for converting genome data of population to pedigree data format and estimating coefficient of inbreeding.
if(!require(sequoia)) install.packages("sequoia")
GPEDI <- function(
                   data, 
				   target0 = NULL
				   ){
					 C = ncol(gen)
					 if(C <= 4){
					            if(is.null(target0)){
								                   print("-------------------result = GPEDI.pedOUT(data)-----------------")
												   
					                               result = GPEDI.pedOUT(data)
									               }else{
												         print("-------------------result =GPEDI.inbreeding(data,target0)-----------------")
												         result =GPEDI.inbreeding(data,target0) 
												        }
							   }else{
		      if(all(c(0,1,2) %in% data)){
									print("-------------------result = GPEDI.genealogy(data)-----------------")
									result = GPEDI.genealogy(data)
									}else{
									     print("-------------------kong-----------------")
										 }
									}
				   return(result)
				   }
GPEDI.genealogy <- function(
                   data 
				   ){

                    file = data
					 pairsss=rownames(file)
					 pairs = data.frame(ID1 = character(),
									   ID2 = character(),
									   stringsAsFactors = FALSE
									   )
					 for (i in 1:(length(pairsss)-1)){
					  
														  for (j in (i+1):length(pairsss))  {
																								 pair<- data.frame(ID1 = pairsss[i],
																								 ID2 = pairsss[j]) 
																								 pairs<- rbind(pairs, pair)
																								}
														}
					 PairL = CalcPairLL(Pairs=pairs,GenoM = file, Err = 1e-04, Plot=FALSE)
                     ped = NULL
                     PairL = PairL[c(which(PairL$TopRel!="U")),]
					 PO=PairL[-c(which(PairL$TopRel!="PO")),]
					 if(nrow(PO)!=0){
					 PO_1 = PO
					 familytree = 0
                     H_all = NULL
					 while(nrow(PO_1) != 0){
					 H = PO_1[1,]
					 H1 = H[[1]]
					 H2 = H[[2]]
					 H_c = H
					 
					 h = NULL
                     familytree = familytree+1
					 H_all[[familytree]] = H_c
]
						    while(nrow(H)!= 0){
							                h = NULL
											for(i in 1:nrow(H)){
											H1 = H[i,1]
											H2 = H[i,2]
											
											PO_1 = PO_1[!(PO_1$ID1==H1 & PO_1$ID2==H2),]
											PO_1 = PO_1[!(PO_1$ID2==H1 & PO_1$ID1==H2),]
											H_3=PO_1[c(which(PO_1$ID1==H1)),]
											H_4=PO_1[c(which(PO_1$ID2==H1)),]
											H_5=PO_1[c(which(PO_1$ID1==H2)),]
											H_6=PO_1[c(which(PO_1$ID2==H2)),]
											h[[i]]= rbind(H_3,H_4,H_5,H_6)
											h[[i]] = unique(h[[i]])
											
																}
											 H= do.call(rbind, h)
											 H=unique(H)
											 
											 H_all[[familytree]]= rbind(H,H_all[[familytree]])
											 H_all[[familytree]]=unique(H_all[[familytree]])
											 
											}
										}
                     
                     for(v in 1:length(H_all)){
					 PO = H_all[[v]]
					 PO_M = c(PO$ID1,PO$ID2)
					 PO_M = unique(PO_M)
					 PairL_1 = PairL[PairL$ID1 %in% PO_M, ]
					 PairL_1 = PairL_1[PairL_1$ID2 %in% PO_M, ]
					 FS_1 = PairL_1$TopRel
					 PO_T = PO
					 FS=PairL_1[-c(which(PairL_1$TopRel!="FS")),]
					 while(nrow(FS)!=0){
                     FS = FS[,1:2]
				for(i in 1:length(FS$ID1)){
					                             sson_1=FS[i,]$ID1
												 sson_3=c(FS[c(which(FS$ID1==sson_1)),])
                                                 sson_4=c(FS[c(which(FS$ID2==sson_1)),])
                                                 sson_2 = c(unlist(sson_3),unlist(sson_4))
												 sson_2= unique(sson_2)
												 for(j in 1:length(sson_2)){
												 PO_T$ID1[PO_T$ID1 == sson_2[j]] = sson_1
												 PO_T$ID2[PO_T$ID2 == sson_2[j]] = sson_1
												  }
												  FS = FS[!(FS$ID1 %in% sson_2), ]
												  FS = FS[!(FS$ID2 %in% sson_2), ]
												  }
												  }
												  
					 PO_T=PO_T[!duplicated(PO_T[, c("ID1", "ID2")]), ]
					 ID_1 = c(PO_T$ID1)
					 ID_2 = c(PO_T$ID2)
					 ID_0=c(ID_1,ID_2)
					 name=NULL
					 min_ID=table(ID_0)
					 min_I=1
					 F1=names(min_ID[min_ID == min_I])
					 ID = F1
					 k=1
					 name [[k]]=ID

					 while(length(ID)!=0){
								  k=k+1
								  ID_1= PO_T[PO_T$ID1 %in% ID, ]
								  ID_1= c(ID_1$ID2)
								  
								  ID_2= PO_T[PO_T$ID2 %in% ID, ]
								  ID_2= c(ID_2$ID1)
								  
								  ID_0=c(ID_1,ID_2)
								  ID_0=unique(ID_0)
								  
								  ID_0=ID_0[!(ID_0 %in% ID)]
								  ID_all=unlist(name)
								  ID = ID_0[!(ID_0 %in% ID_all)]
								  
								  if(length(ID)!=0){
												  name[[k]]=ID
												  }else {
									                    }
								  
								  } 
					 HS=PairL_1[-c(which(PairL_1$TopRel!="HS")),]
					 name_p=NULL
					 ped_m=NULL
					 j = length(name)
					 B = name[[j]]
					 p=1
					 name_p[[p]]=B
					 while (!all(is.na(B))){						  
											 B = na.omit(B)
											 B = B[B != "NA"]
											 name_p2=NULL
											 for (j in 1:length(B)){
											                          C = B[j]
																      HS_1= HS[which(HS$ID1==C), ]
																	  
																	  HS_1=unique(HS_1)
																	  
																	  HS_1=unique(HS_1)
																	  
																	  HS_2= HS[which(HS$ID2==C), ]
																	  
																	  HS_2=unique(HS_2)
																	  
																	  HS_2= c(HS_2$ID1)
																	  
																      FS_1= FS[which(FS$ID1==C), ]

																	  FS_1=unique(FS_1)
																	  
																	  FS_1= c(FS_1$ID2)
																	  
																	  FS_2= FS[which(FS$ID2==C), ]
																	  
																	  FS_2=unique(FS_2)
																	  
																	  FS_2= c(FS_2$ID1)
																	  
																      B = c(B,HS_1,HS_2,FS_1,FS_2)
																	  B = unlist(B)
																      B=unique(B)
																  
																  }
											 for (k in 1:length(B)){
																	  C = B[k]
																	  ID_1= PO[which(PO$ID1==C), ]
																	  ID_1=unique(ID_1)
																	  ID_1= c(ID_1$ID2)
																	  ID_2= PO[which(PO$ID2==C), ]
																	  ID_2=unique(ID_2)
																	  ID_2= c(ID_2$ID1)
																												 
																	 ID_all2=unlist(name_p)
																	 ID_1=ID_1[!(ID_1 %in% ID_all2)]
																	 ID_2=ID_2[!(ID_2 %in% ID_all2)]												 						                         
																	 m=c(C,ID_1,ID_2)
																	 if(length(m)!=3){
																	                 if(length(m)<3){																					 
																									 m[(length(m)+1):3]=NA
																					                 m=matrix(m, nrow = 1,ncol=3 )
																									 ped_m=rbind(ped_m,m)
																									 while(length(ID_1)==0 & length(ID_2)==0){
																																	 ID_1=NA
																																	 ID_2=NA
																																	 }
																									 name_p2[[k]]=c(ID_1,ID_2)
																									}else{
																					                      m=c(C,NA,NA)
																										  m=matrix(m, nrow = 1,ncol=3 )
																										  ped_m=rbind(ped_m,m)
																										  name_p2[[k]]=NA
																										  }																
																					  
																					  }else{
																					     
																							 m=matrix(m, nrow = 1,ncol=3 )
																							 ped_m=rbind(ped_m,m)
																							 while(length(ID_1)==0 & length(ID_2)==0){
																																	 ID_1=NA
																																	 ID_2=NA
																																	 }
																							 name_p2[[k]]=c(ID_1,ID_2)					
																						 }				   
																	}
											 p=p+1
											 name_p[[p]]=c(unlist(name_p2))
											 B=name_p[[p]]		 
											}					 
					 ped_m=unique(ped_m)
					 sex_ID = matrix(0,length(ped_m[,1]),1)
					 S = c(PO$Sex1,PO$Sex2)
		
		if (!is.null(ped_m)){
					 for(i in 1:nrow(ped_m)){
								                    ID_sex = ped_m[i,1]
													ID_sex_row1 = PO[which(PO$ID1==ID_sex), ]
													ID_sex_row2 = PO[which(PO$ID2==ID_sex), ]
													ID_sex1 = unique(c((ID_sex_row1$Sex1),(ID_sex_row2$Sex2)))		  
													if(length(ID_sex1)==1){
													sex_ID[i,1] = ID_sex1
													}else{
													sex_ID[i,1] = NA
													}
												}
					
					 sex_ID[sex_ID == 3] = NA
					 sex_ID[sex_ID == 1] = 0
					 sex_ID[sex_ID == 2] = 1
					 ped_m = cbind(ped_m,sex_ID)
					 if(3 %in% S){
					              colnames(ped_m)=c("ID","parents","parents","gender")
								  }else{
										colnames(ped_m)=c("ID","sire","dam","gender")
					 
										for(i in 1:length(ped_m[,2])){
													m = ped_m[i,2]
													sex_m1 = PO[which(PO$ID1==m), ] 
													sex_m2 = PO[which(PO$ID2==m), ]
													ID_sex1 = unique(c(unlist(ID_sex_row1$Sex1),unlist(ID_sex_row2$Sex2)))
													if(ID_sex1==2){
													              ped_m[i,2] = ped_m[i,4]
																  ped_m[i,3] = m
																  }else{
																	   }

													}

					                     }
                      }else{
					       ped_m = rbind(ped_m,c(NA,NA,NA,NA))
					       colnames(ped_m)=c("ID","parents","parents","gender")
					       }            										 
					 ped_m=unique(ped_m)
					 familynames=matrix(v, nrow =nrow(ped_m),ncol=1 )
					 colnames(familynames) = "familytrees"
					 ped_m = cbind(ped_m,familynames)
					 ped[[v]] = ped_m
					}
                    ped_V = do.call(rbind, ped)
}else{
     ped_m = NULL
     ped_m = rbind(ped_m,c(NA,NA,NA,NA))
	 colnames(ped_m)=c("ID","parents","parents","gender")
	 ped_V = ped_m
	 
	 }
return(ped_V)
					}
GPEDI.pedOUT <- function(
                          data
						  ){
			file = data
			target=file[,1]
			target = na.omit(target)
			coi=matrix(NA,length(target),3)
			for(m in 1:length(target))
				{
				 target0=target[m]
				 j=0
				 E=target0
				 while(!all(is.na(E)))
					{
					 search.result=file[file[,1]%in%E,]
					 j=j+1 
					 E=c(search.result[,3],search.result[,4])#
						E = unique(E)
						}
					 coi[m,1]=target0 
					 
					 coi[m,3]=j 
				}
				g = max(coi[,3])
				file_1 = file[,-2]
				file_1 = na.omit(file_1)
				son = file_1$ID
				son=unique(son)
				parents = c(file_1$sire,file_1$dam)
				parents=unique(parents)
				born_son = length(son)/length(parents)
				parents_1 = c(file_1$dam,file_1$sire)
				parents_1 = na.omit(parents_1)
				parents_1_1 = table(parents_1)
				A = max(parents_1_1)
				parents_1_1 = as.matrix(parents_1_1)
                C = sum(parents_1_1[,1])			
				D = length(parents_1_1)
				B = C/D
				B = round(B, 1)
				gender = file$gender
				sire = sum(gender=="0")
				dam = sum(gender=="1")
				unknow = sum(gender=="NA")
				information=matrix(NA,6,1)
				rownames(information)=c("maxbirths","Avebirths","sire","dam","unknowgender","generations")
				information[1,1] = A
				information[2,1] = B
				information[3,1] = sire
				information[4,1] = dam
				information[5,1] = unknow
				information[6,1] = g
				return(information)
				}	
GPEDI.inbreeding <-function(data,
	target0
	)
	{

	file = data
	target=file[,1]
	target = na.omit(target)
	g=0
	E=target0
	for(i in 1:length(target))
	{
	s=E
	if(s%in%file[,1]&!is.na(target0))
		{
		 search.result=file[file[,1]%in%s,]
		 g=g+1 
		 E=search.result[1,3]#
		}else{
			 g=g
			 }
	}
	gener = g
	tar.s=file[file[,1]%in%target0,3]
	tar.d=file[file[,1]%in%target0,4]
	cross.generation=FALSE
	
	inbreed=matrix(NA,gener,2^(gener-1))
	inbreed[1,1]=target0
	
	for(i in 1:(gener-1))
	{
		search.store=inbreed[i,!is.na(inbreed[i,])]
		search.num=length(search.store)
		for(j in 1:search.num)
		{
			if(search.store[j]%in%file[,1])
			{	
			   search.result=file[file[,1]%in%search.store[j],]
			}else{
			   search.result=matrix(NA,1,4)
			}
			inbreed[i+1,2*j-1]=search.result[1,3]
			inbreed[i+1,2*j]=search.result[1,4]
		
		}
	}

	line=list()
	mm=0
	if(cross.generation)
	{
	 print("This function is pending ...")
	}else{
	   
	   for(i in 2:gener)
	   {
		  search.store=as.character(inbreed[i,!is.na(inbreed[i,])])
		  dup.id=unique(search.store[duplicated(search.store)])
		  if(length(dup.id)>0)
		  {
			 for(dup in 1:length(dup.id))
			 {
				dup.index=grep(dup.id[dup],search.store)
				dup.num=length(dup.index)
				for(m in 1:dup.num)
				{
					mm=mm+1
					line[[mm]]=inbreed[i,dup.index[m]]
				}
				for(j in (i-1):2)
				{
				   dup.index.2=ceiling(dup.index/2)
				   for(m in 1:dup.num)
				   {
					line[[mm-m+1]]=append(line[[mm-m+1]],inbreed[j,dup.index.2[m]])
				   }
				   dup.index=dup.index.2
				}
			 }
			
			 
		  }else{
			 next
		  } 

	   }


	}


	n=length(line)
	    if(n!=0){
				 first.id=NULL
				 for(i in 1:n)
				 {
				  first.id=append(first.id,line[[i]][1])
				 }
				 last.id=NULL
				 for(i in 1:n)
				 {
				  last.id=append(last.id,line[[i]][length(line[[i]])])
				 }
				 table.matrix=as.matrix(table(last.id))
				if(length(table.matrix)!=1)
					 {
					 as.matrix(table.matrix[order(table.matrix[,1]),])    
					 min.id=rownames(table.matrix)[1]
					 max.id=rownames(table.matrix)[2]
					 min.num=as.numeric(table.matrix[min.id,])
					 max.num=as.numeric(table.matrix[max.id,])
					 whole.line=list()
					 mm=0
					 for(i in 1:min.num)
					 {
					  min.index=grep(min.id,last.id)[i]
					  min.first=line[[min.index]][1]
					  
						last.id.2=last.id
						last.id.2[!first.id%in%min.first]=NA
						last.id.2[min.index]=NA
						max.index=grep(max.id,last.id.2)
						for(k in max.index)
						{
						  max.line=line[[k]]
						  max.line2=NULL
						  for(f in length(max.line):1)
						  {
							max.line2=append(max.line2,max.line[f])
						  }
						  mm=mm+1
						  whole.line[[mm]]=append(max.line2,line[[min.index]][-1])
						}
					  
					 }
					if(!is.null(whole.line)){
						 x=0
						 i=1
						 for(i in (1:length(whole.line))){
										  T1=whole.line[[i]]
										  x1=(1/2)^length(T1)
										  i=1+1
										  x=x+x1
										  }
										  return(x)
						 }else{
							   x=0
							  }
					
					
					}else{
						 x=0
						 }
				}else{
						 x=0
						 }
return(x)
}

