library(dplyr)
library(MASS)
library(network) # network
library(GGally) # ggnet2
library(ggplot2) # ggnet2
library(ggnet) # ggnet2
library(sna) # ggnet2 - degree
library(ppcor) # partial cor mat
library(glmnet) # lasso
library(gridExtra)
library(extrafont) #font

full.data = read.csv("C:\\Users\\User\\Desktop\\jb\\network\\full.data.csv")
data = full.data %>% dplyr::select(-year,-kedcd,-BZC_CD_num,-China,-USA,-OECD,-bsi,-Korea) %>% na.omit()

#precision matrix
prec.mat = solve(cov(data))
#partial correlation matrix
part.cor.mat = -prec.mat/sqrt(outer(diag(prec.mat),diag(prec.mat),"*"))

#data scale
data = scale(data,center = TRUE, scale = FALSE)


#lasso and linear regression
n = dim(data)[1] ; p=dim(data)[2]
lasso.P.mat = matrix(0,p,p)
lm.P.mat = matrix(0,p,p)
for(i in 1:p){
  print(i)
  s = sqrt(data[,i]%*%data[,i]/n) 
  alpha = 0.05
  phi.x = alpha/(2*p^2)
  phi = qnorm(1-phi.x)
  lam = drop(2*s)/sqrt(n)*phi
  
  glm.fit = glmnet(data[,-i],data[,i],lambda = lam)
  lasso.P.mat[i,-i]= -prec.mat[i,i]*coef(glm.fit)[-1]
  
  
  formula = as.formula(paste(colnames(data)[i],"~.",sep = ""))
  lm.fit = lm(formula,data = as.data.frame(data))
  lm.P.mat[i,-i] = -prec.mat[i,i]*coef(lm.fit)[-1]
}

#diagnal = 0
diag(lm.P.mat) = diag(prec.mat)
diag(lasso.P.mat) = 0


# ### a is neighborhood of b and b is neghborhood of a
# for(j in 1:p){
#   lasso.P.mat[j,][lasso.P.mat[j,]*lasso.P.mat[,j] ==0] = 0
#   lasso.P.mat[,j][lasso.P.mat[j,]*lasso.P.mat[,j] ==0] = 0
# }

### a is neighborhood of b or b is neghborhood of a
for(j in 1:p){
  lasso.P.mat[j,][lasso.P.mat[j,]+lasso.P.mat[,j] ==0] = 0
  lasso.P.mat[,j][lasso.P.mat[j,]+lasso.P.mat[,j] ==0] = 0
}

###################################
lasso.P.mat=t(lasso.P.mat)######### 매우 중요, 회귀계수와 그래프 화살표 방향 잘 살펴보기!!!!!
###################################
ind = colSums(lasso.P.mat)!=0|rowSums(lasso.P.mat)!=0
lasso.net.mat = lasso.P.mat[ind,ind]
colnames(lasso.net.mat) = colnames(data)[ind]
colnames(lm.P.mat) = colnames(data)



### Network plot of lm
lm.net = network(lm.P.mat,directed = TRUE)


g1=ggnet2(lm.net,size = 6,arrow.size = 6,arrow.gap = 0.025, edge.size = 1,edge.color = "darkgray",label = TRUE) + 
  ggtitle("Network plot of lm") + 
  theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))



#lasso.net
lasso.net = network(lasso.net.mat,directed = TRUE)

# degree : freeman, outdegree, indegree
lasso.net %v% "node.color" = ifelse(sna::degree(lasso.net,cmode = "freeman")>10, "tomato", "steelblue")
lasso.net %v% "label.size" = ifelse(sna::degree(lasso.net,cmode = "freeman")>10, 5, 3)
lasso.net %v% "node.size" = ifelse(sna::degree(lasso.net,cmode = "freeman")>10, 5, 3)

### Network plot of lasso1
g2=ggnet2(lasso.net,arrow.size = 6,arrow.gap = 0.025, edge.size = 1,edge.color = "darkgray",
       legend.size = 10, legend.position = "bottom",
       label = TRUE, label.size = "label.size",
       node.color = "node.color") + 
  ggtitle("Network plot of lasso1") + 
  theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))


### Network plot of lasso2
g3=ggnet2(lasso.net,node.size = 7,arrow.size = 6,arrow.gap = 0.025, edge.size = 1,edge.color = "darkgray",
       legend.size = 10, legend.position = "bottom",
       label = TRUE, label.size = "label.size",
       node.color = "grey15",
       label.color = "node.color")+
  theme(panel.background = element_rect(fill = "grey15")) + 
  ggtitle("Network plot of lasso2") +
  theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))


### Network plot of lasso3
g4=ggnet2(lasso.net,color="node.color",size=0,arrow.size = 6,arrow.gap = 0.04, edge.size = 1,edge.color = "gray35")+
  geom_point(aes(color = color),size=12,color="white")+
  geom_point(aes(color = color), size = 12, alpha = 0.5)+
  geom_point(aes(color = color), size = 8) +
  geom_text(aes(label = colnames(lasso.net.mat)), color = "black", fontface = "bold") + 
  ggtitle("Network plot of lasso3")+
  theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))



grid.arrange(g1,g2,g3,g4)


#######################################################################################
library(glasso)
library(CVglasso) #tune

#lamda tune
set.seed(2019)
cvglasso=CVglasso(data,K=5,crit.cv = "loglik",nlam = 100,lam.min.ratio = 0.0000001)
cvglasso$Tuning

glasso = glasso::glasso(cov(data), rho = cvglasso$Tuning[2])
colnames(glasso$wi) = colnames(data)
gl.ind=colSums(glasso$wi!=0)!=1
glasso.net = network(glasso$wi[gl.ind,gl.ind],directed = TRUE)


glasso.net %v% "node.color" = ifelse(sna::degree(glasso.net,cmode = "freeman")>10, "tomato", "steelblue")
ggnet2(glasso.net,color="node.color",size=0,arrow.size = 6,arrow.gap = 0.04, edge.size = 1,edge.color = "gray35")+
  geom_point(aes(color = color),size=12,color="white")+
  geom_point(aes(color = color), size = 12, alpha = 0.5)+
  geom_point(aes(color = color), size = 8) +
  geom_text(aes(label = colnames(glasso$wi[gl.ind,gl.ind])), color = "black", fontface = "bold") + 
  ggtitle("Network plot of lasso3")+
  theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 15, color = "darkblue"))

