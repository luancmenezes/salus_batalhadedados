setwd("/Users/scp93/Desktop/data/")
#---------------------------------------------------------------------------------------
# leitura dos dados
#---------------------------------------------------------------------------------------
library(readr)
data <- read_delim("/Users/scp93/Desktop/Pasta Sem Título/saude/HSL/AutIntHos-DataSUS-2013-2018.csv", 
                   ";", escape_double = FALSE, trim_ws = TRUE)
lvl<-levels(as.factor(data$CID_PRINCIPAL))

#---------------------------------------------------------------------------------------
# selecao dos cids referentes ao cancer de mama
#---------------------------------------------------------------------------------------
cids<-apply(as.matrix(0:9), 1, function(n){
  txt<-paste("C50", n, sep = "")
  
  return(txt)
})

idx<-apply(cbind(data, 1:dim(data)[1]), 1, function(x){
  cid<-x[26]
  if(sum(cid==cids)){
    out<-x[50]
  } else {
    out<-NA
  }
  
  return(out)
})

idx<-as.numeric(idx[!is.na(idx)])
data<-data[idx,]

#---------------------------------------------------------------------------------------
# incluindo a variavel idade
#---------------------------------------------------------------------------------------
idade<-2018-floor(data$PACIENTE_DATA_NASCIMENTO/10000)
data[,which(names(data)=="PACIENTE_DATA_NASCIMENTO")]<-idade
names(data)[which(names(data)=="PACIENTE_DATA_NASCIMENTO")]<-"IDADE"

#---------------------------------------------------------------------------------------
# definicao do valor maximo a investir que impacte na chance do paciente sobreviver
#---------------------------------------------------------------------------------------
n<-102
quant<-seq(0,1, length.out = n)
valores<-quantile(data$VALOR_TOTAL, quant)

prop.morte<-NULL
# proporcao de vivos por intervalo de investimetno
for(i in 1:(length(valores)-2)){
  idx<-which(data$VALOR_TOTAL>=valores[i] & data$VALOR_TOTAL<=valores[i+1])
  prop.morte[i]<-mean(data$MORTE[idx])
}

valores<-apply(as.matrix(1:length(prop.morte)), 1, function(i){
  mean(valores[i:(i+1)])
})

x<-valores[-(n-2)]
y<-prop.morte[-(n-2)]

# regressao para definir a partir de qual ponto o investimento nao interfere mais na chance de morte do paciente
require(segmented)
neg.x<--x
fit<-lm(y~1)
estim_a0<-segmented(fit, seg.Z= ~neg.x, psi=NA, control = seg.control(K =1, stop.if.error = T))
plot(x, y, pch=16, xlab="valor gasto", ylab="chance de morte")
lines(x, predict(estim_a0), col=2)

# ponto maximo de investimento
-estim_a0$psi[,2]
# intervalo de 95% para o investimento maximo
sort(-confint(estim_a0)$neg.x)
abline(v=sort(-confint(estim_a0)$neg.x)[c(1,3)], col="blue")

# ---------------------------------------------------------------------------------
# simulacao das caracterizacoes de um tumor, a partir de um BD real
# ---------------------------------------------------------------------------------
idx.maligno<-which(data$MORTE==1 | data$CID_PRINCIPAL=="C508" | data$CID_PRINCIPAL=="C509")
n.morte<-length(idx.maligno)
data_tumo<-read.csv("~/Desktop/data/data.csv")[,-c(1,33)]
idx_mal<-which(data_tumo$diagnosis=="M")
range<-rbind(apply(data_tumo[idx_mal,-1], 2, min),
             apply(data_tumo[idx_mal,-1], 2, max))

data_input_mal<-t(apply(as.matrix(1:n.morte), 1, function(vazio){
  out<-apply(range, 2, function(range){
    runif(1, range[1], range[2])
  })
  return(out)
}))

idx_beg<-which(data_tumo$diagnosis=="B")
range_b<-rbind(apply(data_tumo[idx_beg,-1], 2, min),
               apply(data_tumo[idx_beg,-1], 2, max))

prev<-prop.table(table(data_tumo$diagnosis))[2]
data_input_ben<-t(apply(as.matrix(1:(dim(data)[1]-n.morte)), 1, function(vazio){
  if(rbinom(1, 1, prev)){
    out<-apply(range, 2, function(range){
      runif(1, range[1], range[2])
    })
  } else {
    out<-apply(range_b, 2, function(range){
      runif(1, range[1], range[2])
    })
  }
  return(out)
}))

# treinamento do modelo para classificar os tumores a partir de um BD real --> acuracia 0.9415
setwd("/Users/scp93/Desktop/data/")
library(data.table)

dados_caracterizacaoTumor<-read.csv("~/Desktop/data/data.csv")[,-c(1,33)]
dados_caracterizacaoTumor$diagnosis<-as.character(dados_caracterizacaoTumor$diagnosis)
idx<-which(dados_caracterizacaoTumor$diagnosis=="B")
dados_caracterizacaoTumor$diagnosis[idx]<-0
dados_caracterizacaoTumor$diagnosis[-idx]<-1

dados_caracterizacaoTumor$diagnosis<-as.numeric(dados_caracterizacaoTumor$diagnosis)

set.seed(1)
n<-dim(dados_caracterizacaoTumor)[1]
train<-sample(1:n,  .7*n)

# QDA
require(MASS)
fit2<-qda(diagnosis~., data = dados_caracterizacaoTumor[train,])
pred<-predict(fit2, dados_caracterizacaoTumor[-train,])

mat.conf<-table(pred$class, dados_caracterizacaoTumor$diagnosis[-train])
sens<-mat.conf[2,2]/max(sum(mat.conf[,2]), 1)
espec<-mat.conf[1,1]/max(sum(mat.conf[,1]), 1)
acuracia<-sum(diag(mat.conf))/sum(mat.conf)

c(sens, espec, acuracia)

# organizacao dos dados
data_input<-data.frame(rbind(data_input_mal, data_input_ben))

pred<-predict(fit2, data_input)
pred_final<-pred$class

# idx.maligno

data_input_aux<-data_input
data_input_aux[idx.maligno,]<-data_input[1:n.morte,]
data_input_aux[-idx.maligno,]<-data_input[-(1:n.morte),]

# tumo_class<-data_input_aux
pred_final_aux<-rep(NA, dim(data)[1])
pred_final_aux[idx.maligno]<-pred_final[1:n.morte]
pred_final_aux[-idx.maligno]<-pred_final[-(1:n.morte)]

pred_final_aux<-pred_final_aux-1

data<-cbind(data, pred_final, data_input_aux)

#---------------------------------------------------------------------------------------
# definicao de custo
#---------------------------------------------------------------------------------------
require(glmnet)
grid=10^seq(10,-2,length=100)
fit3<-glmnet(as.matrix(data_input_aux), data$VALOR_TOTAL, alpha=.5, lambda=grid)
set.seed (1)
cv.out<-cv.glmnet(as.matrix(data_input_aux), data$VALOR_TOTAL,alpha=.5)
bestlam<-cv.out$lambda.min
pred_val<-predict(fit3,s=bestlam ,newx=as.matrix(data_input_aux))
