library(readr)
data <- read_csv("~/GitHub/SDS385-course-work/Excercise 1/logistic regression/wdbc.csv",col_names = FALSE)
X=as.matrix(data[3:12])
y=as.vector(matrix(nrow=nrow(data),ncol=1))
for(i in 1:nrow(data)){
  if(data[i,2]=="M")y[i]=1
  else y[i]=0
}