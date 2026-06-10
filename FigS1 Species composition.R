setwd("C:\\Users\\lenovo\\Desktop")
#library(leaps)
#sink("output.txt")#将分析结果输出
########################################
library(rtry)
library(patchwork)
library(MuMIn)
library(piecewiseSEM)#structure equation model
#library(semPlot)#绘制结构方程模型图
library(car)# for vif test
library(glmmTMB)#generalise linear mixed model
library(MASS)
library(ggplot2)
library(moments)# 计算偏度系数# 计算峰度系数 
library(ppcor)#偏相关

########################################
#####figS1 物种组成的图形
df=read.csv("gr3.csv");Lf3="Perennial";Fre3=length(na.omit(df$AGBmax))
df=read.csv("gr5.csv");Lf5="Annual";Fre5=length(na.omit(df$AGBmax))
df=read.csv("gr6.csv");Lf6="Biennial";Fre6=length(na.omit(df$AGBmax))

PER3=Fre3/(Fre3+Fre5+ Fre6);PER5=Fre5/(Fre3+Fre5+ Fre6);PER6=Fre6/(Fre3+Fre5+ Fre6)
df4 <- data.frame(Category = c(Lf3,Lf5,Lf6),Value = c(PER3, PER5, PER6))
df4$Category=factor(df4$Category,levels = c(Lf3,Lf5,Lf6))
PER3=sprintf("%.2f%%", PER3 * 100)#保留两位小数并显示为百分比形式
PER5=sprintf("%.2f%%", PER5 * 100);PER6=sprintf("%.2f%%", PER6 * 100)

composition=ggplot(df4, aes(x = Category, y = Value)) + geom_bar(stat = "identity",fill = "light blue") +
  labs(title = "All sites", x = "", y = "Relative frequency") +
  annotate("text",x=1, y=0.55, label=PER3,angle = 0,size=10,family="serif")+
  annotate("text",x=2, y=df4$Value[2]+0.1, label=PER5,angle = 0,size=10,family="serif")+
  annotate("text",x=3, y=df4$Value[3]+0.1, label=PER6,angle = 0,size=10,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0,size=20),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm")) +
  scale_y_continuous(breaks=seq(0,1,by=0.25),limits=c(0,1))
#ggsave("Species composition.pdf",composition,width = 20, height = 20, units = "cm", dpi = 300) 


#按地点分开做物种组成图
df=read.csv("gr2.csv");compos=list()
sitename=unique(df$Location...Site.Name);lifesty=c("perennial" ,"annual" , "biennial" )

for (i in 1:12){
  df2=rtry_select_row(df,
                      (Location...Site.Name %in% sitename[i]),
                      getAncillary = TRUE)
  
  perennial_df=rtry_select_row(df2,
                               (lifehistory %in% lifesty[1]),
                               getAncillary = TRUE)
  annual_df=rtry_select_row(df2,
                            (lifehistory %in% lifesty[2]),
                            getAncillary = TRUE)
  biennial_df=rtry_select_row(df2,
                              (lifehistory %in% lifesty[3]),
                              getAncillary = TRUE)
  
  
  Fre_perennial=length(na.omit(perennial_df$Onsetflower))
  Fre_annual=length(na.omit(annual_df$Onsetflower))
  Fre_biennial=length(na.omit(biennial_df$Onsetflower))
  
  tot=Fre_perennial+Fre_annual+Fre_biennial
  lifesty2= c("Perennial" ,"Annual" , "Biennial" )
  
  df5 <- data.frame(Category = lifesty2,Value = c(Fre_perennial,Fre_annual,Fre_biennial)/tot)
  df5$Category=factor(df5$Category,levels = lifesty2)
  labba=sprintf("%.2f%%", df5$Value * 100)#保留两位小数并显示为百分比形式
  
  compos[[i]]=ggplot(df5, aes(x = Category, y = Value)) + geom_bar(stat = "identity",fill="light blue") +
    labs(title = sitename[i], x = "", y = "Relative frequency") +
    annotate("text",x=1, y=0.55, label=labba[1],angle = 0,size=10,family="serif")+
    annotate("text",x=2, y=df5$Value[2]+0.1, label=labba[2],angle = 0,size=10,family="serif")+
    annotate("text",x=3, y=df5$Value[3]+0.1, label=labba[3],angle = 0,size=10,family="serif")+
    theme(axis.text=element_text(size=30,family="serif",colour = "black"),
          axis.text.x = element_text(vjust=0.0,size=20),
          axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
          axis.title.x = element_text(vjust=0.0),
          axis.line=element_line(size = 0.5, colour = "black", linetype=1),
          axis.ticks.length = unit(0.25, "cm"),
          axis.ticks = element_line(size = 1, color="black"),
          plot.background = element_rect(fill = "white", color = NA), 
          panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                    margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
          plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm")) +
    scale_y_continuous(breaks=seq(0,1,by=0.25),limits=c(0,1))
    #scale_y_continuous(breaks=seq(0,1200,by=300),limits=c(0,1200))
  
  
}
ggsave("Species composition all and each site.pdf",composition/(compos[[1]]|compos[[2]]|compos[[3]])/
         (compos[[4]]|compos[[5]]|compos[[6]])/(compos[[7]]|compos[[8]]|compos[[9]])/(compos[[10]]|compos[[11]]|compos[[12]]),
       width = 40, height = 60, units = "cm", dpi = 300) 


#################################
