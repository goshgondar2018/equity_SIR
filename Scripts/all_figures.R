library(lattice)
library(ggplot2)

set.seed(96)
output_transmission1=ESIR_transmission_matrix_scenarios('half within quintile')

df_all=data.frame(output_transmission1[1])
betas=as.matrix(data.frame(output_transmission1[2]))

# plot of a sample transmission matrix
betas_plot=levelplot(betas,col.regions = heat.colors(100)[length(heat.colors(100)):1],main="Transmission rate matrix",
                    xlab="Quintile i", ylab="Quintile j")
update(betas_plot, aspect=1)       

# plot the standardized S, I, and R proportions and their means across iterations for each quintile under the main transmission matrix scenario
df_all_grouped_quint1_and5 <- df_all%>%
  group_by(days) %>%
  summarise(S1_stand_mean = mean(S1_stand), S2_stand_mean=mean(S2_stand), S3_stand_mean=mean(S3_stand), 
                S4_stand_mean=mean(S4_stand),S5_stand_mean=mean(S5_stand),
                I1_stand_mean = mean(I1_stand), I2_stand_mean=mean(I2_stand), I3_stand_mean=mean(I3_stand), 
                I4_stand_mean=mean(I4_stand),I5_stand_mean=mean(I5_stand),
                R1_stand_mean = mean(R1_stand), R2_stand_mean=mean(R2_stand), R3_stand_mean=mean(R3_stand), 
                R4_stand_mean=mean(R4_stand),R5_stand_mean=mean(R5_stand),
                D1_stand_mean=mean(D1_stand), D2_stand_mean=mean(D2_stand), D3_stand_mean=mean(D3_stand), 
                D4_stand_mean=mean(D4_stand), D5_stand_mean=mean(D5_stand),
            S1_conf_lower=quantile(S1_stand,c(0.025,0.975))[[1]],
            S1_conf_upper=quantile(S1_stand,c(0.025,0.975))[[2]],
            I1_conf_lower=quantile(I1_stand,c(0.025,0.975))[[1]],
            I1_conf_upper=quantile(I1_stand,c(0.025,0.975))[[2]],
            R1_conf_lower=quantile(R1_stand,c(0.025,0.975))[[1]],
            R1_conf_upper=quantile(R1_stand,c(0.025,0.975))[[2]],
            D1_conf_lower=quantile(D1_stand,c(0.025,0.975))[[1]],
            D1_conf_upper=quantile(D1_stand,c(0.025,0.975))[[2]],
            S5_conf_lower=quantile(S5_stand,c(0.025,0.975))[[1]],
            S5_conf_upper=quantile(S5_stand,c(0.025,0.975))[[2]],
            I5_conf_lower=quantile(I5_stand,c(0.025,0.975))[[1]],
            I5_conf_upper=quantile(I5_stand,c(0.025,0.975))[[2]],
            R5_conf_lower=quantile(R5_stand,c(0.025,0.975))[[1]],
            R5_conf_upper=quantile(R5_stand,c(0.025,0.975))[[2]],
            D5_conf_lower=quantile(D5_stand,c(0.025,0.975))[[1]],
            D5_conf_upper=quantile(D5_stand,c(0.025,0.975))[[2]])

colors=c("S1"="blue", "I1"="red","R1"="green","D1"='brown')
colors2=c("S5"="blue","I5"="red","R5"="green", "D5"='brown')


ggplot(df_all_grouped_quint1_and5, aes(x=days))+
  geom_line(aes(y=S1_stand_mean, col="S1"),linetype='solid')+
  geom_ribbon(aes(ymin=S1_conf_lower, ymax=S1_conf_upper),fill="blue", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=I1_stand_mean, col="I1"),linetype='solid')+
  geom_ribbon(aes(ymin=I1_conf_lower, ymax=I1_conf_upper),fill="red", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=R1_stand_mean, col="R1"),linetype='solid')+
  geom_ribbon(aes(ymin=R1_conf_lower, ymax=R1_conf_upper),fill="green", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=D1_stand_mean, col="D1"),linetype='solid')+
  geom_ribbon(aes(ymin=D1_conf_lower, ymax=D1_conf_upper),fill="brown", alpha=0.2,linetype='dashed')+
  scale_y_continuous(name="Proportion of quintile population")+
  scale_color_manual(values = colors, name="")+
  theme_classic()+xlim(c(0,80))+#c(0,75)
  theme(axis.text.x =element_text(size=15),axis.title.x=element_blank(),
        axis.text.y =element_text(size=15),axis.title.y=element_blank(),
        legend.position='none')

ggsave("./figures/ESIR_main_quintile1_feb12-21.pdf")

ggplot(df_all_grouped_quint1_and5,aes(x=days))+
  geom_line(aes(y=S5_stand_mean,  col="S5"),linetype='solid')+
  geom_ribbon(aes(ymin=S5_conf_lower, ymax=S5_conf_upper),fill="blue", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=I5_stand_mean,  col="I5"),linetype='solid')+
  geom_ribbon(aes(ymin=I5_conf_lower, ymax=I5_conf_upper),fill="red", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=R5_stand_mean, col="R5"),linetype='solid')+
  geom_ribbon(aes(ymin=R5_conf_lower, ymax=R5_conf_upper),fill="green", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=D5_stand_mean, col="D5"),linetype='solid')+
  geom_ribbon(aes(ymin=D5_conf_lower, ymax=D5_conf_upper),fill="brown", alpha=0.2,linetype='dashed')+
  scale_y_continuous(name="Proportion of quintile population")+
  scale_color_manual(values = colors2, name="")+
  theme_classic()+xlim(c(0,80))+ #c(0,75)
  theme(axis.text.x =element_text(size=15),axis.title.x=element_blank(),
        axis.text.y =element_text(size=15),axis.title.y=element_blank(),
        legend.position='none')

ggsave("./figures/ESIR_main_quintile5_feb12-21.pdf")

## identify max I and final D for quintiles 1 and 5
max_I_quintile1=max(df_all_grouped_quint1_and5$I1_stand_mean)
max_I_quintile5=max(df_all_grouped_quint1_and5$I5_stand_mean)
final_D_quintile1=df_all_grouped_quint1_and5$D1_stand_mean[400]
final_D_quintile5=df_all_grouped_quint1_and5$D5_stand_mean[400]

# Boxplot of # infected by quintiles
all_infs_time_period<-df_all%>%
  mutate(time_period=ifelse(days%in%seq(1,20),'first',ifelse(days%in%seq(21,40),'second','third')))%>%
  filter(time_period=='first'|time_period=='second')%>%
  group_by(iteration,time_period)%>%
  mutate(max(I1_stand),max(I2_stand),max(I3_stand),max(I4_stand),max(I5_stand))%>%
  select(time_period,`max(I1_stand)`,`max(I2_stand)`,`max(I3_stand)`,`max(I4_stand)`,`max(I5_stand)`)

colnames(all_infs_time_period)[3:7]=c('I1_stand_max','I2_stand_max','I3_stand_max','I4_stand_max','I5_stand_max')

all_infs_time_period_final=all_infs_time_period%>%gather(quintile,
                                                         max_prop_infected, I1_stand_max:I5_stand_max)%>%
  mutate(quintile=ifelse(quintile=='I1_stand_max','1',ifelse(quintile=='I2_stand_max','2',
                         ifelse(quintile=='I3_stand_max','3',ifelse(quintile=='I4_stand_max','4','5')))))
ggplot(data=all_infs_time_period_final,aes(x=quintile, y=max_prop_infected,col=time_period))+
  geom_boxplot()+xlab("Quintile")+
  ylab("Max proportion infected")+
  theme(plot.title = element_text(hjust = 0.5,size=10))+
  theme(panel.background = element_blank(),axis.line=element_line(color="black"))

ggsave("./figures/boxplots_infected_two_periods_feb12-21.pdf")

## extract medians
medians_all_infs<-all_infs_time_period_final%>%
  group_by(time_period,quintile)%>%
  summarise(median_max_prop_infected=median(max_prop_infected))

## vaccination
ALL_scenarios_mean_LQ_UQ_final=ALL_scenarios_mean_LQ_UQ%>%filter(scenario!='homogeneous-base case')%>%
  filter(scenario!='only top quintile')
plot_all_vacc=ggplot(ALL_scenarios_mean_LQ_UQ_final)+
  geom_bar(aes(x=quintile,y=mean,fill=quintile),stat="identity")+
  geom_errorbar(aes(x=quintile,ymin=LQ,ymax=UQ))+
  theme_classic()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.text.y = element_text(size=6),axis.title.y=element_blank())

plot_all_vacc+facet_wrap(~scenario)

ggsave("./figures/all_vacc_final_version_feb12-21.pdf")

# Figures for the supplementary appendix

## plot all quintiles model dynamics for main transmission matrix scenario

df_all_quint1 <- df_all%>%
  group_by(days) %>%
  summarise(S_stand_mean = mean(S1_stand), 
            I_stand_mean = mean(I1_stand), 
            R_stand_mean = mean(R1_stand),
            D_stand_mean=mean(D1_stand),
            S_conf_lower=quantile(S1_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S1_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I1_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I1_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R1_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R1_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D1_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D1_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(1))

df_all_quint2 <- df_all%>%
  group_by(days) %>%
  summarise(S_stand_mean = mean(S2_stand), 
            I_stand_mean = mean(I2_stand), 
            R_stand_mean = mean(R2_stand),
            D_stand_mean=mean(D2_stand),
            S_conf_lower=quantile(S2_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S2_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I2_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I2_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R2_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R2_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D2_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D2_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(2))

df_all_quint3 <- df_all%>%
  group_by(days) %>%
  summarise(S_stand_mean = mean(S3_stand), 
            I_stand_mean = mean(I3_stand), 
            R_stand_mean = mean(R3_stand),
            D_stand_mean=mean(D3_stand),
            S_conf_lower=quantile(S3_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S3_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I3_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I3_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R3_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R3_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D3_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D3_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(3))

df_all_quint4 <- df_all%>%
  group_by(days) %>%
  summarise(S_stand_mean = mean(S4_stand), 
            I_stand_mean = mean(I4_stand), 
            R_stand_mean = mean(R4_stand),
            D_stand_mean=mean(D4_stand),
            S_conf_lower=quantile(S4_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S4_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I4_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I4_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R4_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R4_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D4_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D4_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(4))

df_all_quint5 <- df_all%>%
  group_by(days) %>%
  summarise(S_stand_mean = mean(S5_stand), 
            I_stand_mean = mean(I5_stand), 
            R_stand_mean = mean(R5_stand),
            D_stand_mean=mean(D5_stand),
            S_conf_lower=quantile(S5_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S5_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I5_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I5_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R5_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R5_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D5_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D5_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(5))

df_all_quint=rbind.data.frame(df_all_quint1,df_all_quint2,df_all_quint3,df_all_quint4,df_all_quint5)
colors3=c("S"="blue","I"="red","R"="green", "D"='brown')

plot_all_quints=ggplot(df_all_quint, aes(x=days))+
  geom_line(aes(y=S_stand_mean, col="S"),linetype='solid')+
  geom_ribbon(aes(ymin=S_conf_lower, ymax=S_conf_upper),fill="blue", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=I_stand_mean, col="I"),linetype='solid')+
  geom_ribbon(aes(ymin=I_conf_lower, ymax=I_conf_upper),fill="red", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=R_stand_mean, col="R"),linetype='solid')+
  geom_ribbon(aes(ymin=R_conf_lower, ymax=R_conf_upper),fill="green", alpha=0.2,linetype='dashed')+
  geom_line(aes(y=D_stand_mean, col="D"),linetype='solid')+
  geom_ribbon(aes(ymin=D_conf_lower, ymax=D_conf_upper),fill="brown", alpha=0.2,linetype='dashed')+
  scale_color_manual(values = colors3, name="")+
  theme_classic()+xlim(c(0,80))+
  theme(axis.text.x =element_text(size=15),axis.title.x=element_blank(),
        axis.text.y =element_text(size=15),axis.title.y=element_blank(),
        legend.text= element_text(size=15))

plot_all_quints+facet_wrap(~quintile,scales='free_x')

ggsave("./figures/FigS1_all_quintiles-Feb12-21.pdf")

## run all transmission matrix scenarios and plot model dynamics for each, across quintiles

df_all$scenario<-'half within quintile'

output_transmission2=ESIR_transmission_matrix_scenarios('high within quintile')
df_all2=data.frame(output_transmission2[1])
df_all2$scenario<-"high within quintile"

output_transmission3=ESIR_transmission_matrix_scenarios('low within quintile')
df_all3=data.frame(output_transmission3[1])
df_all3$scenario<-"low within quintile"

output_transmission4=ESIR_transmission_matrix_scenarios('Mexico contact survey')
df_all4=data.frame(output_transmission4[1])
df_all4$scenario<-"Mexico contact survey"

output_transmission5=ESIR_transmission_matrix_scenarios('homogeneous mixing')
df_all5=data.frame(output_transmission5[1])
df_all5$scenario<-"homogeneous mixing"

df_all_scenarios=rbind.data.frame(df_all,df_all2,df_all3,df_all4,df_all5)

df_all_scenarios_quint1 <- df_all_scenarios%>%
  group_by(days,scenario) %>%
  summarise(S_stand_mean = mean(S1_stand), 
            I_stand_mean = mean(I1_stand), 
            R_stand_mean = mean(R1_stand),
            D_stand_mean=mean(D1_stand),
            S_conf_lower=quantile(S1_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S1_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I1_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I1_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R1_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R1_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D1_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D1_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(1))

df_all_scenarios_quint2 <- df_all_scenarios%>%
  group_by(days,scenario) %>%
  summarise(S_stand_mean = mean(S2_stand), 
            I_stand_mean = mean(I2_stand), 
            R_stand_mean = mean(R2_stand),
            D_stand_mean=mean(D2_stand),
            S_conf_lower=quantile(S2_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S2_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I2_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I2_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R2_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R2_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D2_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D2_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(2))

df_all_scenarios_quint3 <- df_all_scenarios%>%
  group_by(days,scenario) %>%
  summarise(S_stand_mean = mean(S3_stand), 
            I_stand_mean = mean(I3_stand), 
            R_stand_mean = mean(R3_stand),
            D_stand_mean=mean(D3_stand),
            S_conf_lower=quantile(S3_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S3_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I3_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I3_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R3_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R3_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D3_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D3_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(3))

df_all_scenarios_quint4 <- df_all_scenarios%>%
  group_by(days,scenario) %>%
  summarise(S_stand_mean = mean(S4_stand), 
            I_stand_mean = mean(I4_stand), 
            R_stand_mean = mean(R4_stand),
            D_stand_mean=mean(D4_stand),
            S_conf_lower=quantile(S4_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S4_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I4_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I4_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R4_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R4_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D4_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D4_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(4))

df_all_scenarios_quint5 <- df_all_scenarios%>%
  group_by(days,scenario) %>%
  summarise(S_stand_mean = mean(S5_stand), 
            I_stand_mean = mean(I5_stand), 
            R_stand_mean = mean(R5_stand),
            D_stand_mean=mean(D5_stand),
            S_conf_lower=quantile(S5_stand,c(0.025,0.975))[[1]],
            S_conf_upper=quantile(S5_stand,c(0.025,0.975))[[2]],
            I_conf_lower=quantile(I5_stand,c(0.025,0.975))[[1]],
            I_conf_upper=quantile(I5_stand,c(0.025,0.975))[[2]],
            R_conf_lower=quantile(R5_stand,c(0.025,0.975))[[1]],
            R_conf_upper=quantile(R5_stand,c(0.025,0.975))[[2]],
            D_conf_lower=quantile(D5_stand,c(0.025,0.975))[[1]],
            D_conf_upper=quantile(D5_stand,c(0.025,0.975))[[2]])%>%
  mutate(quintile=rep(5))

df_all_scenarios_all_quint=rbind.data.frame(df_all_scenarios_quint1,df_all_scenarios_quint2,
                                            df_all_scenarios_quint3,df_all_scenarios_quint4,df_all_scenarios_quint5)
colors3=c("S"="blue","I"="red","R"="green", "D"='brown')

plot_all_scenarios_all_quints=ggplot(df_all_scenarios_all_quint,aes(x=days))+
  geom_line(aes(y=I_stand_mean,  col=scenario),linetype='solid')+
  geom_ribbon(aes(ymin=I_conf_lower, ymax=I_conf_upper,fill=scenario), alpha=0.2,linetype='dashed')+
  theme_classic()+xlim(c(0,100))+
  theme(axis.text.x =element_text(size=10),axis.title.x=element_blank(),
        axis.text.y =element_text(size=10),axis.title.y=element_blank(),
        legend.text=element_text(size=10))+ylab("Proportion of quintile population infected")
  

plot_all_scenarios_all_quints+facet_wrap(~quintile,scales='free_x')

ggsave("./figures/FigS2_all_scenarios-Feb12-21.pdf")

