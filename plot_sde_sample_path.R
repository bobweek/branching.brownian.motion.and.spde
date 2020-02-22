require("ggplot2")
require("gridExtra")

smpl_wc = read.csv("/home/bob/Research/WNEE/sample_path_wc.csv")

x_wc = ggplot(smpl_wc,aes(x=time,y=x,group=spp,color=spp))+theme_bw()+theme_minimal()+xlab("")+ylab("x")+
  geom_line(show.legend=F, lwd=0.1)+scale_fill_viridis_c()+ggtitle("Weak Competition")
G_wc = ggplot(smpl_wc,aes(x=time,y=G,group=spp,color=spp))+theme_bw()+theme_minimal()+xlab("")+ylab("G")+
  geom_line(show.legend=F, lwd=0.1)+scale_fill_viridis_c()
N_wc = ggplot(smpl_wc,aes(x=time,y=N,group=spp,color=spp))+theme_bw()+theme_minimal()+xlab("time")+ylab("N")+
  geom_line(show.legend=F, lwd=0.1)+scale_fill_viridis_c()
smpl_plt_wc = grid.arrange(x_wc,G_wc,N_wc,nrow=3)
ggsave("~/Research/WNEE/smpl_plt_wc.png",smpl_plt_wc,width=5,height=2.5)

smpl_mc = read.csv("/home/bob/Research/WNEE/sample_path_mc.csv")

x_mc = ggplot(smpl_mc,aes(x=time,y=x,group=spp,color=spp))+theme_bw()+theme_minimal()+xlab("")+ylab("x")+
  geom_line(show.legend=F, lwd=0.1)+scale_fill_viridis_c()+ggtitle("Moderate Competition")
G_mc = ggplot(smpl_mc,aes(x=time,y=G,group=spp,color=spp))+theme_bw()+theme_minimal()+xlab("")+ylab("G")+
  geom_line(show.legend=F, lwd=0.1)+scale_fill_viridis_c()
N_mc = ggplot(smpl_mc,aes(x=time,y=N,group=spp,color=spp))+theme_bw()+theme_minimal()+xlab("time")+ylab("N")+
  geom_line(show.legend=F, lwd=0.1)+scale_fill_viridis_c()
smpl_plt_mc = grid.arrange(x_mc,G_mc,N_mc,nrow=3)
ggsave("~/Research/WNEE/smpl_plt_mc.png",smpl_plt_mc,width=5,height=2.5)
