pdf(filename, width=11, height=8.5) 

fig3=ggplot(dat_fig3, aes(x=beta_x, 
                          y=power,
                          group=Approach2, 
                          color=Approach, 
                          linetype=power_out)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)) +
  geom_line(size=1) +
  geom_point(shape=20, size=3, show_guide=FALSE) +
  facet_grid(. ~ ratio2) +
  ggtitle(title_fig3) + 
  labs(x=x_label_beta,y="Power") +
  theme_bw() +
  theme(plot.margin=unit(c(3,2,2,2), "cm"), 
        legend.position="bottom", 
        legend.box="horizontal",
        legend.key=element_rect(colour="white"),
        plot.title=element_text(vjust=6),
        axis.title.x=element_text(vjust=2.1, hjust=-0.125), 
        axis.title.y=element_text(vjust=0)) +
  labs(linetype="Time to Infection", color="Test") +
  guides(linetype=guide_legend(keywidth=1.4)) + 
  geom_text(data=dat_fig3_annotate_txt1, 
            aes(label=dat_fig3_annotate_txt1$mean_n_cases_total, x=dat_fig3_annotate_txt1$beta_x, y=Inf), 
            size=4, 
            color="black", 
            vjust=-2.5) +
  geom_text(data=dat_fig3_annotate_txt1, 
            aes(label=dat_fig3_annotate_txt1$mean_pe, x=dat_fig3_annotate_txt1$beta_x, y=-Inf), 
            size=4, 
            color="black", 
            vjust=4.05) +
  geom_text(data=dat_annotate_lab1,
            aes(label="Mean     ", x=-Inf, y=Inf), 
            size=4, 
            color="black", 
            vjust=-4,
            hjust=1) +
  geom_text(data=dat_annotate_lab1,
            aes(label="No. Cases:", x=-Inf, y=Inf), 
            size=4, 
            color="black", 
            vjust=-2.5,
            hjust=1) +
  geom_text(data=dat_annotate_lab1,
            aes(label="Mean % PE:", x=-Inf, y=-Inf), 
            size=4, 
            color="black", 
            vjust=4.05,
            hjust=1.04) +
  geom_text(data=dat_annotate_lab1,
            aes(label="Power is computed as the proportion of the ~1000 one-sided 2.5%-level tests to reject the null hypotheses.", 
                x=-Inf, y=-Inf), 
            size=3.5, 
            color="black", 
            vjust=11,
            hjust=0) +
  geom_text(data=dat_annotate_lab1,
            aes(label="'Wald/Conc.' uses a daily grid for controls; 'Cox PH Model/Conc.' uses a single record for cases and a daily grid for controls.", 
                x=-Inf, y=-Inf), 
            size=3.5, 
            color="black", 
            vjust=12.5,
            hjust=0) 
gt=ggplot_gtable(ggplot_build(fig3))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)

grid.newpage()
