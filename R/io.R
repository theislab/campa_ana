ggsave_cairo <- function(filename,plot = last_plot(),width=NA,height=NA,units="cm",dpi=300) {
  require(Cairo)
  CairoPNG(filename=filename,height=height,width=width,units=units,res=dpi)
  print(plot)
  dev.off()
}
