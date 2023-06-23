naturemed_reqs <- function() 
{
  # call this before doing your plots
  reqs <<- theme(
    axis.title = element_text(size = 7),
    axis.text = element_text(size=5, family='sans'),
    text=element_text(size=7,family='sans'),
    legend.title=element_text(size=7, family='sans'),
    legend.text=element_text(size=6, family='sans'),
    legend.key.size = unit(10, 'pt'),
    strip.text = element_text(size = 5),
    plot.tag = element_text(size=8, face='bold', color = "black"),
    plot.tag.position = c(0, 1),
    NULL
  )
}

ggarrange_nature <- function(
  ...,
  plotlist = NULL,
  ncol = NULL,
  nrow = NULL,
  labels = NULL,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5, align = c("none", "h", "v", "hv"),
  widths = 1,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL, 
  add_reqs=TRUE
){
  if(add_reqs)
    naturemed_reqs()
  
  plots <- c(list(...), plotlist)
  
  if(add_reqs)
    plots <- lapply(plots, function(p){p + reqs})
  
  out <- ggarrange(plotlist = plots,
                   ncol = ncol,
                   nrow = nrow,
                   labels = labels,
                   label.x = label.x,
                   label.y = label.y,
                   hjust = hjust,
                   vjust = vjust,
                   font.label = list(size = 8, color = "black", face = "bold", family = 'sans'),
                   align = align,
                   widths = widths,
                   heights = heights,
                   legend = legend,
                   common.legend = common.legend,
                   legend.grob = legend.grob)
  return(out)
}

ggsave_nature <- function(filename, p, w=18,h=24, add_reqs=TRUE)
{
  # check size
  tmp <- sort(c(w,h))
  if(tmp[1] > 18 | tmp[2] > 24)
    warning('Plot is bigger than allowed for EDFs. Maximum size is 18cm x 24cm\n')
  if( tmp[1] < 10)
    warning('w and h represent cm units, not inches. Are you sure you want to save such a small plot?\n')
  
  # apply changes: for the moment only works for simple plots, but not for ggarrange.
  # Let's see if anyone answers this:
  # https://stackoverflow.com/questions/74379207/simultaneously-applying-same-modification-to-all-ggarrange-subplots
  if(add_reqs)
  {
    naturemed_reqs()
    p <- p + reqs
  }
  
  # save
  ggsave(filename=filename, plot=p, width=w, height=h, units='cm', dpi=310)

    # return cmd command to open in viewer
    prog <- data.table::fcase(
        filename %like% '.png$', 'gthumb', 
        filename %like% '.pdf$', 'zathura', 
        default = "xdg-open"
    )
    cmd <- sprintf("! %s %s &", prog , filename)
    return(cmd)
}

rotate_x_axis <- function(angle){
    theme( axis.text.x = element_text(angle=angle, hjust = 1, vjust = 1) )
}

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

