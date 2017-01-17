# Latent Growth Mixture Model
Doug Tommet  
`r Sys.Date()`  






```r
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(MplusAutomation)


# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
```



# Create Data
Since we don't have original data to use in this example, we are going to simulate some data using Mplus' Monte Carlo option.

MplusAutomation doesn't work with the Monte Carlo option, so we created the Mplus input file separately and run it with the runModels() command. *Note:* This will run all **.inp** files in the working directory. 

The model simulates longitudinal data over six time points for 500 subjects.  There are four latent classes.


```r
file.remove("gmm001.inp")
```

```
## [1] FALSE
```

```r
runModels()
```

```
## 
## Running model: monte-carlo.inp 
## System command: C:\Windows\system32\cmd.exe /c cd "S:\work\amphibian-lgmm" && "Mplus" "monte-carlo.inp"
```






```r
data <- getSavedata_Data("monte-carlo.out") 
colnames(data) <- tolower(names(data))  

# data <- read.fwf(file="mcdata.dat", header=FALSE, 
#                  widths = c(rep(13, 6), 3), 
#                  col.names = c("y1", "y2", "y3", "y4", "y5", "y6", "c"))

data <- data %>%
  mutate(id = row_number(),
         subsample = runif(nrow(.)))


data.long <- data %>% 
  gather(foo, y, y1:y6) %>%
  mutate(time = as.numeric(str_sub(foo, 2))) %>%
  arrange(id, time) %>%
  select(-foo)

head(data)
```

```
##          y1        y2        y3        y4        y5        y6 c id
## 1 10.648888 11.877818 12.655686 15.146736 17.713483 17.588851 4  1
## 2  7.950499  9.425913  8.590537 10.074895  9.289337 10.080144 2  2
## 3 12.371036 13.421617 17.478057 19.415944 22.843179 26.829369 4  3
## 4  1.968726  1.901051  0.513132  1.417859  3.821579  7.512429 3  4
## 5 10.623194 12.873786 15.352928 17.648866 20.551219 23.809973 4  5
## 6  9.836673 12.727578 15.644474 15.997809 19.392901 20.507929 4  6
##    subsample
## 1 0.36187559
## 2 0.19509423
## 3 0.47261990
## 4 0.37084734
## 5 0.95150092
## 6 0.03593966
```

```r
head(data.long, 18)
```

```
##    c id subsample         y time
## 1  4  1 0.3618756 10.648888    1
## 2  4  1 0.3618756 11.877818    2
## 3  4  1 0.3618756 12.655686    3
## 4  4  1 0.3618756 15.146736    4
## 5  4  1 0.3618756 17.713483    5
## 6  4  1 0.3618756 17.588851    6
## 7  2  2 0.1950942  7.950499    1
## 8  2  2 0.1950942  9.425913    2
## 9  2  2 0.1950942  8.590537    3
## 10 2  2 0.1950942 10.074895    4
## 11 2  2 0.1950942  9.289337    5
## 12 2  2 0.1950942 10.080144    6
## 13 4  3 0.4726199 12.371036    1
## 14 4  3 0.4726199 13.421617    2
## 15 4  3 0.4726199 17.478057    3
## 16 4  3 0.4726199 19.415944    4
## 17 4  3 0.4726199 22.843179    5
## 18 4  3 0.4726199 26.829369    6
```






```r
data.long %>%
  filter(subsample <.2) %>%
  ggplot(aes(x=time, y=y, group=id)) +
    geom_line() +
    geom_smooth(aes(group = 1), size=2)
```

![](000-master_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
# In real life we can't fit this plot because we don't know the true latent class membership.
data.long %>%
  filter(subsample <.25) %>%
  ggplot(aes(x=time, y=y, group=id)) +
    geom_line() +
    facet_grid(. ~c)
```

![](000-master_files/figure-html/unnamed-chunk-10-2.png)<!-- -->





```r
gmm.model.1 <- " %OVERALL% 
                 i s | y1@0 y2@1 y3@2 y4@3 y5@4 y6@5 ;
               "

gmm.body <- mplusObject(
  ANALYSIS = "type=mixture; starts = 100 20;",
  VARIABLE = "classes = c(4); idvariable = id;",
  MODEL = gmm.model.1,
  OUTPUT = "",
  SAVEDATA = "save = cprob; file=model1prob.dat;",
  usevariables = c("id", "y1", "y2", "y3", "y4", "y5", "y6"),
  rdata = data
)
gmm.fit <- mplusModeler(gmm.body,
                          modelout = "gmm001.inp",
                          run = TRUE)
```

```
## 
## Running model: gmm001.inp 
## System command: C:\Windows\system32\cmd.exe /c cd "S:\work\amphibian-lgmm" && "Mplus" "gmm001.inp" 
## Reading model:  gmm001.out
```

```r
param <- extractModelParameters("gmm001.out")$unstandardized %>% 
  filter(grepl("Means", paramHeader)) %>%
  filter(!grepl("Categorical.Latent.Variables", LatentClass))

foo <- getSavedata_Data("gmm001.out") %>%
  rename(model1.c = C) %>%
  select(ID, model1.c) 
colnames(foo) <- tolower(names(foo)) 

data <- data %>%
  full_join(foo, by= "id" )


data.long <- data %>% 
  gather(foo, y, y1:y6) %>%
  mutate(time = as.numeric(str_sub(foo, 2))) %>%
  arrange(id, time) %>%
  select(-foo)
```






```r
param.matrix <- param %>%
  select(param, est, LatentClass) %>%
  spread(param, est) %>%
  select(-LatentClass)


fun.1 <- function(x) 11.707 + 0.171*x
fun.2 <- function(x)  2.282 + 3.486*x
fun.3 <- function(x) 17.348 - 0.059*x
fun.4 <- function(x)  9.261 + 1.923*x

myplot <- function(x, f) {
  data.long %>%
  filter(subsample <.25) %>%
  filter(model1.c==x) %>%
    ggplot(aes(x=time, y=y, group=id)) +
      scale_y_continuous(limits = c(0, 30)) +
      geom_line() + 
      stat_function(fun = f, color="blue", size=2) +
      ggtitle(paste("class", x))
  
}
plot1 <- myplot(1, fun.1)
plot2 <- myplot(2, fun.2)
plot3 <- myplot(3, fun.3)
plot4 <- myplot(4, fun.4)

multiplot(plotlist = list(plot1, plot2, plot3, plot4), cols = 2)
```

![](000-master_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

