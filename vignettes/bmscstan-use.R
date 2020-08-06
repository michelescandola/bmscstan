## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.height= 6, fig.width= 6---------------------------------------------
library(bmscstan)

data(BSE)

str(data.pt)

str(data.ctrl)

ggplot(data.pt, aes(y = RT, x = Body.District:Side , fill = Congruency))+
  geom_boxplot()

ggplot(data.ctrl, aes(y = RT, x = Body.District:Side , fill = Congruency))+
  geom_boxplot()+
  facet_wrap( ~ ID , ncol = 4)

## ---- fig.show='hold'---------------------------------------------------------
qqnorm(data.ctrl$RT, main = "Controls")
qqline(data.ctrl$RT)

qqnorm(data.pt$RT, main = "Single Case")
qqline(data.pt$RT)

## ---- fig.show='hold'---------------------------------------------------------
out <- boxplot.stats( data.ctrl$RT )$out
data.ctrl <- droplevels( data.ctrl[ !data.ctrl$RT %in% out , ] )

out <- boxplot.stats( data.pt$RT )$out
data.pt <- droplevels( data.pt[ !data.pt$RT %in% out , ] )

qqnorm(data.ctrl$RT, main = "Controls")
qqline(data.ctrl$RT)

qqnorm(data.pt$RT, main = "Single Case")
qqline(data.pt$RT)

## -----------------------------------------------------------------------------
contrasts( data.ctrl$Side )          <- contr.treatment( n = 2 )
contrasts( data.ctrl$Congruency )    <- contr.treatment( n = 2 )
contrasts( data.ctrl$Body.District ) <- contr.treatment( n = 2 )

contrasts( data.pt$Side )            <- contr.treatment( n = 2 )
contrasts( data.pt$Congruency )      <- contr.treatment( n = 2 )
contrasts( data.pt$Body.District )   <- contr.treatment( n = 2 )

## -----------------------------------------------------------------------------
data.ctrl$BD_ID <- interaction( data.ctrl$Body.District , data.ctrl$ID )

