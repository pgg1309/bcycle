library(readr)
us_ip <- read_csv("us_ip.csv", col_types = cols(Name = col_skip()))
x <- ts(na.omit(us_ip), start = c(1950,1), frequency = 12)

bc2 <- brybos(x,1,3,12)

bc2 <- brybos(window(x, start = c(1970,1)),1,3,12)

plot(x)
points(time(x)[bc$peak], x[bc$peak], col=4)
points(time(x)[bc$trough], x[bc$trough], col=5)


points(time(x)[bc2$peak+240], x[bc2$peak+240], col=4)
points(time(x)[bc2$trough+240], x[bc2$trough+240], col=5)

us_ism <- read_csv("us_ism.csv", col_types = cols(Name = col_skip()))
ism <- ts(na.omit(us_ism), start=c(1950,1), frequency = 12)


bc <- brybos(ism,0,0,12)
plot(ism, lwd=2)
points(time(ism)[bc$peak], ism[bc$peak], col=2, pch=16)
points(time(ism)[bc$trough], ism[bc$trough], col=3, pch=16)

points(time(ism)[bc2$peak], ism[bc2$peak], col=4, pch='#')
points(time(ism)[bc2$trough], ism[bc2$trough], col = 5, pch='#')
