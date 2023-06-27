library("EpiModel") ## Libreria EPIMODEL
library("ndtv")     ## Libreria para la animación
set.seed(12345)     ## seed para los valores aleatorios

## Se crea una población 
num.g1 <- 300
num.g2 <- 250  ## Número total por grupo
nw <- network_initialize(n = num.g1 + num.g2) ## Se inicia la red

group <- rep(1:2, times = c(num.g1, num.g2)) ## Se determina que hay 2 grupos 
nw <- set_vertex_attribute(nw, attrname = "group", value = group) ## 
## Se especifica que la red cuenta con ambos grupos

deg.dist.g1 <- c(0.46, 0.49, 0.04, 0.01) ## Distribución de grado del nodo grupo1
deg.dist.g2 <- c(0.50, 0.35, 0.10, 0.05) ## Distribución de grado del nodo grupo2

check_degdist_bal(num.g1, num.g2,          ## Comprobación -> sigue una distribución poisson
                deg.dist.g1, deg.dist.g2)

formation <- ~edges + degree(0:1, by = "group") + nodematch("group") 
target.stats <- c(180,138, 147, 125, 87.5, 0)
## Estadística de la red.

coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 10)
coef.diss   ## Disolución de los nodos.

est <- netest(nw, formation, target.stats, coef.diss) ## Se crea la red.
dx <- netdx(est, nsims = 5, nsteps = 500, ncores = 5,
            nwstats.formula = ~edges + degree(0:3, by = "group")) 
## Se simula el comportamiento de 5 redes de forma de control. para ver comportamiento
## de nodos y demás
param <- param.net(inf.prob = 0.2, inf.prob.g2 = 0.2,
                   rec.rate = 0.01, rec.rate.g2 = 0.01)
## parametros de la epidemía, probabilidad de infección y de recuperación para ambos grupos
init <- init.net(i.num = 15, i.num.g2 = 8, 
                 r.num = 0, r.num.g2 = 0)
## condiciones iniciales de infectados
control <- control.net(type = "SIS", nsims = 10, nsteps = 500, ncores = 5)
## se especifica el tipo de modelo de epidemia, y los pasos de tiempo.
sim1 <- netsim(est, param, init, control) ## Se inicía la simulación de la epídemia


## ANIMACIÓN  



sex <- get_vertex_attribute(nw, "group")  
sex.shape <- ifelse(sex == 1,50,4)
## lo anterior especifica que cada nodo tenga una figura distinta de acuerdo a su grupo

## captura de frames y creación de archivo.
nw <- get_network(sim1) 
tm <- get_transmat(sim1)
nw <- color_tea(nw, verbose = FALSE)
slice.par <- list(start = 1, end = 110, interval = 2, aggregate.dur = 1, rule = "any")  ## 10 STEPS 
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))
compute.animation(nw, slice.par = slice.par, verbose = TRUE)
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  vertex.sides = sex.shape,
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/epidemiaF2.html"))

## Plots
par(mfrow = c(1, 2), mar = c(2, 0, 2, 0), mgp = c(1, 1, 0))
plot(sim1, type = "network", at = 1, sims = "mean",
     col.status = TRUE, shp.g2 = "square", main = "Red en t = 1")
plot(sim1, type = "network", at = 500, sims = "mean",
     col.status = TRUE,shp.g2 = "square", main = "Red en t = 500")

plot(dx) ## control de la simulación de red



