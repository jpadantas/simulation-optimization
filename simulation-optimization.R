#Libraries
library(metaheuristicOpt)
library(DEoptim)
library(SPOT)
library(soobench)

#Implement Zakharov function
do_f_zakharov <- function(x) {
  x = 2 * x - 1
  dim = length(x)
  
  ii <- c(1:dim)
  sum1 <- sum(x ^ 2)
  sum2 <- sum(0.5 * ii * x)
  
  y <- sum1 + sum2 ^ 2 + sum2 ^ 4
  return(y)
}

generate_zakharov_function <- function(dimensions)
  soo_function(name="Zakharov",
               id=sprintf("zakharov-%id", dimensions),
               fun=do_f_zakharov,
               dimensions=dimensions,
               lower_bounds=rep(-5, dimensions),
               upper_bounds=rep(10, dimensions),
               best_par=rep(0, dimensions),
               best_value=0)

#Wrap functions for spot() function utilization (dimension 10)

ack10 = wrapFunction(generate_ackley_function(10))
gri10 = wrapFunction(generate_griewank_function(10))
sph10 = wrapFunction(generate_sphere_function(10))
zak10 = wrapFunction(generate_zakharov_function(10))

#Wrap functions for spot() function utilization (dimension 2)
ack2 = wrapFunction(generate_ackley_function(2))
gri2 = wrapFunction(generate_griewank_function(2))
sph2 = wrapFunction(generate_sphere_function(2))
zak2 = wrapFunction(generate_zakharov_function(2))

#Define function domains
ackb = 32.768
grib = 600.0
sphb = 100.0
zakl = -5.0
zaku = 10.0

#Call DEoptim() (400 evals)
ackDE = DEoptim(generate_ackley_function(2),rep(-ackb,2),rep(ackb,2),DEoptim.control(NP = 20,itermax = 20))
griDE = DEoptim(generate_griewank_function(2),rep(-grib,2),rep(grib,2),DEoptim.control(NP = 20,itermax = 20))
sphDE = DEoptim(generate_sphere_function(2),rep(-sphb,2),rep(sphb,2),DEoptim.control(NP = 20,itermax = 20))
zakDE = DEoptim(generate_zakharov_function(2),rep(zakl,2),rep(zaku,2),DEoptim.control(NP = 20,itermax = 20))

#Call DEoptim() (100 evals)
ackDE2 = DEoptim(generate_ackley_function(2),rep(-ackb,2),rep(ackb,2),DEoptim.control(NP = 20,itermax = 5))
griDE2 = DEoptim(generate_griewank_function(2),rep(-grib,2),rep(grib,2),DEoptim.control(NP = 20,itermax = 5))
sphDE2 = DEoptim(generate_sphere_function(2),rep(-sphb,2),rep(sphb,2),DEoptim.control(NP = 20,itermax = 5))
zakDE2 = DEoptim(generate_zakharov_function(2),rep(zakl,2),rep(zaku,2),DEoptim.control(NP = 20,itermax = 5))

#CROSSCHECK
#Call DEoptim() (400 evals)
ackDEO = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20))
griDEO = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20))
sphDEO = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20))
zakDEO = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20))

#Call DEoptim() (100 evals)
ackDEO2 = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5))
griDEO2 = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5))
sphDEO2 = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5))
zakDEO2 = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "DE", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5))

ackDEO2$optimumValue
griDEO2$optimumValue
sphDEO2$optimumValue
zakDEO2$optimumValue

#BHO (400 evals)
ackBHO = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20))
griBHO = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20))
sphBHO = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20))
zakBHO = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20))

#BHO (100 evals)
ackBHO2 = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5))
griBHO2 = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5))
sphBHO2 = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5))
zakBHO2 = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5))

#PSO (400 evals)
ackPSO = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
griPSO = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
sphPSO = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
zakPSO = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

#PSO (100 evals)
ackPSO2 = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
griPSO2 = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
sphPSO2 = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
zakPSO2 = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

#Call spot() (400 evals)
ackSPOT = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=100)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=100)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=100)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=100)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))

#Call spot() (100 evals)
ackSPOT2 = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=25)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT2 = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=25)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT2 = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=25)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT2 = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=25)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))

#Call spot() (20 evals)
ackSPOT3 = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=5)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT3 = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=5)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT3 = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=5)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT3 = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=5)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))

ackSPOT3b = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=5)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT3b = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=5)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT3b = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=5)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT3b = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=5)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))

ackSPOT3c = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=5)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT3c = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=5)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT3c = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=5)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT3c = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=5)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=20,modelControl=list(target="ei"),model=buildKrigingDACE))

#Create results table
results <- t(matrix(c(ackDE$optim$bestval,griDE$optim$bestval,rasDE$optim$bestval,
                      schDE$optim$bestval,sphDE$optim$bestval,zakDE$optim$bestval,
                      ackSPOT$ybest,griSPOT$ybest,rasSPOT$ybest,schSPOT$ybest,
                      sphSPOT$ybest,zakSPOT$ybest,ackDE2$optim$bestval,
                      griDE2$optim$bestval,rasDE2$optim$bestval,schDE2$optim$bestval,
                      sphDE2$optim$bestval,zakDE2$optim$bestval,ackSPOT2$ybest,
                      griSPOT2$ybest,rasSPOT2$ybest,schSPOT2$ybest,sphSPOT2$ybest,
                      zakSPOT2$ybest,ackBHO$optimumValue,griBHO$optimumValue,rasBHO$optimumValue,
                      schBHO$optimumValue,sphBHO$optimumValue,zakBHO$optimumValue,
                      ackBHO2$optimumValue,griBHO2$optimumValue,rasBHO2$optimumValue,
                      schBHO2$optimumValue,sphBHO2$optimumValue,zakBHO2$optimumValue,
                      ackPSO$optimumValue,griPSO$optimumValue,rasPSO$optimumValue,
                      schPSO$optimumValue,sphPSO$optimumValue,zakPSO$optimumValue,
                      ackPSO2$optimumValue,griPSO2$optimumValue,rasPSO2$optimumValue,
                      schPSO2$optimumValue,sphPSO2$optimumValue,zakPSO2$optimumValue),ncol = 6,byrow = TRUE))
colnames(results) <- c("DE 400","SPOT 400","DE 100","SPOT 100","BHO 400","BHO 100","PSO 400","PSO 100")
rownames(results) <- c("Ackley", "Griewank", "Rastrigin", "Schwefel", "Sphere", "Zakharov")

results <- t(matrix(c(ackDE$optim$bestval,griDE$optim$bestval,rasDE$optim$bestval,
                    schDE$optim$bestval,sphDE$optim$bestval,zakDE$optim$bestval,
                    ackSPOT$ybest,griSPOT$ybest,rasSPOT$ybest,schSPOT$ybest,
                    sphSPOT$ybest,zakSPOT$ybest,ackDE2$optim$bestval,
                    griDE2$optim$bestval,rasDE2$optim$bestval,schDE2$optim$bestval,
                    sphDE2$optim$bestval,zakDE2$optim$bestval,ackSPOT2$ybest,
                    griSPOT2$ybest,rasSPOT2$ybest,schSPOT2$ybest,sphSPOT2$ybest,
                    zakSPOT2$ybest,ackBHO$optimumValue,griBHO$optimumValue,rasBHO$optimumValue,
                    schBHO$optimumValue,sphBHO$optimumValue,zakBHO$optimumValue,
                    ackBHO2$optimumValue,griBHO2$optimumValue,rasBHO2$optimumValue,
                    schBHO2$optimumValue,sphBHO2$optimumValue,zakBHO2$optimumValue,
                    ackPSO$optimumValue,griPSO$optimumValue,rasPSO$optimumValue,
                    schPSO$optimumValue,sphPSO$optimumValue,zakPSO$optimumValue,
                    ackPSO2$optimumValue,griPSO2$optimumValue,rasPSO2$optimumValue,
                    schPSO2$optimumValue,sphPSO2$optimumValue,zakPSO2$optimumValue,
                    ackSPOT3$ybest,griSPOT3$ybest,rasSPOT3$ybest,schSPOT3$ybest,
                    sphSPOT3$ybest,zakSPOT3$ybest),ncol = 6,byrow = TRUE))

#Repeat experiment
#Call DEoptim() (400 evals)
ackDEb = DEoptim(generate_ackley_function(2),rep(-ackb,2),rep(ackb,2),DEoptim.control(NP = 20,itermax = 20))
griDEb = DEoptim(generate_griewank_function(2),rep(-grib,2),rep(grib,2),DEoptim.control(NP = 20,itermax = 20))
sphDEb = DEoptim(generate_sphere_function(2),rep(-sphb,2),rep(sphb,2),DEoptim.control(NP = 20,itermax = 20))
zakDEb = DEoptim(generate_zakharov_function(2),rep(zakl,2),rep(zaku,2),DEoptim.control(NP = 20,itermax = 20))

#Call DEoptim() (100 evals)
ackDE2b = DEoptim(generate_ackley_function(2),rep(-ackb,2),rep(ackb,2),DEoptim.control(NP = 20,itermax = 5))
griDE2b = DEoptim(generate_griewank_function(2),rep(-grib,2),rep(grib,2),DEoptim.control(NP = 20,itermax = 5))
sphDE2b = DEoptim(generate_sphere_function(2),rep(-sphb,2),rep(sphb,2),DEoptim.control(NP = 20,itermax = 5))
zakDE2b = DEoptim(generate_zakharov_function(2),rep(zakl,2),rep(zaku,2),DEoptim.control(NP = 20,itermax = 5))

#BHO (400 evals)
ackBHOb = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20))
griBHOb = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20))
sphBHOb = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20))
zakBHOb = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20))

#BHO (100 evals)
ackBHO2b = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5))
griBHO2b = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5))
sphBHO2b = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5))
zakBHO2b = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5))

#PSO (400 evals)
ackPSOb = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
griPSOb = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
sphPSOb = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
zakPSOb = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

#PSO (100 evals)
ackPSO2b = metaOpt(generate_ackley_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
griPSO2b = metaOpt(generate_griewank_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
sphPSO2b = metaOpt(generate_sphere_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
zakPSO2b = metaOpt(generate_zakharov_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

#Call spot() (400 evals)
ackSPOTb = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=100)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOTb = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=100)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOTb = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=100)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOTb = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=100)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))

#Call spot() (100 evals)
ackSPOT2b = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=25)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT2b = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=25)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT2b = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=25)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT2b = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=25)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))

#Call spot() (400 evals)
ackSPOTc = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=100)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOTc = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=100)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOTc = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=100)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOTc = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=100)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))

#Call spot() (100 evals)
ackSPOT2c = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=25)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT2c = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=25)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT2c = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=25)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT2c = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=25)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))

#Create results table
results2 <- t(matrix(c(ackDEb$optim$bestval,griDEb$optim$bestval,rasDEb$optim$bestval,
                    schDEb$optim$bestval,sphDEb$optim$bestval,zakDEb$optim$bestval,
                    ackSPOTb$ybest,griSPOTb$ybest,rasSPOTb$ybest,schSPOTb$ybest,
                    sphSPOTb$ybest,zakSPOTb$ybest,ackDE2b$optim$bestval,
                    griDE2b$optim$bestval,rasDE2b$optim$bestval,schDE2b$optim$bestval,
                    sphDE2b$optim$bestval,zakDE2b$optim$bestval,ackSPOT2b$ybest,
                    griSPOT2b$ybest,rasSPOT2b$ybest,schSPOT2b$ybest,sphSPOT2b$ybest,
                    zakSPOT2b$ybest),ncol = 6,byrow = TRUE))
colnames(results2) <- c("DE 400","SPOT 400","DE 100","SPOT 100")
rownames(results2) <- c("Ackley", "Griewank", "Rastrigin", "Schwefel", "Sphere", "Zakharov")

#Repeat experiment again
#Call DEoptim() (400 evals)
ackDEc = DEoptim(ackley_function(2),rep(-ackb,2),rep(ackb,2),DEoptim.control(NP = 20,itermax = 20))
griDEc = DEoptim(griewank_function(2),rep(-grib,2),rep(grib,2),DEoptim.control(NP = 20,itermax = 20))
sphDEc = DEoptim(sphere_function(2),rep(-sphb,2),rep(sphb,2),DEoptim.control(NP = 20,itermax = 20))
zakDEc = DEoptim(zakharov_function(2),rep(zakl,2),rep(zaku,2),DEoptim.control(NP = 20,itermax = 20))

#Call DEoptim() (100 evals)
ackDE2c = DEoptim(ackley_function(2),rep(-ackb,2),rep(ackb,2),DEoptim.control(NP = 20,itermax = 5))
griDE2c = DEoptim(griewank_function(2),rep(-grib,2),rep(grib,2),DEoptim.control(NP = 20,itermax = 5))
sphDE2c = DEoptim(sphere_function(2),rep(-sphb,2),rep(sphb,2),DEoptim.control(NP = 20,itermax = 5))
zakDE2c = DEoptim(zakharov_function(2),rep(zakl,2),rep(zaku,2),DEoptim.control(NP = 20,itermax = 5))

#BHO (400 evals)
ackBHOc = metaOpt(ackley_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20))
griBHOc = metaOpt(griewank_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20))
sphBHOc = metaOpt(sphere_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20))
zakBHOc = metaOpt(zakharov_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20))

#BHO (100 evals)
ackBHO2c = metaOpt(ackley_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5))
griBHO2c = metaOpt(griewank_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5))
sphBHO2c = metaOpt(sphere_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5))
zakBHO2c = metaOpt(zakharov_function(2),optimType = "MIN", algorithm = "BHO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5))

#PSO (400 evals)
ackPSOc = metaOpt(ackley_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
griPSOc = metaOpt(griewank_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
sphPSOc = metaOpt(sphere_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
zakPSOc = metaOpt(zakharov_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=20, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

#PSO (100 evals)
ackPSO2c = metaOpt(ackley_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-ackb,ackb),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
griPSO2c = metaOpt(griewank_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-grib,grib),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
sphPSO2c = metaOpt(sphere_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(-sphb,sphb),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))
zakPSO2c = metaOpt(zakharov_function(2),optimType = "MIN", algorithm = "PSO", numVar = 2, rangeVar = matrix(c(zakl,zaku),nrow=2), control = list(numPopulation = 20, maxIter=5, Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

#Call spot() (400 evals)
ackSPOTc = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=100)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOTc = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=100)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOTc = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=100)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOTc = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=100)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=400,modelControl=list(target="ei"),model=buildKrigingDACE))

#Call spot() (100 evals)
ackSPOT2c = spot(designLHD(,rep(-ackb,2),rep(ackb,2),control=list(size=25)),ack2,rep(-ackb,2),rep(ackb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
griSPOT2c = spot(designLHD(,rep(-grib,2),rep(grib,2),control=list(size=25)),gri2,rep(-grib,2),rep(grib,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
sphSPOT2c = spot(designLHD(,rep(-sphb,2),rep(sphb,2),control=list(size=25)),sph2,rep(-sphb,2),rep(sphb,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))
zakSPOT2c = spot(designLHD(,rep(zakl,2),rep(zaku,2),control=list(size=25)),zak2,rep(zakl,2),rep(zaku,2),control=list(funEvals=100,modelControl=list(target="ei"),model=buildKrigingDACE))

#Create results table
results3 <- t(matrix(c(ackDEc$optim$bestval,griDEc$optim$bestval,rasDEc$optim$bestval,
                    schDEc$optim$bestval,sphDEc$optim$bestval,zakDEc$optim$bestval,
                    ackSPOTc$ybest,griSPOTc$ybest,rasSPOTc$ybest,schSPOTc$ybest,
                    sphSPOTc$ybest,zakSPOTc$ybest,ackDE2c$optim$bestval,
                    griDE2c$optim$bestval,rasDE2c$optim$bestval,schDE2c$optim$bestval,
                    sphDE2c$optim$bestval,zakDE2c$optim$bestval,ackSPOT2c$ybest,
                    griSPOT2c$ybest,rasSPOT2c$ybest,schSPOT2c$ybest,sphSPOT2c$ybest,
                    zakSPOT2c$ybest),ncol = 6,byrow = TRUE))
colnames(results3) <- c("DE 400","SPOT 400","DE 100","SPOT 100")
rownames(results3) <- c("Ackley", "Griewank", "Rastrigin", "Schwefel", "Sphere", "Zakharov")


results1 <- t(matrix(c(
ackDE$optim$bestval,
griDE$optim$bestval,
rasDE$optim$bestval,
schDE$optim$bestval,
sphDE$optim$bestval,
zakDE$optim$bestval,
ackBHO$optimumValue,
griBHO$optimumValue,
rasBHO$optimumValue,
schBHO$optimumValue,
sphBHO$optimumValue,
zakBHO$optimumValue,
ackPSO$optimumValue,
griPSO$optimumValue,
rasPSO$optimumValue,
schPSO$optimumValue,
sphPSO$optimumValue,
zakPSO$optimumValue,
ackSPOT$ybest,
griSPOT$ybest,
rasSPOT$ybest,
schSPOT$ybest,
sphSPOT$ybest,
zakSPOT$ybest,
ackDE2$optim$bestval,
griDE2$optim$bestval,
rasDE2$optim$bestval,
schDE2$optim$bestval,
sphDE2$optim$bestval,
zakDE2$optim$bestval,
ackBHO2$optimumValue,
griBHO2$optimumValue,
rasBHO2$optimumValue,
schBHO2$optimumValue,
sphBHO2$optimumValue,
zakBHO2$optimumValue,
ackPSO2$optimumValue,
griPSO2$optimumValue,
rasPSO2$optimumValue,
schPSO2$optimumValue,
sphPSO2$optimumValue,
zakPSO2$optimumValue,
ackSPOT2$ybest,
griSPOT2$ybest,
rasSPOT2$ybest,
schSPOT2$ybest,
sphSPOT2$ybest,
zakSPOT2$ybest,
ackSPOT3$ybest,
griSPOT3$ybest,
rasSPOT3$ybest,
schSPOT3$ybest,
sphSPOT3$ybest,
zakSPOT3$ybest),ncol = 6,byrow = TRUE))

results2 <- t(matrix(c(
  ackDEb$optim$bestval,
  griDEb$optim$bestval,
  rasDEb$optim$bestval,
  schDEb$optim$bestval,
  sphDEb$optim$bestval,
  zakDEb$optim$bestval,
  ackBHO$optimumValue,
  griBHOb$optimumValue,
  rasBHOb$optimumValue,
  schBHOb$optimumValue,
  sphBHOb$optimumValue,
  zakBHOb$optimumValue,
  ackPSOb$optimumValue,
  griPSOb$optimumValue,
  rasPSOb$optimumValue,
  schPSOb$optimumValue,
  sphPSOb$optimumValue,
  zakPSOb$optimumValue,
  ackSPOTb$ybest,
  griSPOTb$ybest,
  rasSPOTb$ybest,
  schSPOTb$ybest,
  sphSPOTb$ybest,
  zakSPOTb$ybest,
  ackDE2b$optim$bestval,
  griDE2b$optim$bestval,
  rasDE2b$optim$bestval,
  schDE2b$optim$bestval,
  sphDE2b$optim$bestval,
  zakDE2b$optim$bestval,
  ackBHO2b$optimumValue,
  griBHO2b$optimumValue,
  rasBHO2b$optimumValue,
  schBHO2b$optimumValue,
  sphBHO2b$optimumValue,
  zakBHO2b$optimumValue,
  ackPSO2b$optimumValue,
  griPSO2b$optimumValue,
  rasPSO2b$optimumValue,
  schPSO2b$optimumValue,
  sphPSO2b$optimumValue,
  zakPSO2b$optimumValue,
  ackSPOT2b$ybest,
  griSPOT2b$ybest,
  rasSPOT2b$ybest,
  schSPOT2b$ybest,
  sphSPOT2b$ybest,
  zakSPOT2b$ybest,
  ackSPOT3b$ybest,
  griSPOT3b$ybest,
  rasSPOT3b$ybest,
  schSPOT3b$ybest,
  sphSPOT3b$ybest,
  zakSPOT3b$ybest),ncol = 6,byrow = TRUE))

results3 <- t(matrix(c(
  ackDEc$optim$bestval,
  griDEc$optim$bestval,
  rasDEc$optim$bestval,
  schDEc$optim$bestval,
  sphDEc$optim$bestval,
  zakDEc$optim$bestval,
  ackBHOc$optimumValue,
  griBHOc$optimumValue,
  rasBHOc$optimumValue,
  schBHOc$optimumValue,
  sphBHOc$optimumValue,
  zakBHOc$optimumValue,
  ackPSOc$optimumValue,
  griPSOc$optimumValue,
  rasPSOc$optimumValue,
  schPSOc$optimumValue,
  sphPSOc$optimumValue,
  zakPSOc$optimumValue,
  ackSPOTc$ybest,
  griSPOTc$ybest,
  rasSPOTc$ybest,
  schSPOTc$ybest,
  sphSPOTc$ybest,
  zakSPOTc$ybest,
  ackDE2c$optim$bestval,
  griDE2c$optim$bestval,
  rasDE2c$optim$bestval,
  schDE2c$optim$bestval,
  sphDE2c$optim$bestval,
  zakDE2c$optim$bestval,
  ackBHO2c$optimumValue,
  griBHO2c$optimumValue,
  rasBHO2c$optimumValue,
  schBHO2c$optimumValue,
  sphBHO2c$optimumValue,
  zakBHO2c$optimumValue,
  ackPSO2c$optimumValue,
  griPSO2c$optimumValue,
  rasPSO2c$optimumValue,
  schPSO2c$optimumValue,
  sphPSO2c$optimumValue,
  zakPSO2$optimumValue,
  ackSPOT2c$ybest,
  griSPOT2c$ybest,
  rasSPOT2c$ybest,
  schSPOT2c$ybest,
  sphSPOT2c$ybest,
  zakSPOT2c$ybest,
  ackSPOT3c$ybest,
  griSPOT3c$ybest,
  rasSPOT3c$ybest,
  schSPOT3c$ybest,
  sphSPOT3c$ybest,
  zakSPOT3c$ybest),ncol = 6,byrow = TRUE))

#Average
final <- (results1+results2+results3)/3
final
