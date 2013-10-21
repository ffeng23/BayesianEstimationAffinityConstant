
dt_simN<-read.table("simulation_attach_noise.txt", header=T, sep="\t");
plot(c(0,max(dt_simN[,1]+50)), c(min(dt_simN[,2]),max(dt_simN[,2])), main="simulated", col=2, type="n");

lines(dt_simN[,1], dt_simN[,2], col=8);
dt_sim<-read.table("simulation_attach.txt", header=T, sep="\t");
lines(dt_sim[,1], dt_sim[,2], col=6);




#####for displaying the testing file for linear regression

#simulation detached
dt_simN_de<-read.table("simulation_detach_noise.txt", header=T, sep="\t");
dt_simN_de_short<-read.table("simulation_detach_noise_short.txt", header=T, sep="\t");
dt_simN_de_smoothed<-read.table("simulation_detach_noise_smoothed.txt", header=T, sep="\t");
dt_simN_de_slope<-read.table("simulation_detach_noise_slope.txt", header=T, sep="\t");

invDt_dt=1/dt_simN_de_slope[,2];
invR=1/dt_simN_de_smoothed[,2]

jpeg(filename="plot.jpg", width = 1200, height = 1200, pointsize = 16)
 
plot(c(0,max(dt_simN_de[,1]+50)), c(min(dt_simN_de[,2]),max(dt_simN_de[,2])), main="simulated", col=2, type="l");
lines(dt_simN_de[,1], dt_simN_de[,2], col=8);
lines(dt_simN_de_short[,1], dt_simN_de_short[,2], col=9);
lines(dt_simN_de_smoothed[,1], dt_simN_de_smoothed[,2], col=2);

plot(dt_simN_de_slope[,1], dt_simN_de_slope[,2], main="slope vs time", col=2, type="p");
plot(dt_simN_de_smoothed[,2], dt_simN_de_slope[,2], main="slope vs response", col=2, type="p");
plot(invR, invDt_dt, main="inv dr/dt vs inv RU", col=2, type="p");

par(op)
dev.off();
regOut<-lm(invDt_dt~ invR)
summary(regOut);
n=1;
regOut<-lm(invDt_dt[c(n:330)]~ invR[c(n:330)])
summary(regOut);




setwd("E:/MSI_software/SPR/simulation/MH/MH 2000000");
setwd("E:/MSI_software/SPR/step0.005")
setwd("E:/MSI_software/SPR/BayesianEstimationAffinityConstant/Testing/bin/Debug")
setwd("E:/MSI_software/SPR/BayesianEstimationAffinityConstant/bayestimateCon/bin/Debug")
dt_mc<-read.table("MCMC_run.txt", header=T, sep="\t", );

dt_mc<-dt_mc[c(40000:length(dt_mc[,1])),];
jpeg(filename="trace.jpg", width = 800, height = 2400, pointsize = 16)
op<-par( mfrow = c( 7, 1) )
plot(dt_mc[,1], dt_mc$ka, main="ka", col=2, type="l");
plot(dt_mc[,1], dt_mc$kd, main="kd", col=2, type="l");
plot(dt_mc[,1], dt_mc$kM, main="kM", col=2, type="l");
plot(dt_mc[,1], dt_mc$conc, main="conc", col=2, type="l");
plot(dt_mc[,1], dt_mc$Rmax, main="Rmax", col=2, type="l");
plot(dt_mc[,1], dt_mc$R0, main="R0", col=2, type="l");
plot(dt_mc[,1], dt_mc$curLLD, main="LLD", col=2, type="l");



par(op)
dev.off();

kc<-dt_mc$ka*dt_mc$conc;

mean(dt_mc$ka)
mean(dt_mc$kd)
mean(dt_mc$kM)
mean(dt_mc$conc)
mean(dt_mc$Rmax)
mean(dt_mc$R0)
mean(dt_mc$Sigma)
mean(kc);

plot(dt_mc[,1], kc, main="ka*conc", col=2, type="l");

k1<- mean(dt_mc$conc)* mean(dt_mc$ka)/ mean(dt_mc$kM)
k2<-1.5E3*1.9E-6/1E6

############displaying the fitting results for comparison
jpeg(filename="fitting.jpg", width = 800, height = 2400, pointsize = 16)
op<-par( mfrow = c( 2, 1) )
dt_simN<-read.table("simulation_attach_noise.txt", header=T, sep="\t");
plot(c(0,max(dt_simN[,1]+600)), c(min(dt_simN[,2]),max(dt_simN[,2])), main="simulated", col=2, type="n");

lines(dt_simN[,1], dt_simN[,2], col=8);
dt_simN<-read.table("simulation_detach_noise.txt", header=T, sep="\t");
lines(dt_simN[,1]+400, dt_simN[,2], col=8);

dt_sim<-read.table("simulation_attach.txt", header=T, sep="\t");
#lines(dt_sim[,1], dt_sim[,2], col=2);
points(dt_sim[,1], dt_sim[,2], col=2);
dt_sim<-read.table("simulation_detach.txt", header=T, sep="\t");
#lines(dt_sim[,1], dt_sim[,2], col=2);
points(dt_sim[,1]+400, dt_sim[,2], col=2);

###now for fitting, as means
dt_simMean<-read.table("simulationMean.txt", header=T, sep="\t");
points(dt_simMean[,1], dt_simMean[,2], col=3);

dt_simMean<-read.table("simulationMean_detach.txt", header=T, sep="\t");
points(dt_simMean[,1]+400, dt_simMean[,2], col=3);


###now for fitting, as last set
dt_simLast<-read.table("simulationLast.txt", header=T, sep="\t");
lines(dt_simLast[,1], dt_simLast[,2], col=4, lty=2);

dt_simLast<-read.table("simulationLast_detach.txt", header=T, sep="\t");
lines(dt_simLast[,1]+400, dt_simLast[,2], col=4, lty=2);


lines(t, R[c(1:(length(R)-1))], col=5);
dev.off();


