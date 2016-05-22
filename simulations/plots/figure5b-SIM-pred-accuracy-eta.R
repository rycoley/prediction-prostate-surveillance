### NAME AND EMAIL REDACTED
### Code for "A Bayesian Hierarchical Model for Prediction of Latent Health States from Multiple Data Sources with Application to Active Surveillance of Prostate Cancer"
# Plot predictive accuracy for Figure 5(b)


#Get results from Table 15 (appendix), from appendix-Table15-SIM-pred-accuracy.R


pdf("simulations/plots/figure5b-SIM-pred-accuracy-eta.pdf", width=7, height=5)
par(mfrow=c(2,1), mar=c(3,8,2,1), las=1)

plot(0, type="n", xlim=c(0.5,1), ylim=c(0.5,5.5), xlab="", ylab="", yaxt="n", main="Simulation: True State Observed Post-Surgery")
mtext("AUC (95% Interval)", side=1, line=2 )
text(x=rep(0.48, 5), y=c(5:1), labels=c("Biopsy, Surgery IOP", "Biopsy IOP only", "Surgery IOP only", "Unadjusted", "Logistic"), srt=0, pos=2, xpd=TRUE)

lines(x=c(0.5, 0.83), y=c(5,5), lty="dotted")
points(x=0.83, y=5, pch=19, cex=2)
lines(x=c(0.77, 0.88), y=c(5,5), lty="solid", lwd=3)

lines(x=c(0.5, 0.81), y=c(4,4), lty="dotted")
points(x=0.81, y=4, pch=19, cex=2)
lines(x=c(0.76, 0.86), y=c(4,4), lty="solid", lwd=3)

lines(x=c(0.5, 0.81), y=c(3,3), lty="dotted")
points(x=0.81, y=3, pch=19, cex=2)
lines(x=c(0.76, 0.86), y=c(3,3), lty="solid", lwd=3)

lines(x=c(0.5, 0.80), y=c(2,2), lty="dotted")
points(x=0.80, y=2, pch=19, cex=2)
lines(x=c(0.74, 0.85), y=c(2,2), lty="solid", lwd=3)

lines(x=c(0.5, 0.77), y=c(1,1), lty="dotted")
points(x=0.77, y=1, pch=19, cex=2)
lines(x=c(0.70, 0.83), y=c(1,1), lty="solid", lwd=3)


plot(0, type="n", xlim=c(0.5,1), ylim=c(0.5,5.5), xlab="", ylab="", yaxt="n", main="Simulation: True State Unobserved")
mtext("AUC (95% Interval)", side=1, line=2 )
text(x=rep(0.48, 5), y=c(5:1), labels=c("Biopsy, Surgery IOP", "Biopsy IOP only", "Surgery IOP only", "Unadjusted", "Logistic"), srt=0, pos=2, xpd=TRUE)

lines(x=c(0.5, 0.77), y=c(5,5), lty="dotted")
points(x=0.77, y=5, pch=19, cex=2)
lines(x=c(0.72, 0.81), y=c(5,5), lty="solid", lwd=3)

lines(x=c(0.5, 0.74), y=c(4,4), lty="dotted")
points(x=0.74, y=4, pch=19, cex=2)
lines(x=c(0.69, 0.78), y=c(4,4), lty="solid", lwd=3)

lines(x=c(0.5, 0.73), y=c(3,3), lty="dotted")
points(x=0.73, y=3, pch=19, cex=2)
lines(x=c(0.67, 0.78), y=c(3,3), lty="solid", lwd=3)

lines(x=c(0.5, 0.71), y=c(2,2), lty="dotted")
points(x=0.71, y=2, pch=19, cex=2)
lines(x=c(0.66, 0.75), y=c(2,2), lty="solid", lwd=3)

lines(x=c(0.5, 0.68), y=c(1,1), lty="dotted")
points(x=0.68, y=1, pch=19, cex=2)
lines(x=c(0.63, 0.73), y=c(1,1), lty="solid", lwd=3)

dev.off()

