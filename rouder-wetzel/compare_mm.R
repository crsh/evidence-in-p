load('//users/mmarsman/dropbox/shared/maartenejcollab/NHSTEducationalPsychologySpecialIssue/r/rouder/ej-graph.RData')

one=as.matrix(dat1$one)
pair=as.matrix(dat1$paired)
two=as.matrix(dat1$two)

p=c(one[,6],pair[,8],two[,8])
bf=log10(c(one[,5],pair[,7],two[,7]))

good=(p>.001& p<.15)
mcol=c(rgb(1,0,0,.4),rgb(1,1,0,.4),rgb(.3,1,.8,.4),rgb(0,1,0,.4))

par(cex.main = 1.5, mar = c(5, 5, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(log10(p[good]/2),-bf[good],ylab= "", xlab="")
points(log10(p[good]/2),-bf[good],pch=21,bg=mcol[1],cex=.9)

mtext(expression(log(P[1])),1, line=3.5, cex=1.5, font=2)
mtext(expression(log(BF["01"])), 2, line=3.5, cex=1.5, font=2, las=0)

dev.off()
pdf (file = "//users/mmarsman/dropbox/shared/maartenejcollab/NHSTEducationalPsychologySpecialIssue/tex/figures/rouder.pdf", useDingbats=FALSE)

par(cex.main = 1.5, mar = c(5, 5, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(log10(p[good]/2),-bf[good],ylab= "", xlab="")
points(log10(p[good]/2),-bf[good],pch=21,bg=mcol[1],cex=.9)

mtext(expression(log(P[1])),1, line=3.5, cex=1.5, font=2)
mtext(expression(log(BF["01"])), 2, line=3.5, cex=1.5, font=2, las=0)

dev.off()

