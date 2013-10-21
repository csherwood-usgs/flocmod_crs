% plot_flocmod - Script to look at FlocMod results
load('flocmod.dat')
figure(1); clf
subplot(211)
plot(flocmod(:,1)./60,flocmod(:,2))
ylabel('Shear rate G (s-1)')
subplot(212)
plot(flocmod(:,1)./60,flocmod(:,3))
ylabel('D_{50} (um)')
xlabel('Time (min)')