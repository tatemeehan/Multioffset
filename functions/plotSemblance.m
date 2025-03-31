index = 170;
figure();
imagesc(GPR.D.velocityAxis{ii},GPR.D.stackingTimeAxis{ii},GPR.D.velocityCoherence{ii}{index})
hold on;
plot(GPR.D.stackingVelocity{ii}(:,index.*5),GPR.D.TimeAxis{ii},'k','linewidth',2)
colormap(yetBlack);clim([0.4 1]);
xlabel('Velocity (m/ns)');ylabel('Travel-time (ns)');