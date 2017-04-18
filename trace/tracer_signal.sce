test = read("test.txt",1,512);
resultat_97_lift = read("resultat_97_lift.txt",1,512);
resultat_arm2 = read("resultat_arm_niveau2.txt",1,512);
resultat_arm9 = read("resultat_arm_niveau9.txt",1,512);


subplot(221);
plot(test);
subplot(222);
plot(resultat_97_lift);
subplot(223);
plot(resultat_arm2);
subplot(224);
plot(resultat_arm9);
