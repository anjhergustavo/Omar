ss=getenv('computername');
compname=str2num([ss(length(ss)-1) ss(length(ss))])
fringeprocessing(75*(compname-1)+1,75*compname,'HMD85')
disp('HMD85')
fringeprocessing(75*(compname-1)+1,75*compname,'HMD75')
disp('T60V')
fringeprocessing(75*(compname-1)+1,75*compname,'T60V')
disp('T90v')
fringeprocessing(75*(compname-1)+1,75*compname,'T90V')
disp('AMB')
fringeprocessing(75*(compname-1)+1,75*compname,'AMB')
disp('T120V')
fringeprocessing(75*(compname-1)+1,75*compname,'T120V')
disp('AIRE3.5')
fringeprocessing(75*(compname-1)+1,75*compname,'AIRE3_5')
disp('AIRE4.5')
fringeprocessing(75*(compname-1)+1,75*compname,'AIRE4_5')
disp('AIRE5.5')
fringeprocessing(75*(compname-1)+1,75*compname,'AIRE5_5')
disp('AIRE6.5')
fringeprocessing(75*(compname-1)+1,75*compname,'AIRE6_5')