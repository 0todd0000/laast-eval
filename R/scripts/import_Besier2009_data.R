
# This script pre-imports the Besier (2009) dataset for faster reading
#
# The "readPatakyData" function from ./LAASTv2/buildFiguresFinal.R
#   runs rather slowly. Here it is replicated with a minor modification:
#   the file name is specified as an input argument. Other than that the
#   "readPatakyData" remains unaltered.
#
# All LAASTv2 code remains unaltered in its original, downloaded form. 



rm( list = ls() )


readPatakyData <- function(file.name) {
	# Read in data used for Pataky et al, 2015 paper
	# Pataky TC, Vanrenterghem J, Robinson MA. Zero- vs. one-dimensional, parametric vs. 
	#	non-parametric, and confidence interval vs. hypothesis testing procedures in 
	#	one-dimensional biomechanical trajectory analysis. J Biomech 2015;48:1277-85. 
	#	doi:10.1016/j.jbiomech.2015.02.051. 
	# These data are available publically as supplemental materials 
	# Besier TF, Fredericson M, Gold GE, Beaupr? GS, Delp SL. Knee Muscle Forces during 
	#	Walking and Running in Patellofemoral Pain Patients and Pain-Free Controls. 
	#	J Biomech 2009;42:898-905.
	# in the file KinematicsEMGmomentarmForces-walking.csv
	 
    # # This next line should be adjusted based on where YOU put these data.
    # kinEMGDat = readLines('KinematicsEMGmomentarmForces-walking.csv');
    kinEMGDat = readLines( file.name );
	
	Subject = c();
	Group = c();
	Gender = c();
	Time = c();
	
	HipFlexion = c();
	KneeFlexion = c();
	AnkleDorsiFlexion = c();
	
	SemiMemEMG = c()
	BicepsFemEMG = c();
	RectusFemEMG = c();
	VastusMedEMG = c();
	VastusLatEMG = c();
	MedGastrocEMG = c();
	LatGastrocEMG = c();
	
	SemiMemMA = c();
	SemiTenMA = c();
	BicepsFemLH_MA = c();
	BicepsFemSH_MA = c();
	RectusFemMA = c();
	VastusMedMA = c();
	VastusIntMA = c();
	VastusLatMA = c();
	MedGastrocMA = c();
	LatGastrocMA = c();
	
	SemiMemF = c();
	SemiTenF = c();
	BicepsFemLH_F = c();
	BicepsFemSH_F = c();
	RectusFemF = c();
	VastusMedF = c();
	VastusIntF = c();
	VastusLatF = c();
	MedGastrocF = c();
	LatGastrocF = c();
	
	
	j = 1;
	for (i in 1:length(kinEMGDat)) {
		ind = i%%106;
		if (ind == 1) {
			tmp = unlist(strsplit(kinEMGDat[i],","));
			subject = tmp[2];
			group = tmp[3];
			if (group != "Control") {
				group = "Pain";	 # replace "Patellofemoral Pain"
			}
			gender = tmp[4];
			t = 1;
		}
		if (ind > 4 & ind <105) {
			tmp = unlist(strsplit(kinEMGDat[i],","));
			
			Subject[j] = subject;
			Group[j] = group;
			Gender[j] = gender;
			Time[j] = t;
			
			HipFlexion[j] = as.numeric(tmp[2])
			KneeFlexion[j] = as.numeric(tmp[3])
			AnkleDorsiFlexion[j] = as.numeric(tmp[4])
			
			SemiMemEMG[j] = as.numeric(tmp[6])
			BicepsFemEMG[j] = as.numeric(tmp[7])
			RectusFemEMG[j] = as.numeric(tmp[8])
			VastusMedEMG[j] = as.numeric(tmp[9])
			VastusLatEMG[j] = as.numeric(tmp[10])
			MedGastrocEMG[j] = as.numeric(tmp[11])
			LatGastrocEMG[j] = as.numeric(tmp[12])
			
			SemiMemMA[j] = as.numeric(tmp[14])
			SemiTenMA[j] = as.numeric(tmp[15])
			BicepsFemLH_MA[j] = as.numeric(tmp[16])
			BicepsFemSH_MA[j] = as.numeric(tmp[17])
			RectusFemMA[j] = as.numeric(tmp[18])
			VastusMedMA[j] = as.numeric(tmp[19])
			VastusIntMA[j] = as.numeric(tmp[20])
			VastusLatMA[j] = as.numeric(tmp[21])
			MedGastrocMA[j] = as.numeric(tmp[22])
			LatGastrocMA[j] = as.numeric(tmp[23])
			
			SemiMemF[j] = as.numeric(tmp[25])
			SemiTenF[j] = as.numeric(tmp[26])
			BicepsFemLH_F[j] = as.numeric(tmp[27])
			BicepsFemSH_F[j] = as.numeric(tmp[28])
			RectusFemF[j] = as.numeric(tmp[29])
			VastusMedF[j] = as.numeric(tmp[30])
			VastusIntF[j] = as.numeric(tmp[31])
			VastusLatF[j] = as.numeric(tmp[32])
			MedGastrocF[j] = as.numeric(tmp[33])
			LatGastrocF[j] = as.numeric(tmp[34])
			
			j = j+1;
			t = t+1;
		}
	}


	
	dataOut = data.frame( Subject, Group, Gender, Time, 
						HipFlexion, KneeFlexion, AnkleDorsiFlexion, 
						SemiMemEMG, BicepsFemEMG, RectusFemEMG, VastusMedEMG, VastusLatEMG, MedGastrocEMG, LatGastrocEMG, 
						SemiMemMA, SemiTenMA, BicepsFemLH_MA, BicepsFemSH_MA, RectusFemMA, VastusMedMA, VastusIntMA, VastusLatMA, MedGastrocMA, LatGastrocMA, 
						SemiMemF, SemiTenF, BicepsFemLH_F, BicepsFemSH_F, RectusFemF, VastusMedF, VastusIntF, VastusLatF, MedGastrocF, LatGastrocF);

	dataOut = dataOut[order(Time),];
	return (dataOut);
}







# assemble directories 
dirREPO      <- dirname( dirname( dirname( sys.frame(1)$ofile ) ) )  # repository path
dirDATA      <- file.path( dirREPO, 'Data' )  # path to the Data directory in this repository
dirR         <- file.path( dirREPO, 'R' )     # path to the R directory in this repository
file.name0   <- file.path( dirR,    'LAASTv2', 'KinematicsEMGmomentarmForces-walking.csv' )
file.name1   <- file.path( dirDATA, 'Besier2009-MedGastrocF.csv' )





# load data:
df      <- readPatakyData( file.name0 )
i1      <- which(df$Group == 'Control')
i2      <- which(df$Group == 'Pain')
# x1      <- df$Time[i1]
# x2      <- df$Time[i2]
y1      <- df$MedGastrocF[i1]
y2      <- df$MedGastrocF[i2]

# reshape in to 2D arrays and stack
dim(y1) <- c(16,100)
dim(y2) <- c(26,100)
y       <- rbind(y1, y2)

# # optionally plot:
# graphics.off()
# matplot( t(y), type='l', lty=1, lwd=0.5, col='black')


# write:
group   <- c(rep(0, 16), rep(1, 26))
a       <- cbind( format( group , digits = 1 ) , format( round( y , 3 ) , nsmall = 3 ) )
write.table(a, file=file.name1, row.names=F, col.names=F, sep=',')



