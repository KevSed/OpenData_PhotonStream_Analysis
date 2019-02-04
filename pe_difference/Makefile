# INPUTFILES
CRAB_FILE=crab
GAMMA_FILE1=gamma
GAMMA_FILE2=/gustav
PROTON_FILE1=proton
PROTON_FILE2=/uwe

P_RUN1=066
P_RUN2=009
G_RUN1=010
G_RUN2=910
C_RUN1=04
C_RUN2=162

DL1=/net/big-tank/POOL/projects/fact/photon-stream/facttools


all: dbscan thresholds facttools

hist_dbscan:
#	TEXINPUTS=../header_matplotlib.tex: \
#	MATPLOTLIBRC=../matplotlibrc_half \
#	python plot_histograms.py DBSCAN $(PROTON_FILE1) $(PROTON_FILE2) $(P_RUN1)$(P_RUN2) delta --number=19
#
#	TEXINPUTS=../header_matplotlib.tex: \
#	MATPLOTLIBRC=../matplotlibrc_half \
#	python plot_histograms.py DBSCAN $(GAMMA_FILE1) $(GAMMA_FILE2) $(G_RUN1)$(G_RUN2) delta --number=120

	TEXINPUTS=../header_matplotlib.tex: \
	MATPLOTLIBRC=../matplotlibrc_full \
	python plot_histograms.py DBSCAN $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) pe

image_dbscan:
	python plot_cleaned_images.py DBSCAN $(PROTON_FILE1) $(PROTON_FILE2) $(P_RUN1)$(P_RUN2) $(FEAT) --number=19

	python plot_cleaned_images.py DBSCAN $(GAMMA_FILE1) $(GAMMA_FILE2) $(G_RUN1)$(G_RUN2) $(FEAT) --number=120

	python plot_cleaned_images.py DBSCAN $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT)

thresholds:
	python plot_cleaned_images.py thresholds $(PROTON_FILE1) $(PROTON_FILE2) $(P_RUN1)$(P_RUN2) $(FEAT) --number=19
	python plot_cleaned_images.py thresholds $(GAMMA_FILE1) $(GAMMA_FILE2) $(G_RUN1)$(G_RUN2) $(FEAT)
	python plot_cleaned_images.py thresholds $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT)

facttools: $(DL1)/proton/$(P_RUN1)$(P_RUN2)_dl1.fits $(DL1)/gamma/$(G_RUN1)$(G_RUN2)_dl1.fits $(DL1)/data/201311$(C_RUN1)_$(C_RUN2)_dl1.fits
	python plot_cleaned_images.py facttools proton $(P_RUN1)$(P_RUN2) $(FEAT) --number=12
	python plot_cleaned_images.py facttools gamma $(G_RUN1)$(G_RUN2) $(FEAT)
	python plot_cleaned_images.py facttools data 201311$(C_RUN1)_$(C_RUN2) $(FEAT)

data:
#	python plot_cleaned_images.py facttools data 201311$(C_RUN1)_$(C_RUN2) $(FEAT)
#	python plot_cleaned_images.py thresholds $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT)
	# TEXINPUTS=../header_matplotlib.tex: \
	# MATPLOTLIBRC=../matplotlibrc_half
	python plot_cleaned_images.py DBSCAN $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT)

$(DL1)/proton/$(P_RUN1)$(P_RUN2)_dl1.fits:
	java -jar ~/Packages/fact-tools/target/fact-tools-v1.0.3-31-g7c675d6.jar ./save_dl1_mc.xml -Dinfile=file:/net/big-tank/POOL/projects/fact/data/open_crab_sample/raw/simulations/proton/00$(P_RUN1)$(P_RUN2).$(P_RUN2)_D_MonteCarlo$(P_RUN1)_Events.fits.gz -Ddrsfile=file:/net/big-tank/POOL/projects/fact/data/open_crab_sample/raw/simulations/mc_drs_file.drs.fits.gz -Doutfile=file:$(DL1)/proton/$(P_RUN1)$(P_RUN2)_dl1.fits

$(DL1)/gamma/$(G_RUN1)$(G_RUN2)_dl1.fits:
	java -jar ~/Packages/fact-tools/target/fact-tools-v1.0.3-31-g7c675d6.jar ./save_dl1_mc.xml -Dinfile=file:/net/big-tank/POOL/projects/fact/data/open_crab_sample/raw/simulations/gamma_diffuse/00$(G_RUN1)$(G_RUN2).$(G_RUN2)_D_MonteCarlo$(G_RUN1)_Events.fits.gz -Ddrsfile=file:/net/big-tank/POOL/projects/fact/data/open_crab_sample/raw/simulations/mc_drs_file.drs.fits.gz -Doutfile=file:$(DL1)/gamma/$(G_RUN1)$(G_RUN2)_dl1.fits

$(DL1)/data/201311$(C_RUN1)_$(C_RUN2)_dl1.fits:
	java -jar ~/Packages/fact-tools/target/fact-tools-v1.0.3-31-g7c675d6.jar ./save_dl1.xml -Dinfile=file:/net/big-tank/POOL/projects/fact/data/open_crab_sample/raw/2013/11/$(C_RUN1)/201311$(C_RUN1)_$(C_RUN2).fits.fz -Ddrsfile=file:/net/big-tank/POOL/projects/fact/data/open_crab_sample/raw/2013/11/$(C_RUN1)/201311$(C_RUN1)_152.drs.fits.gz -Doutfile=file:$(DL1)/data/201311$(C_RUN1)_$(C_RUN2)_dl1.fits


.PHONY: all dbscan thresholds facttools
