# INPUTFILES
CRAB_FILE=crab
GAMMA_FILE=gamma/gustav
PROTON_FILE=proton/uwe

P_RUN1=066
P_RUN2=009
G_RUN1=010
G_RUN2=910
C_RUN1=04
C_RUN2=162

DL1=/net/big-tank/POOL/projects/fact/photon-stream/facttools


all: dbscan thresholds facttools

dbscan:
#	python plot_cleaned_images.py DBSCAN $(PROTON_FILE) $(P_RUN1)$(P_RUN2) $(FEAT) --number=19
#	python delta_delta_same_pixel.py DBSCAN $(CRAB_FILE) $(G_RUN1)$(G_RUN2) $(FEAT) --number=120
	TEXINPUTS=../header_matplotlib.tex: \
	MATPLOTLIBRC=../matplotlibrc_half \
	python delta_delta_same_pixel.py DBSCAN $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT) --number=6690

thresholds:
#	python plot_cleaned_images.py thresholds $(PROTON_FILE) $(P_RUN1)$(P_RUN2) $(FEAT) --number=19
#	python plot_cleaned_images.py thresholds $(GAMMA_FILE) $(G_RUN1)$(G_RUN2) $(FEAT)
	TEXINPUTS=../header_matplotlib.tex: \
	MATPLOTLIBRC=../matplotlibrc_half \
	python delta_delta_same_pixel.py thresholds $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT) --number=6690

data:
#	python plot_cleaned_images.py facttools data 201311$(C_RUN1)_$(C_RUN2) $(FEAT)
	TEXINPUTS=../header_matplotlib.tex: \
	MATPLOTLIBRC=../matplotlibrc_half \
	python plot_cleaned_images.py thresholds $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT)
#	python plot_cleaned_images.py DBSCAN $(CRAB_FILE) 201311$(C_RUN1)_$(C_RUN2) $(FEAT)

.PHONY: all dbscan thresholds
