# INPUTFILES
CRAB_FILE=/net/big-tank/POOL/projects/fact/photon-stream/stream_data/crab
GAMMA_FILE=/net/big-tank/POOL/projects/fact/photon-stream/stream_data/gamma
PROTON_FILE=/net/big-tank/POOL/projects/fact/photon-stream/stream_data/proton

# OUTPUT
OUTDIR=/home/ksedlaczek/OUT_$(EPS)
DATAOUT=/net/big-tank/POOL/projects/fact/photon-stream/features/$(EPS)

# CONFIG FILES
CONFIG=configs/analysis_config.yaml
CUT_CONFIG=configs/quality_cuts.yaml
DATAMC_CONFIG=configs/data_mc_separator.yaml

# THETA2 CONFIGURATION
PREDICTION_THRESHOLD=0.778
THETA2_CUT=0.08

# EPS
# EPS=0.15

all: $(addprefix $(OUTDIR)/, \
	pdf/theta2_plot.pdf \
	) \
	$(addprefix $(OUTDIR)/, \
	crab_application_done \
	) \
	$(addprefix $(OUTDIR)/, \
	simulation_application_done\
	)

data: proton gamma crab

proton: $(DATAOUT)
	OPENBLAS_MAIN_FREE=1 python3 gendata.py \
		$(DATAOUT)/proton_sample.hdf5 \
		$(PROTON_FILE)/**/*.phs.jsonl.gz \
		-e $(EPS)

gamma: $(DATAOUT)
	OPENBLAS_MAIN_FREE=1 python3 gendata.py \
		$(DATAOUT)/gamma_sample.hdf5 \
		$(GAMMA_FILE)/**/*.phs.jsonl.gz \
		-e $(EPS)

crab: $(DATAOUT)
	OPENBLAS_MAIN_FREE=1 python3 gendata.py \
		$(DATAOUT)/crab_data.hdf5 \
		$(CRAB_FILE)/*.phs.jsonl.gz \
		-e $(EPS)


# Apply precuts to the files
$(OUTDIR)/crab_data_precuts.hdf5: $(DATAOUT)/crab_data.hdf5 $(CUT_CONFIG) | $(OUTDIR)
	aict_apply_cuts ./configs/quality_cuts.yaml \
		$(DATAOUT)/crab_data.hdf5 \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		--chunksize=10000

$(OUTDIR)/gamma_precuts.hdf5: $(DATAOUT)/gamma_sample.hdf5 configs/quality_cuts.yaml | $(OUTDIR)
	aict_apply_cuts ./configs/quality_cuts.yaml \
		$(DATAOUT)/gamma_sample.hdf5 \
		$(OUTDIR)/gamma_precuts.hdf5 \
		--chunksize=10000

$(OUTDIR)/proton_precuts.hdf5: $(DATAOUT)/proton_sample.hdf5 $(CUT_CONFIG) | $(OUTDIR)
	aict_apply_cuts \
		$(CUT_CONFIG)\
		$(DATAOUT)/proton_sample.hdf5 \
		$(OUTDIR)/proton_precuts.hdf5 \
		--chunksize=10000


# Train Separator
$(OUTDIR)/separator.pkl $(OUTDIR)/separator_performance.hdf5: $(CONFIG) $(OUTDIR)/proton_precuts.hdf5 $(OUTDIR)/gamma_precuts.hdf5 | $(OUTDIR)
	aict_train_separation_model \
		$(CONFIG) \
		$(OUTDIR)/gamma_precuts.hdf5 \
		$(OUTDIR)/proton_precuts.hdf5 \
		$(OUTDIR)/separator_performance.hdf5 \
		$(OUTDIR)/separator.pkl

# Train disp model
$(OUTDIR)/disp_model.pkl $(OUTDIR)/sign_model.pkl $(OUTDIR)/cv_disp.hdf5: $(CONFIG) $(OUTDIR)/gamma_precuts.hdf5 | $(OUTDIR)
	aict_train_disp_regressor \
		$(CONFIG) \
		$(OUTDIR)/gamma_precuts.hdf5 \
		$(OUTDIR)/cv_disp.hdf5 \
		$(OUTDIR)/disp_model.pkl \
		$(OUTDIR)/sign_model.pkl


# Train energy regression
$(OUTDIR)/energy_model.pkl $(OUTDIR)/energy_performance.hdf5: $(CONFIG) $(OUTDIR)/gamma_precuts.hdf5 | $(OUTDIR)
	aict_train_energy_regressor \
		$(CONFIG) \
		$(OUTDIR)/gamma_precuts.hdf5 \
		$(OUTDIR)/energy_performance.hdf5 \
		$(OUTDIR)/energy_model.pkl


# Apply models to simulations
$(OUTDIR)/simulation_application_done: $(CONFIG) $(OUTDIR)/separator.pkl $(CONFIG) $(OUTDIR)/disp_model.pkl $(OUTDIR)/sign_model.pkl $(CONFIG) $(OUTDIR)/energy_model.pkl | $(OUTDIR)/gamma_precuts.hdf5 $(OUTDIR)/proton_precuts.hdf5
	aict_apply_separation_model \
		$(CONFIG) \
		$(OUTDIR)/gamma_precuts.hdf5 \
		$(OUTDIR)/separator.pkl \
		--chunksize=100000 --yes

	aict_apply_separation_model \
		$(CONFIG) \
		$(OUTDIR)/proton_precuts.hdf5 \
		$(OUTDIR)/separator.pkl \
		--chunksize=100000 --yes

	aict_apply_disp_regressor \
		$(CONFIG) \
		$(OUTDIR)/gamma_precuts.hdf5 \
		$(OUTDIR)/disp_model.pkl \
		$(OUTDIR)/sign_model.pkl \
		--yes --chunksize=100000

	aict_apply_disp_regressor \
		$(CONFIG) \
		$(OUTDIR)/proton_precuts.hdf5 \
		$(OUTDIR)/disp_model.pkl \
		$(OUTDIR)/sign_model.pkl \
		--yes --chunksize=100000

	fact_calculate_theta $(OUTDIR)/proton_precuts.hdf5 --yes
	fact_calculate_theta $(OUTDIR)/gamma_precuts.hdf5 --yes

	touch $(OUTDIR)/simulation_application_done

# Apply models to crab data
$(OUTDIR)/crab_application_done: $(CONFIG) $(OUTDIR)/separator.pkl $(CONFIG) $(OUTDIR)/disp_model.pkl $(OUTDIR)/sign_model.pkl $(CONFIG) $(OUTDIR)/energy_model.pkl | $(OUTDIR)/crab_data_precuts.hdf5
	aict_apply_separation_model \
		$(CONFIG) \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		$(OUTDIR)/separator.pkl \
		--chunksize=100000 --yes

	aict_apply_disp_regressor \
		$(CONFIG) \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		$(OUTDIR)/disp_model.pkl \
		$(OUTDIR)/sign_model.pkl \
		--yes --chunksize=100000

	aict_apply_energy_regressor \
		$(CONFIG) \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		$(OUTDIR)/energy_model.pkl \
		--chunksize=100000 --yes

	fact_calculate_theta $(OUTDIR)/crab_data_precuts.hdf5 --yes --source CRAB
	fact_calculate_radec $(OUTDIR)/crab_data_precuts.hdf5 --yes

	touch $(OUTDIR)/crab_application_done

# Generate plots
$(OUTDIR)/pdf/theta2_plot.pdf: $(OUTDIR)/pdf $(OUTDIR)/crab_data_precuts.hdf5 | $(OUTDIR)/crab_application_done
	MATPLOTLIBRC=matplotlibrc_half
	TEXINPUTS=$$(pwd): fact_plot_theta_squared \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		--threshold=$(PREDICTION_THRESHOLD) \
		--theta2-cut=$(THETA2_CUT) \
		-o $(OUTDIR)/pdf/theta2_plot.pdf


theta2:
	MATPLOTLIBRC=matplotlibrc_full \
	TEXINPUTS=$$(pwd): fact_plot_theta_squared \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		--threshold=$(PREDICTION_THRESHOLD) \
		--theta2-cut=$(THETA2_CUT) \
		-o $(OUTDIR)/pdf/theta2_plot.pdf

skymap:
	MATPLOTLIBRC=matplotlibrc_half
	TEXINPUTS=$$(pwd): fact_plot_skymap \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		--threshold=$(PREDICTION_THRESHOLD) \
		--key events \
		-o $(OUTDIR)/pdf/data_skymap.pdf

protonskymap:
	# fact_calculate_radec $(OUTDIR)/proton_precuts.hdf5 --yes
	python skymaps.py \
		$(OUTDIR)/proton_precuts.hdf5 \
		$(OUTDIR)/pdf/proton_skymap.pdf \
		proton \
		--threshold=0.4 \

gammaskymap:
	# fact_calculate_radec $(OUTDIR)/proton_precuts.hdf5 --yes
	python skymaps.py \
		$(OUTDIR)/gamma_precuts.hdf5 \
		$(OUTDIR)/pdf/gamma_skymap.pdf \
		gamma \
		--threshold=0.0

theta:
	MATPLOTLIBRC=matplotlibrc_half \
	TEXINPUTS=$$(pwd): fact_plot_theta_squared \
	$(OUTDIR)/gamma_precuts.hdf5 \
		--theta2-cut=$(THETA2_CUT) \
		-o $(OUTDIR)/pdf/theta2_plot_g.pdf

	MATPLOTLIBRC=matplotlibrc_half \
	TEXINPUTS=$$(pwd): fact_plot_theta_squared \
	$(OUTDIR)/proton_precuts.hdf5 \
		--threshold=0 \
		--theta2-cut=$(THETA2_CUT) \
		-o $(OUTDIR)/pdf/theta2_plot_p.pdf


mcdata_p: configs/data_mc.yaml | $(OUTDIR)/pdf
	TEXINPUTS=$$(pwd)/header-matplotlib.tex: MATPLOTLIBRC=./matplotlibrc_half \
	python plot_data_mc_comparison.py \
		$(OUTDIR)/crab_data_precuts.hdf5 \
		$(OUTDIR)/proton_precuts.hdf5 \
		$(OUTDIR)/pdf/proton_data_mc.pdf \
		-c configs/data_mc.yaml \
		-f 1

std_phs: $(OUTDIR)/pdf
	TEXINPUTS=$$(pwd)/header-matplotlib.tex: MATPLOTLIBRC=./matplotlibrc_half \
	python phs_std_comparison.py \
		$(EPS) \
		crab
#
#mcdata2: configs/data_mc.yaml | $(OUTDIR)/pdf
#	MATPLOTLIBBACKEND=agg python3 dmc_com.py \
#		$(OUTDIR)
#
mcdata: configs/data_mc.yaml | $(OUTDIR)/pdf
	TEXINPUTS=$$(pwd)/header-matplotlib.tex: \
	MATPLOTLIBRC=./matplotlibrc_full \
	python ks_test.py \
		--cuts \
		--threshold=0.0 \
		$(OUTDIR)

#	bash rm_models.sh
#
#	aict_train_separation_model \
#		$(DATAMC_CONFIG) \
#		$(OUTDIR)/crab_data_precuts.hdf5 \
#		$(OUTDIR)/proton_precuts.hdf5 \
#		$(OUTDIR)/data_mc_separator_performance.hdf5 \
#		$(OUTDIR)/data_mc_separator.pkl
#
#	MATPLOTLIBRC=matplotlibrc_half TEXINPUTS=$$(pwd):
#	aict_plot_separator_performance \
#	       $(DATAMC_CONFIG) \
#	       $(OUTDIR)/data_mc_separator_performance.hdf5 \
#	       $(OUTDIR)/data_mc_separator.pkl \
#	       -o $(OUTDIR)/pdf/data_mc_separation.pdf


# Make performance plots of energy regression and separation
plot_performance: # $(CONFIG) $(OUTDIR)/separator_performance.hdf5 $(OUTDIR)/separator.pkl $(CONFIG) $(OUTDIR)/energy_performance.hdf5	$(OUTDIR)/energy_model.pkl | $(OUTDIR) $(OUTDIR)/pdf
#	TEXINPUTS=$$(pwd)/header-matplotlib.tex: \
#	MATPLOTLIBRC=./matplotlibrc_full
	aict_plot_regressor_performance \
	       $(CONFIG) \
	       $(OUTDIR)/energy_performance.hdf5 \
	       $(OUTDIR)/energy_model.pkl \
	       -o $(OUTDIR)/pdf/energy_performance.pdf


#	TEXINPUTS=$$(pwd)/header-matplotlib.tex: \
#	MATPLOTLIBRC=./matplotlibrc_full
	aict_plot_separator_performance \
	       $(CONFIG) \
	       $(OUTDIR)/separator_performance.hdf5 \
	       $(OUTDIR)/separator.pkl \
	       -o $(OUTDIR)/pdf/separation_performance.pdf

meta:
	python add_meta_data.py $(DATAOUT)/crab_data.hdf5

rate:
	TEXINPUTS=$$(pwd)/header-matplotlib.tex: \
	MATPLOTLIBRC=./matplotlibrc_full \
	python plot_ratescan.py

delta:
	TEXINPUTS=$$(pwd)/header-matplotlib.tex: \
	MATPLOTLIBRC=./matplotlibrc_half \
	python delta_delta/plot_delta_delta.py $(OUTDIR)/crab_data_precuts.hdf5 -o $(OUTDIR)/pdf/delta_delta.pdf -t $(PREDICTION_THRESHOLD)

sdelta:
#	TEXINPUTS=$$(pwd)/header-matplotlib.tex: \
#	MATPLOTLIBRC=./matplotlibrc_half
	python delta_delta/plot_specific_delta_delta.py $(OUTDIR)/crab_data_precuts.hdf5 -o $(OUTDIR)/pdf/delta_delta_specific -t $(PREDICTION_THRESHOLD)


$(OUTDIR):
	mkdir -p $(OUTDIR)

$(DATAOUT):
	mkdir -p $(DATAOUT)

$(OUTDIR)/pdf: $(OUTDIR)
	mkdir -p $(OUTDIR)/pdf

cuts:
	python3 detection_gridsearch.py \
		$(OUTDIR)/crab_data_precuts.hdf5

clean:
	rm -rf $(OUTDIR)


.PHONY: all clean cuts plot_performance data
