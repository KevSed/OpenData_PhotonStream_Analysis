all: eps_*.pdf

eps_*.pdf: plot_eps_scores.py
	TEXINPUTS=../header-matplotlib.tex: MATPLOTLIBRC=../matplotlibrc_full \
	python plot_eps_scores.py
