%%
clc;
close all;
clear all;
set(0,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',20);
set(0, 'DefaultLineLineWidth', 1.5);
addpath(genpath('../utils'));
addpath(genpath('../GGP'));

%%
load('../data/ANC/ANC.mat');
cnts = wordcounts';
clear words wordcounts;

%%
figure;
plot_loglog(cnts);
title('observed data - proportion');
figure;
plot_rank(cnts);
title('observed data - rank');

%%
n_samples = 50000;
results = Model3fit(cnts, n_samples);