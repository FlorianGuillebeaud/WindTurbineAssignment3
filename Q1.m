% Assignment 3 
% Question 1
close all
clear all
clc

%% Read Blade and airfoil Data %%
blade_data = xlsread('Blade_data');

mode_shapes = importdata('modeshapes.txt');
W3_100 = importdata('cylinder_ds.txt'); %100% CILINDER
W3_60  = importdata('FFA-W3-600_ds.txt'); %600
W3_48  = importdata('FFA-W3-480_ds.txt'); %480
W3_36  = importdata('FFA-W3-360_ds.txt'); %360
W3_30  = importdata('FFA-W3-301_ds.txt'); %301
W3_24  = importdata('FFA-W3-241_ds.txt'); %241

%%