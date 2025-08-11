#!/bin/bash

febench --config ./config.yaml --cwd 7net-0 --calc_type sevennet --model 7net-0.pth --model_path $POTENTIAL_DIRNAME 

febench --config ./config.yaml --cwd esen_oam --calc_type esen --model esen_oam --model_path $POTENTIAL_DIRNAME --modal oam


febench --config ./config.yaml --cwd dpa_openlam_omat --calc_type dpa --model dpa_openlam --model_path $POTENTIAL_DIRNAME --modal omat


febench --config ./config.yaml --cwd orb_omat --calc_type orb --model orb_omat --model_path $POTENTIAL_DIRNAMEs --modal omat

febench --config ./config.yaml --cwd orb/mpa --calc_type orb --model orb_mpa --model_path $POTENTIAL_DIRNAMEs --modal mpa
