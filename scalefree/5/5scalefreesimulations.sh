#!/bin/bash


./macgregor_20_scalefree < parameters.txt > output_20_5.txt
mv spikes.txt spikes_20_5.txt

./macgregor_25_scalefree < parameters.txt > output_25_5.txt
mv spikes.txt spikes_25_5.txt

./macgregor_30_scalefree < parameters.txt > output_30_5.txt
mv spikes.txt spikes_30_5.txt

./macgregor_35_scalefree < parameters.txt > output_35_5.txt
mv spikes.txt spikes_35_5.txt

./macgregor_40_scalefree < parameters.txt > output_40_5.txt
mv spikes.txt spikes_40_5.txt