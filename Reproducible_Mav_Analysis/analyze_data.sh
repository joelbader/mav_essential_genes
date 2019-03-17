#!/usr/bin/bash

mkdir -p output
mkdir -p intermed_output

#Process Raw Data
./process_raw_reads

#Identify LPI (Low Permissibility to Insertion) sites in avium genome
./find_lpi_sites.py

#Calculate results using TIFA
./tifa.py
