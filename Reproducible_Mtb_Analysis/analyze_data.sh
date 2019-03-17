#!/usr/bin/bash

#Download Raw Data
Reads_from_SRA/download_raw_reads.sh

#Correct Raw Data
Reads_from_SRA/correct_reads.sh

#Process Raw Data
./process_raw_reads

#Identify LPI (Low Permissibility to Insertion) sites in avium genome
./find_lpi_sites.py

#Calculate results using TIFA
./tifa_mtb.py
